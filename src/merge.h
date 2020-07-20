#ifndef MERGE_H
#define MERGE_H

#include <limits>

#include <boost/icl/split_interval_map.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/progress.hpp>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "version.h"
#include "scan.h"
#include "util.h"


namespace coralns
{

  #ifndef NAN_COUNT
  #define NAN_COUNT -99999999
  #endif

  struct MergeConfig {
    typedef std::map<std::string, uint32_t> TChrMap;
    typedef std::vector<std::string> TRefMap;
    typedef std::pair<uint32_t, uint32_t> TInterval;
    typedef std::map<TInterval, uint32_t> TChrIntervals;
    typedef std::vector<TChrIntervals> TGenomicIntervals;
    typedef std::vector<std::string> TSampleNames;
    
    TChrMap chrMap;
    TRefMap refName;
    TGenomicIntervals gIntervals;
    std::string byType;
    TSampleNames sampleNames;
    boost::filesystem::path outfile;
    std::vector<boost::filesystem::path> countFiles;
  };

  inline void
  collectCounts(MergeConfig const& c) {
    // Open output file
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));

    // Header
    dataOut << "chr\tstart\tend";
    for(uint32_t i = 0; i < c.countFiles.size(); ++i) dataOut << "\t" << c.sampleNames[i];
    dataOut << std::endl;
    
    // Fetch all intervals and sample names
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Collect counts" << std::endl;
    boost::progress_display show_progress( c.gIntervals.size() );
    for(uint32_t refIndex = 0; refIndex < c.gIntervals.size(); ++refIndex) {
      ++show_progress;
      typedef std::vector<float> TChrCounts;
      typedef std::vector<TChrCounts> TFileCounts;
      TFileCounts counts(c.countFiles.size(), TChrCounts());
      for(uint32_t i = 0; i < c.countFiles.size(); ++i) {
	counts[i].resize(c.gIntervals[refIndex].size(), NAN_COUNT);
	std::ifstream file(c.countFiles[i].string().c_str(), std::ios_base::in | std::ios_base::binary);
	boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
	dataIn.push(boost::iostreams::gzip_decompressor());
	dataIn.push(file);
	std::istream instream(&dataIn);
	std::string line;
	while(std::getline(instream, line)) {
	  if (line.size()) {
	    typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	    boost::char_separator<char> sep("\t");
	    Tokenizer tokens(line, sep);
	    Tokenizer::iterator tokIter = tokens.begin();
	    if (tokIter != tokens.end()) {
	      std::string chrom = *tokIter++;
	      if (tokIter != tokens.end()) {
		if (std::string(*tokIter) == "start") {
		  // Header
		  continue;
		}
		int32_t start = boost::lexical_cast<int32_t>(*tokIter++);
		int32_t end = boost::lexical_cast<int32_t>(*tokIter++);
		++tokIter; // mappable positions
		float countval = boost::lexical_cast<float>(*tokIter++);
		float cn = boost::lexical_cast<float>(*tokIter++);
		typename MergeConfig::TChrMap::const_iterator cit = c.chrMap.find(chrom);
		uint32_t ridx = cit->second;
		if (ridx != refIndex) continue;
		typename MergeConfig::TChrIntervals::const_iterator it = c.gIntervals[refIndex].find(std::make_pair(start, end));
		if (c.byType == "CN") counts[i][it->second] = cn;
		else counts[i][it->second] = countval;
	      }
	    }
	  }
	}
	dataIn.pop();
      }
      // Output
      for(MergeConfig::TChrIntervals::const_iterator it = c.gIntervals[refIndex].begin(); it != c.gIntervals[refIndex].end(); ++it) {
	dataOut << c.refName[refIndex] << '\t' << it->first.first << '\t' << it->first.second;
	for(uint32_t i = 0; i < c.countFiles.size(); ++i) {
	  if (counts[i][it->second] == NAN_COUNT) dataOut << "\tNaN";
	  else dataOut << '\t' << counts[i][it->second];
	}
	dataOut << std::endl;
      }
    }
    
    // Close
    dataOut.pop();
  }

  int merge(int argc, char **argv) {
    MergeConfig c;
    
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("byType,b", boost::program_options::value<std::string>(&c.byType)->default_value("CN"), "merge by 'CN' or 'counts'")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("cov.gz"), "coverage output file")
      ;

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<std::vector<boost::filesystem::path> >(&c.countFiles), "input count files")
      ;

    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);

    // Set the visibility
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic);

    // Parse command-line
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);

    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file"))) {
      std::cout << std::endl;
      std::cout << "Usage: coral " << argv[0] << " [OPTIONS] <sample1.gz> <sample2.gz> ..." << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    }

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "coral ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;

    // Check input count files
    c.sampleNames.resize(c.countFiles.size());
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Parse input intervals" << std::endl;
    boost::progress_display show_progress( c.countFiles.size() );
    for(uint32_t i = 0; i < c.countFiles.size(); ++i) {
      ++show_progress;
      if (!(boost::filesystem::exists(c.countFiles[i]) && boost::filesystem::is_regular_file(c.countFiles[i]) && boost::filesystem::file_size(c.countFiles[i]))) {
	std::cerr << "Input coverage file is missing: " << c.countFiles[i].string() << std::endl;
	return 1;
      } else {
	// Fetch all intervals
	std::ifstream file(c.countFiles[i].string().c_str(), std::ios_base::in | std::ios_base::binary);
	boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
	dataIn.push(boost::iostreams::gzip_decompressor());
	dataIn.push(file);
	std::istream instream(&dataIn);
	std::string line;
	while(std::getline(instream, line)) {
	  if (line.size()) {
	    typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	    boost::char_separator<char> sep("\t");
	    Tokenizer tokens(line, sep);
	    Tokenizer::iterator tokIter = tokens.begin();
	     if (tokIter != tokens.end()) {
	       std::string chrom = *tokIter++;
	       if (tokIter != tokens.end()) {
		 if (std::string(*tokIter) == "start") {
		   ++tokIter; ++tokIter;
		   std::string name = *tokIter;
		   std::size_t lastidx = name.find_last_of("_");
		   c.sampleNames[i] = name.substr(0, lastidx);
		   // Header
		   continue;
		 }
		 int32_t start = boost::lexical_cast<int32_t>(*tokIter++);
		 int32_t end = boost::lexical_cast<int32_t>(*tokIter++);
		 if (start < end) {
		   if (c.chrMap.find(chrom) == c.chrMap.end()) {
		     c.chrMap.insert(std::make_pair(chrom, c.chrMap.size()));
		     c.refName.resize(c.chrMap.size(), chrom);
		   }
		   uint32_t refIndex = c.chrMap[chrom];
		   if (refIndex >= c.gIntervals.size()) c.gIntervals.resize(refIndex + 1, MergeConfig::TChrIntervals());
		   c.gIntervals[refIndex].insert(std::make_pair(std::make_pair(start, end), 0));
		 }
	       }
	     }
	  }
	}
	dataIn.pop();
      }
    }

    // Assign positions
    for(uint32_t refIndex = 0; refIndex < c.gIntervals.size(); ++refIndex) {
      uint32_t idx = 0;
      for(MergeConfig::TChrIntervals::iterator it = c.gIntervals[refIndex].begin(); it != c.gIntervals[refIndex].end(); ++it, ++idx) it->second = idx;
    }

    // Get all values
    collectCounts(c);

    // End
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
    return 0;
  }

  
}

#endif
