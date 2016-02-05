/*
============================================================================
Strand-Seq Watson-Crick Region Analyzer
============================================================================
Copyright (C) 2016 Tobias Rausch

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
============================================================================
Contact: Tobias Rausch (rausch@embl.de)
============================================================================
*/

#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <fstream>

#define BOOST_DISABLE_ASSERTS
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <htslib/sam.h>

#include "util.h"

using namespace streq;

struct Config {
  bool hasRegionFile;
  unsigned short minMapQual;
  uint32_t window;
  uint32_t segment;
  boost::filesystem::path watsonRatio;
  boost::filesystem::path regions;
  boost::filesystem::path loadPre;
  std::vector<boost::filesystem::path> files;
};

int main(int argc, char **argv) {
  Config c;

  // Parameter
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("map-qual,q", boost::program_options::value<unsigned short>(&c.minMapQual)->default_value(1), "min. mapping quality")
    ("window,w", boost::program_options::value<uint32_t>(&c.window)->default_value(1000000), "window length")
    ("segment,s", boost::program_options::value<uint32_t>(&c.segment)->default_value(5000), "segment size")
    ("loadpre,l", boost::program_options::value<boost::filesystem::path>(&c.loadPre)->default_value("pre.out"), "preprocessing information")
    ("wcratio,a", boost::program_options::value<boost::filesystem::path>(&c.watsonRatio)->default_value("watson.out"), "output file for WC ratio")
    ;

  boost::program_options::options_description region("Region options");
  region.add_options()
    ("regions,r", boost::program_options::value<boost::filesystem::path>(&c.regions), "bed file with regions to analyze")
    ;

  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input bam file")
    ;

  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(region).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic).add(region);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);


  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("loadpre"))) {
    std::cout << "Usage: " << argv[0] << " [OPTIONS] <strand.seq.bam>" << std::endl;
    std::cout << visible_options << "\n";
    return 1;
  }

  // Check preprocessing file
  if (!(boost::filesystem::exists(c.loadPre) && boost::filesystem::is_regular_file(c.loadPre) && boost::filesystem::file_size(c.loadPre))) {
    std::cerr << "Pre-processing information is missing: " << c.loadPre.string() << std::endl;
    return 1;
  }

  // Check regions file
  if (vm.count("regions")) {
    if (!(boost::filesystem::exists(c.regions) && boost::filesystem::is_regular_file(c.regions) && boost::filesystem::file_size(c.regions))) {
      std::cerr << "Region file is missing: " << c.regions.string() << std::endl;
      return 1;
    }
    c.hasRegionFile = true;
  } else c.hasRegionFile = false;

  // Load bam file
  typedef std::vector<samFile*> TSamFile;
  typedef std::vector<hts_idx_t*> TIndex;
  TSamFile samfile;
  TIndex idx;
  samfile.resize(c.files.size());
  idx.resize(c.files.size());
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
    idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
  }
  bam_hdr_t* hdr = sam_hdr_read(samfile[0]);

  // Load pre-processing information
  typedef uint32_t TPos;
  typedef uint32_t TFileIndex;
  typedef std::pair<TPos, TFileIndex> TPosFilePair;
  typedef std::vector<TPosFilePair> TWindowList;
  typedef std::vector<TWindowList> TGenomicWindows;
  TGenomicWindows watsonWindows;
  watsonWindows.resize(hdr->n_targets);
  TGenomicWindows crickWindows;
  crickWindows.resize(hdr->n_targets);
  TGenomicWindows wcWindows;
  wcWindows.resize(hdr->n_targets);
  std::ifstream preFile(c.loadPre.string().c_str(), std::ifstream::in);
  if (preFile.is_open()) {
    while (preFile.good()) {
      std::string line;
      getline(preFile, line);
      typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
      boost::char_separator<char> sep(" \t,;");
      Tokenizer tokens(line, sep);
      Tokenizer::iterator tokIter = tokens.begin();
      if (tokIter!=tokens.end()) {
	int32_t ind = boost::lexical_cast<int32_t>(*tokIter++);
	int32_t refIndex = boost::lexical_cast<int32_t>(*tokIter++);
	int32_t bin = boost::lexical_cast<int32_t>(*tokIter++);
	int32_t fileInd = boost::lexical_cast<int32_t>(*tokIter++);
	if (ind == 0) watsonWindows[refIndex].push_back(std::make_pair(bin, fileInd));
	else if (ind == 1) crickWindows[refIndex].push_back(std::make_pair(bin, fileInd));
	else if (ind == 2) wcWindows[refIndex].push_back(std::make_pair(bin, fileInd));
      }
    }
    preFile.close();
  }

  // Sort windows
  for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
    std::sort(watsonWindows[refIndex].begin(), watsonWindows[refIndex].end());
    std::sort(crickWindows[refIndex].begin(), crickWindows[refIndex].end());
    std::sort(wcWindows[refIndex].begin(), wcWindows[refIndex].end());
  }

  // Store intervals
  typedef std::pair<uint32_t, uint32_t> TInterval;
  typedef std::map<TInterval, std::string> TIntervalMap;
  typedef std::vector<TIntervalMap> TGenomicIntervals;
  TGenomicIntervals genomicIntervals;
  genomicIntervals.resize(hdr->n_targets);
  if (c.hasRegionFile) {
    std::ifstream regionFile(c.regions.string().c_str(), std::ifstream::in);
    if (regionFile.is_open()) {
      while (regionFile.good()) {
	std::string line;
	getline(regionFile, line);
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep(" \t,;");
	Tokenizer tokens(line, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter!=tokens.end()) {
	  std::string chrName = *tokIter++;
	  int32_t tid = bam_name2id(hdr, chrName.c_str());
	  if (tid >= 0) {
	    uint32_t start = boost::lexical_cast<int32_t>(*tokIter++);
	    uint32_t end = boost::lexical_cast<int32_t>(*tokIter++);
	    std::string intId = *tokIter++;
	    genomicIntervals[tid].insert(std::make_pair(std::make_pair(start, end), intId));
	  }
	}
      }
      regionFile.close();
    }
  } else {
    for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
      if (hdr->target_len[refIndex] < c.window) continue;
      std::string intId(hdr->target_name[refIndex]);
      intId += ":" + boost::lexical_cast<std::string>(0) + "-" + boost::lexical_cast<std::string>(hdr->target_len[refIndex]);
      genomicIntervals[refIndex].insert(std::make_pair(std::make_pair(0, hdr->target_len[refIndex]), intId));
    }
  }

  // Output WC ratios
  std::ofstream ofile(c.watsonRatio.string().c_str());
  ofile << "chr\tstart\tend\twratio\tsupport\ttype\tid" << std::endl;
  for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
    for(TIntervalMap::const_iterator itR = genomicIntervals[refIndex].begin(); itR != genomicIntervals[refIndex].end(); ++itR) {
      std::string intervalName = itR->second;
      uint32_t intervalSize = itR->first.second - itR->first.first;
      uint32_t indS = (int) (itR->first.first / c.window);
      uint32_t bins = intervalSize / c.segment + 1;
      for(std::size_t i = 0; i<3; ++i) {
	typedef std::vector<uint32_t> TCounter;
	TCounter watsonCount;
	TCounter crickCount;
	watsonCount.resize(bins, 0);
	crickCount.resize(bins, 0);
	
	TWindowList::const_iterator itWindows;
	TWindowList::const_iterator itWindowsEnd;
	if (i == 0) {
	  itWindows = std::lower_bound(watsonWindows[refIndex].begin(), watsonWindows[refIndex].end(), std::make_pair(indS, (uint32_t) 0));
	  itWindowsEnd = watsonWindows[refIndex].end();
	} else if (i == 1) {
	  itWindows = std::lower_bound(crickWindows[refIndex].begin(), crickWindows[refIndex].end(), std::make_pair(indS, (uint32_t) 0));
	  itWindowsEnd = crickWindows[refIndex].end();
	} else if (i == 2) {
	  itWindows = std::lower_bound(wcWindows[refIndex].begin(), wcWindows[refIndex].end(), std::make_pair(indS, (uint32_t) 0));
	  itWindowsEnd = wcWindows[refIndex].end();
	}
	for(;itWindows != itWindowsEnd; ++itWindows) {
	  if (itWindows->first * c.window >= itR->first.second) break;
	  int32_t regionStart = std::max(itWindows->first * c.window, itR->first.first);
	  int32_t regionEnd = std::min((itWindows->first + 1) * c.window, itR->first.second);
	  hts_itr_t* iter = sam_itr_queryi(idx[itWindows->second], refIndex, regionStart, regionEnd);
	  bam1_t* rec = bam_init1();
	  while (sam_itr_next(samfile[itWindows->second], iter, rec) >= 0) {
	    if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	    if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
	    
	    int32_t pos = rec->core.pos + halfAlignmentLength(rec);
	    if ((pos >= (int32_t) itR->first.first) && (pos < (int32_t) itR->first.second)) {
	      uint32_t binInd = (uint32_t) ((pos - itR->first.first) / c.segment);
	      if (rec->core.flag & BAM_FREAD1) 
		if (rec->core.flag & BAM_FREVERSE) ++crickCount[binInd];
		else ++watsonCount[binInd];
	      else
		if (rec->core.flag & BAM_FREVERSE) ++watsonCount[binInd];
		else ++crickCount[binInd];
	    }
	  }
	}
	
	TCounter::const_iterator itWatson = watsonCount.begin();
	TCounter::const_iterator itCrick =crickCount.begin();
	for(std::size_t bin = 0; bin < bins; ++itWatson, ++itCrick, ++bin) {
	  uint32_t sup = *itWatson + *itCrick;
	  if (sup > 0) {
	    double wRatio = ((double) *itWatson / (double) (sup));
	    int32_t regionStart = bin * c.segment + itR->first.first;
	    int32_t regionEnd = std::min((uint32_t) ((bin + 1) * c.segment + itR->first.first), itR->first.second);
	    ofile << hdr->target_name[refIndex] << '\t' << regionStart << '\t' << regionEnd << '\t' << wRatio << '\t' << sup << '\t';
	    if (i == 0) ofile << "Watson" << '\t';
	    else if (i == 1) ofile << "Crick" << '\t';
	    else if (i == 2) ofile << "WatsonCrick" << '\t';
	    ofile << intervalName << std::endl;
	  }
	}
      }
    }
  }
  ofile.close();

  // Close bam
  bam_hdr_destroy(hdr);
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    hts_idx_destroy(idx[file_c]);
    sam_close(samfile[file_c]);
  }

  return 0;
}
