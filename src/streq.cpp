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
  unsigned short minMapQual;
  uint32_t segment;
  boost::filesystem::path watsonRatio;
  boost::filesystem::path ww;
  boost::filesystem::path wc;
  boost::filesystem::path region;
};


template<typename TConfig, typename TRefIndex, typename TIntervalIter, typename TOutfile>
inline void
_countSegments(bam_hdr_t* hdr, samFile* wwfile, hts_idx_t* wwidx, TConfig const& c, TRefIndex const refIndex, TIntervalIter const& itR, std::string const& id, TOutfile& ofile) {
  std::string intervalName = itR->second;
  uint32_t intervalSize = itR->first.second - itR->first.first;
  uint32_t bins = intervalSize / c.segment + 1;

  typedef std::vector<uint32_t> TCounter;
  TCounter watsonCount;
  TCounter crickCount;
  watsonCount.resize(bins, 0);
  crickCount.resize(bins, 0);

  hts_itr_t* iter = sam_itr_queryi(wwidx, refIndex, itR->first.first, itR->first.second);
  bam1_t* rec = bam_init1();
  while (sam_itr_next(wwfile, iter, rec) >= 0) {
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
  bam_destroy1(rec);
  hts_itr_destroy(iter);
	
  typename TCounter::const_iterator itWatson = watsonCount.begin();
  typename TCounter::const_iterator itCrick =crickCount.begin();
  for(std::size_t bin = 0; bin < bins; ++itWatson, ++itCrick, ++bin) {
    uint32_t sup = *itWatson + *itCrick;
    if (sup > 0) {
      double wRatio = ((double) *itWatson / (double) (sup));
      int32_t regionStart = bin * c.segment + itR->first.first;
      int32_t regionEnd = std::min((uint32_t) ((bin + 1) * c.segment + itR->first.first), itR->first.second);
      ofile << hdr->target_name[refIndex] << '\t' << regionStart << '\t' << regionEnd << '\t' << wRatio << '\t' << sup << "\t" << id << "\t" << intervalName << std::endl;
    }
  }
}
 

int main(int argc, char **argv) {
  Config c;

  // Parameter
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("map-qual,q", boost::program_options::value<unsigned short>(&c.minMapQual)->default_value(1), "min. mapping quality")
    ("segment,e", boost::program_options::value<uint32_t>(&c.segment)->default_value(5000), "segment size")
    ("samestrand,s", boost::program_options::value<boost::filesystem::path>(&c.ww)->default_value("ww.bam"), "input same strand bam")
    ("diffstrand,d", boost::program_options::value<boost::filesystem::path>(&c.wc)->default_value("wc.bam"), "input different strand bam")
    ("wcratio,a", boost::program_options::value<boost::filesystem::path>(&c.watsonRatio)->default_value("watson.out"), "output file for WC ratio")
    ;


  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&c.region), "input bed file")
    ;

  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);


  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("input-file"))) {
    std::cout << "Usage: " << argv[0] << " [OPTIONS] <regions.bed>" << std::endl;
    std::cout << visible_options << "\n";
    return 1;
  }

  // Check ww file
  if (!(boost::filesystem::exists(c.ww) && boost::filesystem::is_regular_file(c.ww) && boost::filesystem::file_size(c.ww))) {
    std::cerr << "Same strand BAM file is missing: " << c.ww.string() << std::endl;
    return 1;
  }
  
  // Check wc file
  if (!(boost::filesystem::exists(c.wc) && boost::filesystem::is_regular_file(c.wc) && boost::filesystem::file_size(c.wc))) {
    std::cerr << "Different strand BAM file is missing: " << c.wc.string() << std::endl;
    return 1;
  }

  // Check regions file
  if (!(boost::filesystem::exists(c.region) && boost::filesystem::is_regular_file(c.region) && boost::filesystem::file_size(c.region))) {
    std::cerr << "Region file is missing: " << c.region.string() << std::endl;
    return 1;
    }

  // Load bam file
  samFile* wwfile = sam_open(c.ww.string().c_str(), "r");
  hts_idx_t* wwidx = sam_index_load(wwfile, c.ww.string().c_str());
  samFile* wcfile = sam_open(c.wc.string().c_str(), "r");
  hts_idx_t* wcidx = sam_index_load(wcfile, c.wc.string().c_str());
  bam_hdr_t* hdr = sam_hdr_read(wwfile);

  // Store intervals
  typedef std::pair<uint32_t, uint32_t> TInterval;
  typedef std::map<TInterval, std::string> TIntervalMap;
  typedef std::vector<TIntervalMap> TGenomicIntervals;
  TGenomicIntervals genomicIntervals;
  genomicIntervals.resize(hdr->n_targets);
  std::ifstream regionFile(c.region.string().c_str(), std::ifstream::in);
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

  // Output WC ratios
  std::ofstream ofile(c.watsonRatio.string().c_str());
  ofile << "chr\tstart\tend\twratio\tsupport\ttype\tid" << std::endl;
  for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
    for(TIntervalMap::const_iterator itR = genomicIntervals[refIndex].begin(); itR != genomicIntervals[refIndex].end(); ++itR) {
      _countSegments(hdr, wwfile, wwidx, c, refIndex, itR, "Watson", ofile);
      _countSegments(hdr, wcfile, wcidx, c, refIndex, itR, "WatsonCrick", ofile);
    }
  }
  ofile.close();
  

  // Close bam
  bam_hdr_destroy(hdr);
  hts_idx_destroy(wcidx);
  hts_idx_destroy(wwidx);
  sam_close(wwfile);
  sam_close(wcfile);

  return 0;
}
