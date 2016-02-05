/*
============================================================================
Strand-Seq Watson-Crick Classifier
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
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <htslib/sam.h>

#include "util.h"

using namespace streq;

struct Config {
  unsigned short minMapQual;
  uint32_t window;
  boost::filesystem::path ww;
  boost::filesystem::path wc;
  std::vector<boost::filesystem::path> files;
};

int main(int argc, char **argv) {

#ifdef PROFILE
  ProfilerStart("preprocess.prof");
#endif

  Config c;

  // Parameter
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("map-qual,q", boost::program_options::value<unsigned short>(&c.minMapQual)->default_value(1), "min. mapping quality")
    ("window,w", boost::program_options::value<uint32_t>(&c.window)->default_value(1000000), "window length")
    ("samestrand,s", boost::program_options::value<boost::filesystem::path>(&c.ww)->default_value("ww.bam"), "output same strand bam")
    ("diffstrand,d", boost::program_options::value<boost::filesystem::path>(&c.wc)->default_value("wc.bam"), "output different strand bam")
    ;

  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input bam file")
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
    std::cout << "Usage: " << argv[0] << " [OPTIONS] <strand.seq1.bam> <strand.seq2.bam> ... <strand.seqN.bam>" << std::endl;
    std::cout << visible_options << "\n";
    return 1;
  } 

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  // Load bam files
  typedef std::vector<samFile*> TSamFile;
  typedef std::vector<hts_idx_t*> TIndex;
  TSamFile samfile;
  TIndex idx;
  samfile.resize(c.files.size());
  idx.resize(c.files.size());
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
    if (samfile[file_c] == NULL) {
      std::cerr << "Fail to open file " << c.files[file_c].string() << std::endl;
      return 1;
    }
    idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
    if (idx[file_c] == NULL) {
      std::cerr << "Fail to open index for " << c.files[file_c].string() << std::endl;
      return 1;
    }
  }
  bam_hdr_t* hdr = sam_hdr_read(samfile[0]);
  if (hdr == NULL) {
    std::cerr << "Fail to open header for " << c.files[0].string() << std::endl;
    return 1;
  }

  // Calculate valid windows
  typedef uint32_t TPos;
  typedef std::set<TPos> TWindowSet;
  typedef std::vector<TWindowSet> TGenomicWindows;
  typedef std::vector<TGenomicWindows> TFileWindows;
  TFileWindows watsonWindows;
  watsonWindows.resize(c.files.size());
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) watsonWindows[file_c].resize(hdr->n_targets);
  TFileWindows crickWindows;
  crickWindows.resize(c.files.size());
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) crickWindows[file_c].resize(hdr->n_targets);
  TFileWindows wcWindows;
  wcWindows.resize(c.files.size());
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) wcWindows[file_c].resize(hdr->n_targets);

  // Parse bam (contig by contig)
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "BAM file parsing" << std::endl;
  boost::progress_display show_progress( c.files.size() );
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    ++show_progress;
    for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
      if (hdr->target_len[refIndex] < c.window) continue;
      uint32_t bins = hdr->target_len[refIndex] / c.window + 1;
      typedef std::vector<uint32_t> TCounter;
      TCounter watsonCount;
      TCounter crickCount;
      watsonCount.resize(bins, 0);
      crickCount.resize(bins, 0);
      
      hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
	
	int32_t pos = rec->core.pos + halfAlignmentLength(rec);
	if (rec->core.flag & BAM_FREAD1) 
	  if (rec->core.flag & BAM_FREVERSE) ++crickCount[(int) (pos / c.window)];
	  else ++watsonCount[(int) (pos / c.window)];
	else
	  if (rec->core.flag & BAM_FREVERSE) ++watsonCount[(int) (pos / c.window)];
	  else ++crickCount[(int) (pos / c.window)];
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
      
      // Get reads per window
      TCounter support;
      support.resize(bins, 0);
      TCounter::iterator itSupport = support.begin();
      TCounter::const_iterator itWatson = watsonCount.begin();
      TCounter::const_iterator itCrick =crickCount.begin();
      for(std::size_t bin = 0; bin < bins; ++itWatson, ++itCrick, ++itSupport, ++bin) 
	*itSupport = *itWatson + *itCrick;
      std::sort(support.begin(), support.end());
      uint32_t lowerCutoff = support[(int) (bins/100)];
      uint32_t upperCutoff = support[(int) (99*bins/100)];
      support.clear();
      
      // Get Watson Ratio
      itWatson = watsonCount.begin();
      itCrick =crickCount.begin();
      typedef std::vector<double> TRatio;
      TRatio wRatio;
      for(std::size_t bin = 0; bin < bins; ++itWatson, ++itCrick, ++bin) {
	uint32_t sup = *itWatson + *itCrick;
	// At least 1 read every 10,000 bases
	if ((sup > (c.window / 10000)) && (sup>lowerCutoff) && (sup<upperCutoff)) wRatio.push_back(((double) *itWatson / (double) (sup)));
      }
      std::sort(wRatio.begin(), wRatio.end());
      uint32_t wRatioSize = wRatio.size();
      if (!wRatioSize) continue;
      
      // Categorize chromosomes (exclude chromosomes with recombination events)
      double lowerW = wRatio[(int) (wRatioSize/10)];
      double upperW = wRatio[(int) (9*wRatioSize/10)];
      itWatson = watsonCount.begin();
      itCrick =crickCount.begin();
      if ((lowerW > 0.8) && (upperW > 0.8)) {
	// Hom. Watson
	for(std::size_t bin = 0; bin < bins; ++itWatson, ++itCrick, ++bin) {
	  uint32_t sup = *itWatson + *itCrick;
	  if ((sup > (c.window / 10000)) && (sup>lowerCutoff) && (sup<upperCutoff)) {
	    double wTmpRatio = (double) *itWatson / (double) (sup);
	    if (wTmpRatio > 0.8) watsonWindows[file_c][refIndex].insert(bin);
	  }
	}
      } else if ((lowerW < 0.2) && (upperW < 0.2)) {
	// Hom. Crick
	for(std::size_t bin = 0; bin < bins; ++itWatson, ++itCrick, ++bin) {
	  uint32_t sup = *itWatson + *itCrick;
	  if ((sup > (c.window / 10000)) && (sup>lowerCutoff) && (sup<upperCutoff)) {
	    double wTmpRatio = (double) *itWatson / (double) (sup);
	    if (wTmpRatio < 0.2) crickWindows[file_c][refIndex].insert(bin);
	  }
	}
      } else if ((lowerW > 0.3) && (upperW < 0.7)) {
	//WC
	for(std::size_t bin = 0; bin < bins; ++itWatson, ++itCrick, ++bin) {
	  uint32_t sup = *itWatson + *itCrick;
	  if ((sup > (c.window / 10000)) && (sup>lowerCutoff) && (sup<upperCutoff)) {
	    double wTmpRatio = (double) *itWatson / (double) (sup);
	    if ((wTmpRatio > 0.3) && (wTmpRatio < 0.7)) wcWindows[file_c][refIndex].insert(bin);
	  }
	}
      }
    }
  }
  

  // Write bam files
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "BAM writing" << std::endl;
  boost::progress_display spr( c.files.size() );
  samFile* wwbam = sam_open(c.ww.string().c_str(), "wb");
  if (wwbam == NULL) {
    std::cerr << "Fail to open file " << c.ww.string() << std::endl;
    return 1;
  }
  sam_hdr_write(wwbam, hdr);
  samFile* wcbam = sam_open(c.wc.string().c_str(), "wb");
  if (wcbam == NULL) {
    std::cerr << "Fail to open file " << c.wc.string() << std::endl;
    return 1;
  }
  sam_hdr_write(wcbam, hdr);
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    ++spr;
    for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
      if (hdr->target_len[refIndex] < c.window) continue;
      hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
	
	int32_t pos = rec->core.pos + halfAlignmentLength(rec);
	uint32_t bin = (int) (pos / c.window);
	if (watsonWindows[file_c][refIndex].find(bin) != watsonWindows[file_c][refIndex].end()) sam_write1(wwbam, hdr, rec);
	else if (crickWindows[file_c][refIndex].find(bin) != crickWindows[file_c][refIndex].end()) {
	  rec->core.flag ^= BAM_FREVERSE;
	  rec->core.flag ^= BAM_FMREVERSE;
	  sam_write1(wwbam, hdr, rec);
	} 
	else if (wcWindows[file_c][refIndex].find(bin) != wcWindows[file_c][refIndex].end()) sam_write1(wcbam, hdr, rec);
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }
  }
  sam_close(wcbam);
  sam_close(wwbam);


  // Close bam
  bam_hdr_destroy(hdr);
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    hts_idx_destroy(idx[file_c]);
    sam_close(samfile[file_c]);
  }

#ifdef PROFILE
  ProfilerStop();
#endif

  // End
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
  return 0;
}
