/*
============================================================================
Strand-Seq Watson-Crick Counter
============================================================================
Copyright (C) 2015 Tobias Rausch

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


struct Config {
  unsigned short minMapQual;
  uint32_t window;
  boost::filesystem::path preOutput;
  std::vector<boost::filesystem::path> files;
};

inline int32_t halfAlignmentLength(bam1_t* rec) {
  uint32_t* cigar = bam_get_cigar(rec);
  unsigned int alen = 0;
  for (unsigned int i = 0; i < rec->core.n_cigar; ++i)
    if (bam_cigar_op(cigar[i]) == BAM_CMATCH) alen+=bam_cigar_oplen(cigar[i]);
  return (alen / 2);
}


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
    ("preout,p", boost::program_options::value<boost::filesystem::path>(&c.preOutput)->default_value("pre.out"), "output preprocessing info")
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
  if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("preout"))) {
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

  // Parse bam (contig by contig)
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Bam file parsing" << std::endl;
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
	    if (wTmpRatio > 0.8) watsonWindows[refIndex].push_back(std::make_pair(bin, file_c));
	  }
	}
      } else if ((lowerW < 0.2) && (upperW < 0.2)) {
	// Hom. Crick
	for(std::size_t bin = 0; bin < bins; ++itWatson, ++itCrick, ++bin) {
	  uint32_t sup = *itWatson + *itCrick;
	  if ((sup > (c.window / 10000)) && (sup>lowerCutoff) && (sup<upperCutoff)) {
	    double wTmpRatio = (double) *itWatson / (double) (sup);
	    if (wTmpRatio < 0.2) crickWindows[refIndex].push_back(std::make_pair(bin, file_c));
	  }
	}
      } else if ((lowerW > 0.3) && (upperW < 0.7)) {
	//WC
	for(std::size_t bin = 0; bin < bins; ++itWatson, ++itCrick, ++bin) {
	  uint32_t sup = *itWatson + *itCrick;
	  if ((sup > (c.window / 10000)) && (sup>lowerCutoff) && (sup<upperCutoff)) {
	    double wTmpRatio = (double) *itWatson / (double) (sup);
	    if ((wTmpRatio > 0.3) && (wTmpRatio < 0.7)) wcWindows[refIndex].push_back(std::make_pair(bin, file_c));
	  }
	}
      }
    }
  }
  

  // Sort windows
  for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
    std::sort(watsonWindows[refIndex].begin(), watsonWindows[refIndex].end());
    std::sort(crickWindows[refIndex].begin(), crickWindows[refIndex].end());
    std::sort(wcWindows[refIndex].begin(), wcWindows[refIndex].end());
  }

  // Dump pre-processing information
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Output pre-processing information" << std::endl;
  boost::progress_display sp( hdr->n_targets );

  std::ofstream ofile(c.preOutput.string().c_str());
  for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
    ++sp;
    for(TWindowList::iterator iW = watsonWindows[refIndex].begin(); iW != watsonWindows[refIndex].end(); ++iW) ofile << "0\t" << refIndex << '\t' << iW->first << '\t' << iW->second << '\t' << std::endl;
    for(TWindowList::iterator iC = crickWindows[refIndex].begin(); iC != crickWindows[refIndex].end(); ++iC) ofile << "1\t" << refIndex << '\t' << iC->first << '\t' << iC->second << '\t' << std::endl;
    for(TWindowList::iterator iWC = wcWindows[refIndex].begin(); iWC != wcWindows[refIndex].end(); ++iWC) ofile << "2\t" << refIndex << '\t' << iWC->first << '\t' << iWC->second << '\t' << std::endl;
  }
  ofile.close();

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
