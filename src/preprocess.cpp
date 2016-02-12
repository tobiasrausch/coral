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
#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "strandutil.h"
#include "strandsnp.h"
#include "strandhaplo.h"

using namespace streq;

struct Config {
  bool hasVariationFile;
  unsigned short minMapQual;
  uint32_t window;
  boost::filesystem::path ww;
  boost::filesystem::path wc;
  boost::filesystem::path variation;
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
    ("map-qual,q", boost::program_options::value<unsigned short>(&c.minMapQual)->default_value(10), "min. mapping quality")
    ("window,w", boost::program_options::value<uint32_t>(&c.window)->default_value(1000000), "window length")
    ("samestrand,s", boost::program_options::value<boost::filesystem::path>(&c.ww)->default_value("ww.bam"), "output same strand bam")
    ("diffstrand,d", boost::program_options::value<boost::filesystem::path>(&c.wc)->default_value("wc.bam"), "output different strand bam")
    ("variation,v", boost::program_options::value<boost::filesystem::path>(&c.variation), "SNP VCF file to unify WC data (optional)")
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
  
  // Check variation VCF file
  if (vm.count("variation")) {
    if (!(boost::filesystem::exists(c.variation) && boost::filesystem::is_regular_file(c.variation) && boost::filesystem::file_size(c.variation))) {
      std::cerr << "Input SNP VCF file is missing: " << c.variation.string() << std::endl;
      return 1;
    }
    c.hasVariationFile = true;
  } else c.hasVariationFile = false;

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

  // Load variation data
  typedef std::vector<Snp> TSnpVector;
  typedef std::vector<TSnpVector> TGenomicSnps;
  TGenomicSnps snps;
  if (c.hasVariationFile) 
    if (_loadVariationData(c, hdr, snps)) return 1;

  // Variation counts
  typedef std::vector<RACount> TCountVector;
  typedef std::vector<TCountVector> TGenomicCounts;
  typedef std::vector<TGenomicCounts> TFileCounts;
  TFileCounts fCount;
  if (c.hasVariationFile) {
    fCount.resize(c.files.size());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      fCount[file_c].resize(hdr->n_targets);
      for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
	  if (hdr->target_len[refIndex] < c.window) continue;
	  fCount[file_c][refIndex].resize(snps[refIndex].size(), RACount());
      }
    }
  }

  // Watson Ratio counts
  typedef float TWRatio;
  typedef std::vector<TWRatio> TWRatioVector;
  TWRatioVector wRatio;
  typedef std::vector<TWRatioVector> TGenomicWRatio;
  typedef std::vector<TGenomicWRatio> TFileWRatio;
  TFileWRatio fWR;
  fWR.resize(c.files.size());
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) fWR[file_c].resize(hdr->n_targets);

  // Parse bam (contig by contig)
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "BAM file parsing" << std::endl;
  boost::progress_display show_progress( c.files.size() );
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    ++show_progress;
    for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
      if (hdr->target_len[refIndex] < c.window) continue;
      uint32_t bins = hdr->target_len[refIndex] / c.window + 1;
      fWR[file_c][refIndex].resize(bins);
      typedef std::vector<uint32_t> TCounter;
      TCounter watsonCount(bins, 0);
      TCounter crickCount(bins, 0);
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
	if (c.hasVariationFile) _refAltCount(rec, snps[refIndex], fCount[file_c][refIndex]);
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
            
      // Get Watson Ratio
      uint32_t lowerReadSupportCutoff = _getReadSupportPercentile(watsonCount, crickCount, 1);
      TCounter::const_iterator itWatson = watsonCount.begin();
      TCounter::const_iterator itCrick = crickCount.begin();
      for(std::size_t bin = 0; bin < bins; ++itWatson, ++itCrick, ++bin) {
	uint32_t sup = *itWatson + *itCrick;
	// At least 1 read every 10,000 bases
	if ((sup > (c.window / 10000)) && (sup > lowerReadSupportCutoff)) {
	  fWR[file_c][refIndex][bin] = ((float) *itWatson / (float) (sup));
	  wRatio.push_back(((float) *itWatson / (float) (sup)));
	} else fWR[file_c][refIndex][bin] = -1;
      }
    }
  }

  // Get global optimal watson and crick threshold
  float crickCut = 0;
  float watsonCut = 0;
  boost::tie(crickCut, watsonCut) = _biModalMinima(wRatio);
  wRatio.clear();

  // Black-listing variation/repeat/outlier bins across cells
  std::vector<float> watfraction;
  std::vector<float> crifraction;
  for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
    if (hdr->target_len[refIndex] < c.window) continue;
    uint32_t bins = hdr->target_len[refIndex] / c.window + 1;
    for(std::size_t bin = 0; bin < bins; ++bin) {
      uint32_t totalCells = 0;
      float crickFraction = 0;
      float watsonFraction = 0;
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	if (fWR[file_c][refIndex][bin] != -1) {
	  ++totalCells;
	  if (fWR[file_c][refIndex][bin] < crickCut) ++crickFraction;
	  else if (fWR[file_c][refIndex][bin] > watsonCut) ++watsonFraction;
	}
      }
      if (totalCells > (c.files.size() / 2)) {
	crickFraction /= (float) totalCells;
	watsonFraction /= (float) totalCells;
	crifraction.push_back(crickFraction);
	watfraction.push_back(watsonFraction);
      }
    }
  }
  float watMedianFrac = 0;
  float watMedianMAD = 0;
  float criMedianFrac = 0;
  float criMedianMAD = 0;
  boost::tie(watMedianFrac, watMedianMAD) = _getWWStrandBounds(watfraction);
  boost::tie(criMedianFrac, criMedianMAD) = _getWWStrandBounds(crifraction);
  float lowerWatFracBound = watMedianFrac - 3 * watMedianMAD;
  float upperWatFracBound = watMedianFrac + 3 * watMedianMAD;
  float lowerCriFracBound = criMedianFrac - 3 * criMedianMAD;
  float upperCriFracBound = criMedianFrac + 3 * criMedianMAD;
  watfraction.clear();
  crifraction.clear();

  // Black-list bins
  unsigned int numwhitelist=0;
  unsigned int numblacklist=0;  
  unsigned int numhighcov=0;
  unsigned int numlowcov=0;  
  for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
    if (hdr->target_len[refIndex] < c.window) continue;
    uint32_t bins = hdr->target_len[refIndex] / c.window + 1;
    for(std::size_t bin = 0; bin < bins; ++bin) {
      bool validBin = true;
      uint32_t totalCells = 0;
      float crickFraction = 0;
      float watsonFraction = 0;
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	if (fWR[file_c][refIndex][bin] != -1) {
	  ++totalCells;
	  if (fWR[file_c][refIndex][bin] < crickCut) ++crickFraction;
	  else if (fWR[file_c][refIndex][bin] > watsonCut) ++watsonFraction;
	  ++numhighcov;
	} else ++numlowcov;
      }
      if (totalCells > (c.files.size() / 2)) {
	crickFraction /= (float) totalCells;
	watsonFraction /= (float) totalCells;
	if ((crickFraction < lowerCriFracBound) || (watsonFraction < lowerWatFracBound) || (crickFraction > upperCriFracBound) || (watsonFraction > upperWatFracBound)) validBin = false;
      } else validBin = false;
      if (!validBin) {
	++numblacklist;
	for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) fWR[file_c][refIndex][bin] = -1;
      }	else ++numwhitelist;
    }
  }

  // Categorize chromosomes based on white-listed bins
  typedef uint32_t TPos;
  typedef std::set<TPos> TWindowSet;
  typedef std::vector<TWindowSet> TGenomicWindows;
  typedef std::vector<TGenomicWindows> TFileWindows;
  TFileWindows watsonWindows(c.files.size());
  TFileWindows crickWindows(c.files.size());
  TFileWindows wcWindows(c.files.size());
  TFileWindows wcFlipWindows(c.files.size());
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    watsonWindows[file_c].resize(hdr->n_targets);
    crickWindows[file_c].resize(hdr->n_targets);
    wcWindows[file_c].resize(hdr->n_targets);
    wcFlipWindows[file_c].resize(hdr->n_targets);
  }
  uint32_t wWindowCount = 0;
  uint32_t cWindowCount = 0;
  uint32_t wcWindowCount = 0;
  float watsonBound = std::min((float) 0.9, watsonCut + 3 * watMedianMAD);
  float crickBound = std::max((float) 0.1, crickCut - 3 * criMedianMAD);
  float lWCBound = std::min((float) 0.4, std::max(crickCut + 3 * criMedianMAD, (float) 0.25));
  float uWCBound = std::max((float) 0.6, std::min(watsonCut - 3 * watMedianMAD, (float) 0.75));
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
      if (hdr->target_len[refIndex] < c.window) continue;
      uint32_t bins = hdr->target_len[refIndex] / c.window + 1;
      TWRatioVector chrWRatio;
      for(std::size_t bin = 0; bin < bins; ++bin)
	if (fWR[file_c][refIndex][bin] != -1) chrWRatio.push_back(fWR[file_c][refIndex][bin]);
      if ((chrWRatio.empty()) || (chrWRatio.size() < bins/2)) continue;

      // Categorize chromosomes (exclude chromosomes with recombination events)
      std::sort(chrWRatio.begin(), chrWRatio.end());
      TWRatio lowerW = chrWRatio[(int) (chrWRatio.size()/10)];
      TWRatio upperW = chrWRatio[(int) (9*chrWRatio.size()/10)];
      if (lowerW >= watsonBound) {
	// Hom. Watson
	for(std::size_t bin = 0; bin < bins; ++bin) {
	  watsonWindows[file_c][refIndex].insert(bin);
	  ++wWindowCount;
	}
      } else if (upperW <= crickBound) {
	// Hom. Crick
	for(std::size_t bin = 0; bin < bins; ++bin) {
	  crickWindows[file_c][refIndex].insert(bin);
	  ++cWindowCount;
	}
      } else if ((lowerW >= lWCBound) && (upperW <= uWCBound)) {
	//WC
	for(std::size_t bin = 0; bin < bins; ++bin) {
	  wcWindows[file_c][refIndex].insert(bin);
	  ++wcWindowCount;
	}
      }
    }
  }

  // Process variation data
  if (c.hasVariationFile) _processVariationData(c, hdr, snps, fCount, fWR, wcWindows, wcFlipWindows);
  
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
	else if (wcFlipWindows[file_c][refIndex].find(bin) != wcFlipWindows[file_c][refIndex].end()) {
	  rec->core.flag ^= BAM_FREVERSE;
	  rec->core.flag ^= BAM_FMREVERSE;
	  sam_write1(wcbam, hdr, rec);
	}
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

  // End
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;

  // Output statistics
  std::cout << "Threshold statistics: Crick cut=" << crickCut << ", Median Crick fraction=" << criMedianFrac << ", MAD=" << criMedianMAD << std::endl;
  std::cout << "Threshold statistics: Watson cut=" << watsonCut << ", Median Watson fraction=" << watMedianFrac << ", MAD=" << watMedianMAD << std::endl;
  std::cout << "Coverage statistics: #HighCoverageWindows=" << numhighcov << ", #LowCoverageWindows=" << numlowcov << std::endl;
  std::cout << "Bin statistics: #Whitelist=" << numwhitelist << ", #Blacklist=" << numblacklist << std::endl;
  std::cout << "Watson-Watson: Range=[" << watsonBound << ",1], #WatsonWindows=" << wWindowCount << std::endl;
  std::cout << "Crick-Crick: Range=[0," << crickBound << "], #CrickWindows=" << cWindowCount << std::endl;
  std::cout << "Watson-Crick: Range=[" << lWCBound << "," << uWCBound << "], #WatsonCrickWindows=" << wcWindowCount << std::endl;


#ifdef PROFILE
  ProfilerStop();
#endif


  return 0;
}
