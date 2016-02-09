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
#include <htslib/vcf.h>

#include "util.h"

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

struct Snp {
  uint32_t pos;
  char ref;
  char alt;
  
  Snp(uint32_t p, char r, char a) : pos(p), ref(r), alt(a) {}
};

struct RACount {
  uint16_t watsonRef;
  uint16_t watsonAlt;
  uint16_t crickRef;
  uint16_t crickAlt;
  uint16_t falseBase;
  
  RACount() : watsonRef(0), watsonAlt(0), crickRef(0), crickAlt(0), falseBase(0) {}
  RACount(uint16_t wr, uint16_t wa, uint16_t cr, uint16_t ca, uint16_t f) : watsonRef(wr), watsonAlt(wa), crickRef(cr), crickAlt(ca), falseBase(f) {}
};


template<typename TRecord>
struct SortSnps : public std::binary_function<TRecord, TRecord, bool>
{
  inline bool operator()(TRecord const& s1, TRecord const& s2) const {
    return s1.pos < s2.pos;
  }
};


template<typename TIterator, typename TValue>
inline void
_getMedian(TIterator begin, TIterator end, TValue& median) {
  std::nth_element(begin, begin + (end - begin) / 2, end);
  median = *(begin + (end - begin) / 2);
}

template<typename TIterator, typename TValue>
inline void
_getMAD(TIterator begin, TIterator end, TValue median, TValue& mad) {
  std::vector<TValue> absDev;
  for(;begin<end;++begin)
    absDev.push_back(std::abs((TValue)*begin - median));
  _getMedian(absDev.begin(), absDev.end(), mad);
}

template<typename TConfig, typename TGenomicSnps>
inline int32_t 
_loadVariationData(TConfig const& c, bam_hdr_t const* hdr, TGenomicSnps& snps) {
  typedef typename TGenomicSnps::value_type TSnpVector;
  typedef typename TSnpVector::value_type TSnp;

  // Open VCF file
  htsFile* ifile = bcf_open(c.variation.string().c_str(), "r");
  if (ifile == NULL) {
    std::cerr << "SNP VCF files is missing " << c.variation.string() << std::endl;
    return 1;
  }
  hts_idx_t* bcfidx = bcf_index_load(c.variation.string().c_str());
  if (bcfidx == NULL) {
    std::cerr << "SNP VCF index file is missing " << c.variation.string() << std::endl;
    return 1;
  }
  bcf_hdr_t* vcfh = bcf_hdr_read(ifile);
    
  // Load SNPs
  snps.clear();
  snps.resize(hdr->n_targets);
  for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
    if (hdr->target_len[refIndex] < c.window) continue;
    std::string chrName(hdr->target_name[refIndex]);
    uint32_t chrid = bcf_hdr_name2id(vcfh, chrName.c_str());
    hts_itr_t* itervcf = bcf_itr_queryi(bcfidx, chrid, 0, hdr->target_len[refIndex]);
    bcf1_t* var = bcf_init();
    while (bcf_itr_next(ifile, itervcf, var) >= 0) {
      bcf_unpack(var, BCF_UN_STR);
      std::vector<std::string> alleles;
      for(std::size_t i = 0; i<var->n_allele; ++i) alleles.push_back(std::string(var->d.allele[i]));
      // Only bi-allelic SNPs
      if ((alleles.size() == 2) && (alleles[0].size() == 1) && (alleles[1].size() == 1)) snps[refIndex].push_back(TSnp(var->pos, alleles[0][0], alleles[1][0]));
    }
    bcf_destroy(var);
    hts_itr_destroy(itervcf);

    // Sort Snps
    std::sort(snps[refIndex].begin(), snps[refIndex].end(), SortSnps<TSnp>());
  }
  // Close VCF
  bcf_hdr_destroy(vcfh);
  hts_idx_destroy(bcfidx);
  bcf_close(ifile);

  return 0;
}

template<typename TSnpVector, typename TCountVector>
inline void
_refAltCount(bam1_t const* rec, TSnpVector const& snps, TCountVector& chrCounts) {
  typedef typename TSnpVector::value_type TSnp;
  
  // Annotate SNPs
  typename TSnpVector::const_iterator iSnp = std::lower_bound(snps.begin(), snps.end(), TSnp(rec->core.pos, 'A', 'A'), SortSnps<TSnp>());
  typename TSnpVector::const_iterator iSnpEnd = std::upper_bound(snps.begin(), snps.end(), TSnp(rec->core.pos + alignmentLength(rec), 'A', 'A'), SortSnps<TSnp>());
  if (iSnp != iSnpEnd) {
    std::string sequence;
    sequence.resize(rec->core.l_qseq);
    uint8_t* seqptr = bam_get_seq(rec);
    for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
    uint32_t* cigar = bam_get_cigar(rec);
    for(;iSnp != iSnpEnd; ++iSnp) {
      int32_t gp = rec->core.pos; // Genomic position
      int32_t sp = 0; // Sequence position
      bool foundChar = false;
      for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) sp += bam_cigar_oplen(cigar[i]);
	else if (bam_cigar_op(cigar[i]) == BAM_CINS) sp += bam_cigar_oplen(cigar[i]);
	else if (bam_cigar_op(cigar[i]) == BAM_CDEL) gp += bam_cigar_oplen(cigar[i]);
	else if (bam_cigar_op(cigar[i]) == BAM_CMATCH) {
	  if (gp + bam_cigar_oplen(cigar[i]) < iSnp->pos) {
	    gp += bam_cigar_oplen(cigar[i]);
	    sp += bam_cigar_oplen(cigar[i]);
	  } else {
	    for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]); ++k, ++sp, ++gp) {
	      if (gp == iSnp->pos) {
		foundChar = true;
		break;
	      }
	    }
	    if (foundChar) break;
	  }
	}
      }
      if (foundChar) {
	if (sequence[sp] == iSnp->ref) { 
	  if (rec->core.flag & BAM_FREAD1) 
	    if (rec->core.flag & BAM_FREVERSE) ++chrCounts[(uint32_t) (iSnp - snps.begin())].crickRef;
	    else ++chrCounts[(uint32_t) (iSnp - snps.begin())].watsonRef;
	  else
	    if (rec->core.flag & BAM_FREVERSE) ++chrCounts[(uint32_t) (iSnp - snps.begin())].watsonRef;
	    else ++chrCounts[(uint32_t) (iSnp - snps.begin())].crickRef;
	} else if (sequence[sp] == iSnp->alt) {
	  if (rec->core.flag & BAM_FREAD1) 
	    if (rec->core.flag & BAM_FREVERSE) ++chrCounts[(uint32_t) (iSnp - snps.begin())].crickAlt;
	    else ++chrCounts[(uint32_t) (iSnp - snps.begin())].watsonAlt;
	  else
	    if (rec->core.flag & BAM_FREVERSE) ++chrCounts[(uint32_t) (iSnp - snps.begin())].watsonAlt;
	    else ++chrCounts[(uint32_t) (iSnp - snps.begin())].crickAlt;
	} else ++chrCounts[(uint32_t) (iSnp - snps.begin())].falseBase;
      }
    }
  }
}


template<typename TRefAltCount>
inline bool
_getWatsonAllele(TRefAltCount const& ra, bool& allele) {
  uint32_t watsonCount = ra.watsonRef + ra.watsonAlt;
  if (watsonCount > 0) {
    if (ra.watsonRef > 2*ra.watsonAlt) {
      allele = false;
      return true;
    } else if (2 * ra.watsonRef < ra.watsonAlt) {
      allele = true;
      return true;
    }
  }
  return false;
}

template<typename TRefAltCount>
inline bool
_getCrickAllele(TRefAltCount const& ra, bool& allele) {
  uint32_t crickCount = ra.crickRef + ra.crickAlt;
  if (crickCount > 0) {
    if (ra.crickRef > 2*ra.crickAlt) {
      allele = false;
      return true;
    } else if (2 * ra.crickRef < ra.crickAlt) {
      allele = true;
      return true;
    }
  }
  return false;
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

  // Parse bam (contig by contig)
  typedef float TWRatio;
  typedef std::vector<TWRatio> TWRatioVector;
  typedef std::vector<TWRatioVector> TGenomicWRatio;
  typedef std::vector<TGenomicWRatio> TFileWRatio;
  TFileWRatio fWR;
  fWR.resize(c.files.size());
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) fWR[file_c].resize(hdr->n_targets);
  TWRatioVector wRatio;
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
	if (c.hasVariationFile) _refAltCount(rec, snps[refIndex], fCount[file_c][refIndex]);
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
      support.clear();
      
      // Get Watson Ratio
      itWatson = watsonCount.begin();
      itCrick =crickCount.begin();
      for(std::size_t bin = 0; bin < bins; ++itWatson, ++itCrick, ++bin) {
	uint32_t sup = *itWatson + *itCrick;
	// At least 1 read every 10,000 bases
	if ((sup > (c.window / 10000)) && (sup>lowerCutoff)) {
	  fWR[file_c][refIndex][bin] = ((float) *itWatson / (float) (sup));
	  wRatio.push_back(((float) *itWatson / (float) (sup)));
	} else fWR[file_c][refIndex][bin] = -1;
      }
    }
  }

  // Generate Watson Ratio Statistics
  std::sort(wRatio.begin(), wRatio.end());
  double crickCut = wRatio[(int) (wRatio.size()/4)];
  double watsonCut = wRatio[(int) (3*wRatio.size()/4)];
  std::vector<float> homfraction;
  for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
    if (hdr->target_len[refIndex] < c.window) continue;
    uint32_t bins = hdr->target_len[refIndex] / c.window + 1;
    for(std::size_t bin = 0; bin < bins; ++bin) {
      TWRatioVector binWRatio;
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	if (fWR[file_c][refIndex][bin] != -1) binWRatio.push_back(fWR[file_c][refIndex][bin]);
      }
      if (binWRatio.size() > (c.files.size() / 2)) {
	float crickFraction = 0;
	float watsonFraction = 0;
	for(std::size_t k = 0; k<binWRatio.size(); ++k) {
	  if (binWRatio[k]<crickCut) ++crickFraction;
	  if (binWRatio[k]>watsonCut) ++watsonFraction;
	}
	crickFraction /= (float) binWRatio.size();
	watsonFraction /= (float) binWRatio.size();
	homfraction.push_back(crickFraction);
	homfraction.push_back(watsonFraction);
      }
    }
  }
  float hommedian = 0;
  _getMedian(homfraction.begin(), homfraction.end(), hommedian);
  float hommad = 0;
  _getMAD(homfraction.begin(), homfraction.end(), hommedian, hommad);
  std::cout << "Strand-Seq statistics: Crick cut=" << crickCut << ", Watson cut=" << watsonCut << ", Median Watson and Crick fraction=" << hommedian << ", MAD=" << hommad << std::endl;
  float fracDev = 3 * hommad;
  float upperFracBound = hommedian + fracDev;
  float lowerFracBound = hommedian - fracDev;
  
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
      TWRatioVector binWRatio;
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	if (fWR[file_c][refIndex][bin] != -1) {
	  binWRatio.push_back(fWR[file_c][refIndex][bin]);
	  ++numhighcov;
	} else ++numlowcov;
      }
      if (binWRatio.size() > (c.files.size() / 2)) {
	float crickFraction = 0;
	float watsonFraction = 0;
	for(std::size_t k = 0; k<binWRatio.size(); ++k) {
	  if (binWRatio[k]<crickCut) ++crickFraction;
	  if (binWRatio[k]>watsonCut) ++watsonFraction;
	}
	crickFraction /= (float) binWRatio.size();
	watsonFraction /= (float) binWRatio.size();
	if ((crickFraction < lowerFracBound) || (watsonFraction < lowerFracBound) || (crickFraction > upperFracBound) || (watsonFraction > upperFracBound)) validBin = false;
      } else validBin = false;
      if (!validBin) {
	++numblacklist;
	for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) fWR[file_c][refIndex][bin] = -1;
      }	else ++numwhitelist;
    }
  }
  std::cout << "Bin statistics: #Whitelist=" << numwhitelist << ", #Blacklist=" << numblacklist << ", #HighCoverageWindows=" << numhighcov << ", #LowCoverageWindows=" << numlowcov << std::endl;

  // Categorize chromosomes based on white-listed bins
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
  TFileWindows wcFlipWindows;
  wcFlipWindows.resize(c.files.size());
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) wcFlipWindows[file_c].resize(hdr->n_targets);

  uint32_t wWindowCount = 0;
  uint32_t cWindowCount = 0;
  uint32_t wcWindowCount = 0;
  double watsonBound = std::min(0.9, watsonCut + fracDev);
  double crickBound = std::max(0.1, crickCut - fracDev);
  double lWCBound = std::max(crickCut + fracDev, 0.25);
  lWCBound = std::min(0.4, lWCBound);
  double uWCBound = std::min(watsonCut - fracDev, 0.75);
  uWCBound = std::max(0.6, uWCBound);
  for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
    for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
      if (hdr->target_len[refIndex] < c.window) continue;
      uint32_t bins = hdr->target_len[refIndex] / c.window + 1;
      TWRatioVector chrWRatio;
      for(std::size_t bin = 0; bin < bins; ++bin) {
	if (fWR[file_c][refIndex][bin] != -1) chrWRatio.push_back(fWR[file_c][refIndex][bin]);
      }
      // Debug chromosomal watson ratios
      //std::cerr << refIndex << '\t' << file_c;
      //for(int i = 0; i < chrWRatio.size(); ++i) std::cerr << '\t' << chrWRatio[i];
      //std::cerr << std::endl;
      std::sort(chrWRatio.begin(), chrWRatio.end());
      uint32_t chrWRatioSize = chrWRatio.size();
      if ((!chrWRatioSize) || (chrWRatioSize < bins/2)) continue;

      // Categorize chromosomes (exclude chromosomes with recombination events)
      double lowerW = chrWRatio[(int) (chrWRatioSize/10)];
      double upperW = chrWRatio[(int) (9*chrWRatioSize/10)];
      if (lowerW >= watsonBound) {
	// Hom. Watson
	for(std::size_t bin = 0; bin < bins; ++bin) {
	  //uint32_t sup = *itWatson + *itCrick;
	  //if ((sup > (c.window / 10000)) && (sup>lowerCutoff)) {
	  //double wTmpRatio = (double) *itWatson / (double) (sup);
	  //if (wTmpRatio >= 0.8) 
	  watsonWindows[file_c][refIndex].insert(bin);
	  ++wWindowCount;
	}
      } else if (upperW <= crickBound) {
	// Hom. Crick
	for(std::size_t bin = 0; bin < bins; ++bin) {
	  //uint32_t sup = *itWatson + *itCrick;
	  //if ((sup > (c.window / 10000)) && (sup>lowerCutoff)) {
	  /// double wTmpRatio = (double) *itWatson / (double) (sup);
	  //if (wTmpRatio <= 0.2) crickWindows[file_c][refIndex].insert(bin);
	  //}
	  crickWindows[file_c][refIndex].insert(bin);
	  ++cWindowCount;
	}
      } else if ((lowerW >= lWCBound) && (upperW <= uWCBound)) {
	//WC
	for(std::size_t bin = 0; bin < bins; ++bin) {
      // uint32_t sup = *itWatson + *itCrick;
      //  if ((sup > (c.window / 10000)) && (sup>lowerCutoff)) {
      //    double wTmpRatio = (double) *itWatson / (double) (sup);
      //    if ((wTmpRatio >= 0.4) && (wTmpRatio <= 0.6)) wcWindows[file_c][refIndex].insert(bin);
      //  }
	  wcWindows[file_c][refIndex].insert(bin);
	  ++wcWindowCount;
	}
      }
    }
  }
  std::cout << "Watson-Watson: Range=[" << watsonBound << ",1], #WatsonWindows=" << wWindowCount << std::endl;
  std::cout << "Crick-Crick: Range=[0," << crickBound << "], #CrickWindows=" << cWindowCount << std::endl;
  std::cout << "Watson-Crick: Range=[" << lWCBound << "," << uWCBound << "], #WatsonCrickWindows=" << wcWindowCount << std::endl;
  
  // Process variation data
  if (c.hasVariationFile) {
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Process variation data" << std::endl;
    boost::progress_display sprog( hdr->n_targets);
    for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
      ++sprog;
      if (hdr->target_len[refIndex] < c.window) continue;
      // Sum-up ref and alt counts
      TCountVector sumCounts;
      sumCounts.resize(snps[refIndex].size(), RACount());
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	TCountVector::iterator itSC = sumCounts.begin();
	for(TCountVector::const_iterator itC = fCount[file_c][refIndex].begin(); itC != fCount[file_c][refIndex].end(); ++itC, ++itSC) {
	  itSC->watsonRef += itC->watsonRef;
	  itSC->watsonAlt += itC->watsonAlt;
	  itSC->crickRef += itC->crickRef;
	  itSC->crickAlt += itC->crickAlt;
	  itSC->falseBase += itC->falseBase;
	}
      }
      
      // Find informative het. SNPs
      typedef std::set<std::size_t> THetSnp;
      THetSnp hetSNP;
      for(TCountVector::const_iterator itSC = sumCounts.begin(); itSC != sumCounts.end(); ++itSC) {
	// Ref and alt seen (true variation) + watson and crick seen (no strand bias)?
	if ((itSC->falseBase == 0) && ((itSC->watsonRef + itSC->crickRef) > 0) && ((itSC->watsonAlt + itSC->crickAlt) > 0) && ((itSC->watsonRef + itSC->watsonAlt) > 0) && ((itSC->crickRef + itSC->crickAlt) > 0)) hetSNP.insert((std::size_t) (itSC - sumCounts.begin()));
      }

      // Build Haplotypes for wcWindows
      std::vector<uint32_t> support;
      support.resize(c.files.size(), 0);
      typedef std::vector<bool> TFlipVector;
      typedef std::vector<TFlipVector> TFlipMatrix;
      TFlipMatrix fm;
      fm.resize(c.files.size(), TFlipVector());
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) fm[file_c].resize(c.files.size(), false);
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	if (wcWindows[file_c][refIndex].empty()) continue;
	for(unsigned int file_d = file_c + 1; file_d < c.files.size(); ++file_d) {
	  if (wcWindows[file_d][refIndex].empty()) continue;
	  uint32_t diffCount = 0;
	  uint32_t agreeCount = 0;
	  // Iterate informative SNPs
	  for(THetSnp::const_iterator itHS = hetSNP.begin(); itHS != hetSNP.end(); ++itHS) {
	    uint32_t curBin = snps[refIndex][*itHS].pos / c.window;
	    if ((fWR[file_c][refIndex][curBin] == -1) || (fWR[file_d][refIndex][curBin] == -1)) continue; // Blacklisted bin
	    if ((wcWindows[file_c][refIndex].find(curBin) != wcWindows[file_c][refIndex].end()) && (wcWindows[file_d][refIndex].find(curBin) != wcWindows[file_d][refIndex].end())) {
	      bool cWAllele = false;
	      bool watsonSuccess = _getWatsonAllele(fCount[file_c][refIndex][*itHS], cWAllele);
	      bool cCAllele = false;
	      bool crickSuccess = _getCrickAllele(fCount[file_c][refIndex][*itHS], cCAllele);
	      // At least one allele must be called
	      if ((watsonSuccess) && (!crickSuccess)) cCAllele = !cWAllele;
	      else if ((!watsonSuccess) && (crickSuccess)) cWAllele = !cCAllele;
	      else if ((!watsonSuccess) && (!crickSuccess)) continue;  // Not covered
	      else if (cWAllele == cCAllele) continue; // Incorrect genotyping
	      bool dWAllele = false;
	      watsonSuccess = _getWatsonAllele(fCount[file_d][refIndex][*itHS], dWAllele);
	      bool dCAllele = false;
	      crickSuccess = _getCrickAllele(fCount[file_d][refIndex][*itHS], dCAllele);
	      // At least one allele must be called
	      if ((watsonSuccess) && (!crickSuccess)) dCAllele = !dWAllele;
	      else if ((!watsonSuccess) && (crickSuccess)) dWAllele = !dCAllele;
	      else if ((!watsonSuccess) && (!crickSuccess)) continue;  // Not covered
	      else if (dWAllele == dCAllele) continue; // Incorrect genotyping

	      // Same haplotype?
	      if (cWAllele == dWAllele) ++agreeCount;
	      else ++diffCount;

	      // Debug code
	      //std::cerr << file_c << ',' << refIndex << ',' << snps[refIndex][*itHS].pos << ',' <<  snps[refIndex][*itHS].ref << ',' <<  snps[refIndex][*itHS].alt << ',' << fCount[file_c][refIndex][*itHS].watsonRef << ',' << fCount[file_c][refIndex][*itHS].watsonAlt << ',' << fCount[file_c][refIndex][*itHS].crickRef << ',' << fCount[file_c][refIndex][*itHS].crickAlt << ',' << std::endl;
	      //std::cerr << file_d << ',' << refIndex << ',' << snps[refIndex][*itHS].pos << ',' <<  snps[refIndex][*itHS].ref << ',' <<  snps[refIndex][*itHS].alt << ',' << fCount[file_d][refIndex][*itHS].watsonRef << ',' << fCount[file_d][refIndex][*itHS].watsonAlt << ',' << fCount[file_d][refIndex][*itHS].crickRef << ',' << fCount[file_d][refIndex][*itHS].crickAlt << ',' << std::endl;
	    }
	  }
	  if (diffCount>agreeCount) {
	    support[file_c] += diffCount;
	    support[file_d] += diffCount;
	    fm[file_c][file_d] = true;
	    fm[file_d][file_c] = true;
	  } else {
	    support[file_c] += agreeCount;
	    support[file_d] += agreeCount;
	  }
	}
      }
      // Find best covered sample (most informative het. SNPs)
      uint32_t bestSupport = 0;
      uint32_t bestIndex = 0;
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	if (support[file_c] > bestSupport) {
	  bestSupport = support[file_c];
	  bestIndex = file_c;
	}
      }
      if (bestSupport > 0) {
	for(unsigned int file_d = 0; file_d < c.files.size(); ++file_d) {
	  if (file_d == bestIndex) continue; // No flip required
	  if (fm[bestIndex][file_d]) {
	    // Flip required
	    if (wcWindows[file_d][refIndex].empty()) continue;
	    uint32_t bins = hdr->target_len[refIndex] / c.window + 1;
	    for(uint32_t bin = 0; bin < bins; ++bin) {
	      if (wcWindows[file_d][refIndex].find(bin) != wcWindows[file_d][refIndex].end()) {
		wcFlipWindows[file_d][refIndex].insert(bin);
		wcWindows[file_d][refIndex].erase(bin);
	      }
	    }
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

#ifdef PROFILE
  ProfilerStop();
#endif

  // End
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
  return 0;
}
