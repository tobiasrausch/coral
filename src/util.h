/*
============================================================================
Coral: COpy-numbeR ALterations
============================================================================
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

#ifndef UTIL_H
#define UTIL_H

#include <boost/filesystem.hpp>
#include <boost/multi_array.hpp>
#include <boost/progress.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/zlib.hpp>

#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/sam.h>

#include "matrix.h"
#include "gflars.h"

namespace coralns
{

  #ifndef LAST_BIN
  #define LAST_BIN 65535
  #endif

  struct ScanWindow {
    bool select;
    int32_t start;
    int32_t end;
    uint32_t cov;
    uint32_t uniqcov;

    ScanWindow() : select(false), start(0), end(0), cov(0), uniqcov(0) {}
    explicit ScanWindow(int32_t const s) : select(false), start(s), end(s+1), cov(0), uniqcov(0) {}
  };

  template<typename TScanWindow>
  struct SortScanWindow : public std::binary_function<TScanWindow, TScanWindow, bool>
  {
    inline bool operator()(TScanWindow const& sw1, TScanWindow const& sw2) {
      return ((sw1.start<sw2.start) || ((sw1.start == sw2.start) && (sw1.end < sw2.end)));
    }
  };


  template<typename TConfig>
  inline int32_t
  _findScanWindow(TConfig const& c, uint32_t const reflen, std::vector<uint16_t> const& binMap, int32_t const midPoint) {
    if (c.hasScanFile) {
      if (binMap[midPoint] == LAST_BIN) return -1;
      else return binMap[midPoint];
    } else {
      uint32_t bin = midPoint / c.scanWindow;
      uint32_t allbins = reflen / c.scanWindow;
      if (bin >= allbins) return -1;
      else return bin;
    }
    return -1;
  }
    
  struct LibraryInfo {
    int32_t rs;
    int32_t median;
    int32_t mad;
    int32_t minNormalISize;
    int32_t maxNormalISize;

    LibraryInfo() : rs(0), median(0), mad(0), minNormalISize(0), maxNormalISize(0) {}
  };


  struct SegmentConfig {
    typedef double TPrecision;
    typedef boost::multi_array<TPrecision, 2> TSignalMatrix;
    
    uint32_t k;
    double epsilon;
    double dpthreshold;
    std::string outprefix;
    boost::filesystem::path signal;
  };

  
  struct Interval {
    uint32_t istart;
    uint32_t iend;
    
    Interval() : istart(0), iend(0) {}
    Interval(uint32_t a, uint32_t b) : istart(a), iend(b) {}
  };  
  
  struct NormalizedBinCounts {
    typedef SegmentConfig::TSignalMatrix TSignalMatrix;
    typedef std::vector<Interval> TIntervals;
    std::string chr;
    uint32_t rows;
    uint32_t cols;
    TIntervals itv;
    TSignalMatrix sm;
    
    NormalizedBinCounts() : chr(""), rows(0), cols(0) {}
  };


  struct SmoothSignal {
    typedef Recap::TPrecision TPrecision;
    typedef std::vector<uint32_t> TIndexVector;
    typedef boost::multi_array<TPrecision, 2> TSignalMatrix;

    TIndexVector jumps;
    TSignalMatrix smooth;
    TSignalMatrix updown;
  };
  
  template<typename TSignalMatrix>
  inline void
  smoothsignal(TSignalMatrix const& sm, std::vector<uint32_t> const& jumps, SmoothSignal& res) {
    typedef Recap::TPrecision TPrecision;
    uint32_t nrow = sm.shape()[0];
    uint32_t ncol = sm.shape()[1];
    
    std::vector<int32_t> b(jumps.begin(), jumps.end());
    b.push_back(-1);
    b.push_back(nrow-1);
    std::sort(b.begin(), b.end());

    res.jumps.clear();
    for(uint32_t i = 1; i < b.size(); ++i) res.jumps.push_back(b[i]);
    uint32_t k = res.jumps.size();
    res.smooth.resize(boost::extents[k][ncol]);
    for(uint32_t i = 0; i < k; ++i) {
      uint32_t istart = b[i] + 1;
      uint32_t iend = b[i+1] + 1;
      for(uint32_t j = 0; j < ncol; ++j) {
	TPrecision avg = 0;
	for(uint32_t ki = istart; ki<iend; ++ki) avg += sm[ki][j];
	avg /= (TPrecision) (iend - istart);
	res.smooth[i][j] = avg;
      }
    }
  }

  inline void
  updown(SmoothSignal& res) {
    typedef Recap::TPrecision TPrecision;
    uint32_t nrow = res.smooth.shape()[0];
    uint32_t ncol = res.smooth.shape()[1];
    res.updown.resize(boost::extents[nrow][2]);
    for(uint32_t i = 0; i<nrow; ++i) {
      uint32_t upcount = 0;
      TPrecision up = 0;
      uint32_t downcount = 0;
      TPrecision down = 0;
      for(uint32_t j = 0; j<ncol; ++j) {
	if (res.smooth[i][j] > 0) {
	  up += res.smooth[i][j];
	  ++upcount;
	} else if (res.smooth[i][j] < 0) {
	  down += res.smooth[i][j];
	  ++downcount;
	}
      }
      if (upcount) res.updown[i][0] = up / (TPrecision) upcount;
      else res.updown[i][0] = 0;
      if (downcount) res.updown[i][1] = down / (TPrecision) downcount;
      else res.updown[i][1] = 0;
    }
  }

  inline bool
  nContent(std::string const& s) {
    for(uint32_t i = 0; i < s.size(); ++i) {
      if ((s[i] == 'N') || (s[i] == 'n')) return true;
    }
    return false;
  }

  inline std::string
  compressStr(std::string const& data) {
    std::stringstream compressed;
    std::stringstream origin(data);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> out;
    out.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_speed)));
    out.push(origin);
    boost::iostreams::copy(out, compressed);
    return compressed.str();
  }
  
  inline double
  entropy(std::string const& st) {
    typedef double TPrecision;
    std::vector<char> stvec(st.begin(), st.end());
    std::set<char> alphabet(stvec.begin(), stvec.end());
    TPrecision ent = 0;
    for(std::set<char>::const_iterator c = alphabet.begin(); c != alphabet.end(); ++c) {
      int ctr = 0;
      for (std::vector<char>::const_iterator s = stvec.begin(); s != stvec.end(); ++s)
	if (*s == *c) ++ctr;
      TPrecision freq = (TPrecision) ctr / (TPrecision) stvec.size();
      ent += (freq) * log(freq)/log(2);
    }
    return -ent;
  }
  
  template<typename TSignalMatrix>
  inline void
  expandpiecewiseconstant(std::vector<uint32_t> const& jumps, TSignalMatrix const& val, TSignalMatrix& res) {
    uint32_t ncol = val.shape()[1];
    
    std::vector<int32_t> b(jumps.begin(), jumps.end());
    b.push_back(-1);
    std::sort(b.begin(), b.end());

    res.resize(boost::extents[b[b.size()-1]+1][ncol]);
    for(uint32_t i = 0; i<jumps.size(); ++i) {
      uint32_t istart = b[i] + 1;
      uint32_t iend = b[i+1] + 1;
      for(uint32_t j = 0; j < ncol; ++j)
	for(uint32_t ki = istart; ki<iend; ++ki) res[ki][j] = val[i][j];
    }
  }

  template<typename TVector>
  inline double
  _medianMutVector(TVector& n) {
    if (!n.empty()) {
      std::size_t s = n.size();
      std::sort(n.begin(), n.end());
      if (s % 2 == 0) return (n[s/2 - 1] + n[s/2]) / 2;
      else return n[s/2];
    } else return 0;
  }
      
  template<typename TBamRecord>
  inline uint8_t
  getLayout(TBamRecord const& al) {
    if (al.flag & BAM_FREAD1) {
      if (!(al.flag & BAM_FREVERSE)) {
	if (!(al.flag & BAM_FMREVERSE)) return 0;
	else return (al.pos < al.mpos) ? 2 : 3;
      } else {
	if (!(al.flag & BAM_FMREVERSE)) return (al.pos > al.mpos) ? 2 : 3;
	else return 1;
      }
    } else {
      if (!(al.flag & BAM_FREVERSE)) {
	if (!(al.flag & BAM_FMREVERSE)) return 0;
	else return (al.pos < al.mpos) ? 2 : 3;
      } else {
	if (!(al.flag & BAM_FMREVERSE)) return (al.pos > al.mpos) ? 2 : 3;
	else return 1;
      }
    }
  }

  inline uint32_t
  setMinChrLen(bam_hdr_t const* hdr, double const xx) {
    uint32_t minChrLen = 0;
    std::vector<uint32_t> chrlen(hdr->n_targets, 0);
    uint64_t genomelen = 0;
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      chrlen[refIndex] = hdr->target_len[refIndex];
      genomelen += hdr->target_len[refIndex];
    }
    std::sort(chrlen.begin(), chrlen.end(), std::greater<uint32_t>());
    uint64_t cumsum = 0;
    for(uint32_t i = 0; i < chrlen.size(); ++i) {
      cumsum += chrlen[i];
      minChrLen = chrlen[i];
      if (cumsum > genomelen * xx) break;
    }
    return minChrLen;
  }
  
  
  template<typename TConfig>
  inline bool
  chrNoData(TConfig const& c, uint32_t const refIndex, hts_idx_t const* idx) {
    // Check we have mapped reads on this chromosome
    std::string suffix("cram");
    std::string str(c.bamFile.string());
    if ((str.size() >= suffix.size()) && (str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0)) return false;
    uint64_t mapped = 0;
    uint64_t unmapped = 0;
    hts_idx_get_stat(idx, refIndex, &mapped, &unmapped);
    if (mapped) return false;
    else return true;
  }    
  
  inline unsigned hash_string(const char *s) {
    unsigned h = 37;
    while (*s) {
      h = (h * 54059) ^ (s[0] * 76963);
      s++;
    }
    return h;
  }

  inline std::size_t hash_se(bam1_t* rec) {
    std::size_t seed = hash_string(bam_get_qname(rec));
    boost::hash_combine(seed, rec->core.tid);
    boost::hash_combine(seed, rec->core.pos);
    return seed;
  }
  
  inline std::size_t hash_pair(bam1_t* rec) {
    std::size_t seed = hash_string(bam_get_qname(rec));
    boost::hash_combine(seed, rec->core.tid);
    boost::hash_combine(seed, rec->core.pos);
    boost::hash_combine(seed, rec->core.mtid);
    boost::hash_combine(seed, rec->core.mpos);
    return seed;
  }

  inline std::size_t hash_pair_mate(bam1_t* rec) {
    std::size_t seed = hash_string(bam_get_qname(rec));
    boost::hash_combine(seed, rec->core.mtid);
    boost::hash_combine(seed, rec->core.mpos);
    boost::hash_combine(seed, rec->core.tid);
    boost::hash_combine(seed, rec->core.pos);
    return seed;
  }


  inline uint32_t alignmentLength(bam1_t const* rec) {
    uint32_t* cigar = bam_get_cigar(rec);
    uint32_t alen = 0;
    for (std::size_t i = 0; i < rec->core.n_cigar; ++i)
      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF) || (bam_cigar_op(cigar[i]) == BAM_CDEL) || (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP)) alen += bam_cigar_oplen(cigar[i]);
    return alen;
  }

  inline uint32_t
  lastAlignedPosition(bam1_t const* rec) {
    return rec->core.pos + alignmentLength(rec);
  }

  inline uint32_t halfAlignmentLength(bam1_t const* rec) {
    return (alignmentLength(rec) / 2);
  }

  template<typename TConfig>
  inline void
  getLibraryParams(TConfig const& c, LibraryInfo& li) {
    // Open file handles
    samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.bamFile.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Iterate all samples
    uint32_t maxAlignmentsScreened=10000000;
    uint32_t maxNumAlignments=1000000;
    uint32_t minNumAlignments=1000;
    uint32_t alignmentCount=0;
    uint32_t processedNumPairs = 0;
    uint32_t processedNumReads = 0;
    uint32_t rplus = 0;
    uint32_t nonrplus = 0;
    typedef std::vector<uint32_t> TSizeVector;
    TSizeVector vecISize;
    TSizeVector readSize;

    // Estimate lib params
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Estimate library parameters" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );
    
    // Collect insert sizes
    bool libCharacterized = false;
    for(uint32_t refIndex=0; refIndex < (uint32_t) hdr->n_targets; ++refIndex) {
      ++show_progress;
      if (libCharacterized) continue;
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (!(rec->core.flag & BAM_FREAD2) && (rec->core.l_qseq < 65000)) {
	  if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	  if ((alignmentCount > maxAlignmentsScreened) || ((processedNumReads >= maxNumAlignments) && (processedNumPairs == 0)) || (processedNumPairs >= maxNumAlignments)) {
	    // Paired-end library with enough pairs
	    libCharacterized = true;
	    break;
	  }
	  ++alignmentCount;
	      
	  // Single-end
	  if (processedNumReads < maxNumAlignments) {
	    readSize.push_back(rec->core.l_qseq);
	    ++processedNumReads;
	  }
	  
	  // Paired-end
	  if ((rec->core.flag & BAM_FPAIRED) && !(rec->core.flag & BAM_FMUNMAP) && (rec->core.tid==rec->core.mtid)) {
	    if (processedNumPairs < maxNumAlignments) {
	      vecISize.push_back(abs(rec->core.isize));
	      if (getLayout(rec->core) == 2) ++rplus;
	      else ++nonrplus;
	      ++processedNumPairs;
	    }
	  }
	}
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }
    
    // Get library parameters
    if (processedNumReads >= minNumAlignments) {
      std::sort(readSize.begin(), readSize.end());
      li.rs = readSize[readSize.size() / 2];
    }
    if (processedNumPairs >= minNumAlignments) {
      std::sort(vecISize.begin(), vecISize.end());
      int32_t median = vecISize[vecISize.size() / 2];
      std::vector<uint32_t> absDev;
      for(uint32_t i = 0; i < vecISize.size(); ++i) absDev.push_back(std::abs((int32_t) vecISize[i] - median));
      std::sort(absDev.begin(), absDev.end());
      int32_t mad = absDev[absDev.size() / 2];

      // Get default library orientation
      if ((median >= 50) && (median<=100000)) {
	if (rplus < nonrplus) {
	  std::cerr << "Warning: Sample has a non-default paired-end layout!" << std::endl;
	  std::cerr << "The expected paired-end orientation is   ---Read1--->      <---Read2---  which is the default illumina paired-end layout." << std::endl;
	    
	}
	li.median = median;
	li.mad = mad;
	li.maxNormalISize = median + (c.mad * mad);
	li.minNormalISize = 0;
	if (c.mad * mad < median) li.minNormalISize = median - (c.mad * mad);
      }
    }

    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
  }

}

#endif
