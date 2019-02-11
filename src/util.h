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



namespace coralns
{

  struct ScanWindow {
    bool select;
    uint32_t cov;
    uint32_t rplus;
    uint32_t nonrplus;
    uint32_t uniqcov;
    double layoutratio;
  };
    
  struct LibraryInfo {
    int32_t rs;
    int32_t median;
    int32_t mad;
    int32_t minNormalISize;
    int32_t maxNormalISize;

    LibraryInfo() : rs(0), median(0), mad(0), minNormalISize(0), maxNormalISize(0) {}
  };


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
    int32_t median = 0;
    int32_t mad = 0;
    if (processedNumPairs >= minNumAlignments) {
      std::sort(vecISize.begin(), vecISize.end());
      median = vecISize[vecISize.size() / 2];
      std::vector<uint32_t> absDev;
      for(uint32_t i = 0; i < vecISize.size(); ++i) absDev.push_back(std::abs((int32_t) vecISize[i] - median));
      std::sort(absDev.begin(), absDev.end());
      mad = absDev[absDev.size() / 2];

      // Get default library orientation
      if ((median >= 50) && (median<=100000)) {
	if (rplus < nonrplus) {
	  std::cerr << "Warning: Sample has a non-default paired-end layout!" << std::endl;
	  std::cerr << "The expected paired-end orientation is   ---Read1--->      <---Read2---  which is the default illumina paired-end layout." << std::endl;
	    
	} else {
	  li.median = median;
	  li.mad = mad;
	  li.maxNormalISize = median + (c.mad * mad);
	  li.minNormalISize = 0;
	  if (c.mad * mad < median) li.minNormalISize = median - (c.mad * mad);
	}
      }
    }

    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
  }

}

#endif
