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

#ifndef COUNT_H
#define COUNT_H

#include <limits>

#include <boost/icl/split_interval_map.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/progress.hpp>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "baf.h"
#include "cnv.h"
#include "bed.h"
#include "version.h"
#include "segment.h"
#include "scan.h"
#include "util.h"
#include "gmm.h"

namespace coralns
{

  struct CountDNAConfig {
    bool hasStatsFile;
    bool hasBedFile;
    bool hasScanFile;
    bool hasControlFile;
    bool noScanWindowSelection;
    uint32_t nchr;
    uint32_t meanisize;
    uint32_t window_size;
    uint32_t window_offset;
    uint32_t scanWindow;
    uint32_t minChrLen;
    uint32_t minCnvSize;
    uint16_t minQual;
    uint16_t minBaseQual;
    uint16_t mad;
    uint16_t minCoverage;
    uint16_t minSnps;
    uint16_t ploidy;
    float exclgc;
    float uniqueToTotalCovRatio;
    float fracWindow;
    float fragmentUnique;
    std::string sampleName;
    std::string outprefix;
    boost::filesystem::path genome;
    boost::filesystem::path statsFile;
    boost::filesystem::path mapFile;
    boost::filesystem::path controlFile;
    boost::filesystem::path bamFile;
    boost::filesystem::path bedFile;
    boost::filesystem::path vcffile;
    boost::filesystem::path scanFile;
  };


  
  template<typename TConfig, typename TGenomicVariants>
  inline int32_t
  bamCount(TConfig const& c, LibraryInfo const& li, std::vector<GcBias> const& gcbias, std::pair<uint32_t, uint32_t> const& gcbound, TGenomicVariants const& cvar, TGenomicVariants const& gvar) {
    // Load bam file
    samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.bamFile.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // BED regions
    typedef std::set<std::pair<uint32_t, uint32_t> > TChrIntervals;
    typedef std::vector<TChrIntervals> TRegionsGenome;
    TRegionsGenome bedRegions;
    if (c.hasBedFile) {
      if (!_parsePotOverlappingIntervals(c.bedFile.string(), c.hasBedFile, hdr, bedRegions)) {
	std::cerr << "Couldn't parse BED intervals. Do the chromosome names match?" << std::endl;
	return 1;
      }
    }

    // Estimate SD
    SDAggregator sda(c.minCnvSize);
    
    // Parse BAM file
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Count fragments" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );

    // Open output files
    std::string filename = c.outprefix + ".fixed.cov.gz";
    boost::iostreams::filtering_ostream dataOutFixed;
    dataOutFixed.push(boost::iostreams::gzip_compressor());
    dataOutFixed.push(boost::iostreams::file_sink(filename.c_str(), std::ios_base::out | std::ios_base::binary));
    dataOutFixed << "chr\tstart\tend\t" << c.sampleName << "_mappable\t" << c.sampleName << "_counts\t" << c.sampleName << "_CN\t" << c.sampleName << "_MAF" << std::endl;
    filename = c.outprefix + ".adaptive.cov.gz";
    boost::iostreams::filtering_ostream dataOutAdapt;
    dataOutAdapt.push(boost::iostreams::gzip_compressor());
    dataOutAdapt.push(boost::iostreams::file_sink(filename.c_str(), std::ios_base::out | std::ios_base::binary));
    dataOutAdapt << "chr\tstart\tend\t" << c.sampleName << "_mappable\t" << c.sampleName << "_counts\t" << c.sampleName << "_CN\t" << c.sampleName << "_MAF" << std::endl;
    filename = c.outprefix + ".baf.gz";
    boost::iostreams::filtering_ostream dataOutBaf;
    dataOutBaf.push(boost::iostreams::gzip_compressor());
    dataOutBaf.push(boost::iostreams::file_sink(filename.c_str(), std::ios_base::out | std::ios_base::binary));
    dataOutBaf << "chr\tpos\t" << c.sampleName << "_control_baf\t" << c.sampleName << "_target_baf" << std::endl;
    
    
    // Iterate chromosomes
    faidx_t* faiMap = fai_load(c.mapFile.string().c_str());
    faidx_t* faiRef = fai_load(c.genome.string().c_str());
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      ++show_progress;
      if (chrNoData(c, refIndex, idx)) continue;

      // Check presence in mappability map
      std::string tname(hdr->target_name[refIndex]);
      int32_t seqlen = faidx_seq_len(faiMap, tname.c_str());
      if (seqlen == - 1) continue;
      else seqlen = -1;
      char* seq = faidx_fetch_seq(faiMap, tname.c_str(), 0, faidx_seq_len(faiMap, tname.c_str()), &seqlen);

      // Check presence in reference
      seqlen = faidx_seq_len(faiRef, tname.c_str());
      if (seqlen == - 1) continue;
      else seqlen = -1;
      char* ref = faidx_fetch_seq(faiRef, tname.c_str(), 0, faidx_seq_len(faiRef, tname.c_str()), &seqlen);

      // Get GC and Mappability
      std::vector<uint16_t> uniqContent(hdr->target_len[refIndex], 0);
      std::vector<uint16_t> gcContent(hdr->target_len[refIndex], 0);
      {
	// Mappability map
	typedef boost::dynamic_bitset<> TBitSet;
	TBitSet uniq(hdr->target_len[refIndex], false);
	for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	  if (seq[i] == 'C') uniq[i] = 1;
	}

	// GC map
	typedef boost::dynamic_bitset<> TBitSet;
	TBitSet gcref(hdr->target_len[refIndex], false);
	for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	  if ((ref[i] == 'c') || (ref[i] == 'C') || (ref[i] == 'g') || (ref[i] == 'G')) gcref[i] = 1;
	}

	// Sum across fragment
	int32_t halfwin = (int32_t) (c.meanisize / 2);
	int32_t usum = 0;
	int32_t gcsum = 0;
	for(int32_t pos = halfwin; pos < (int32_t) hdr->target_len[refIndex] - halfwin; ++pos) {
	  if (pos == halfwin) {
	    for(int32_t i = pos - halfwin; i<=pos+halfwin; ++i) {
	      usum += uniq[i];
	      gcsum += gcref[i];
	    }
	  } else {
	    usum -= uniq[pos - halfwin - 1];
	    gcsum -= gcref[pos - halfwin - 1];
	    usum += uniq[pos + halfwin];
	    gcsum += gcref[pos + halfwin];
	  }
	  gcContent[pos] = gcsum;
	  uniqContent[pos] = usum;
	}
      }
      
      // Coverage track
      typedef uint16_t TCount;
      uint32_t maxCoverage = std::numeric_limits<TCount>::max();
      typedef std::vector<TCount> TCoverage;
      TCoverage cov(hdr->target_len[refIndex], 0);

      // Split-read breakpoints
      typedef std::pair<uint32_t, uint32_t> TPosSupport;
      std::vector<TPosSupport> splitBp;
      {
	// Split-read map
	TCoverage splitCovLeft(hdr->target_len[refIndex], 0);
	TCoverage splitCovRight(hdr->target_len[refIndex], 0);
      
	// Mate map
	typedef boost::unordered_map<std::size_t, bool> TMateMap;
	TMateMap mateMap;
	
	// Count reads
	hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	int32_t lastAlignedPos = 0;
	std::set<std::size_t> lastAlignedPosReads;
	while (sam_itr_next(samfile, iter, rec) >= 0) {
	  if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP)) continue;
	  if (rec->core.qual < c.minQual) continue;
	  
	  // Get clippings
	  uint32_t rp = rec->core.pos; // reference pointer
	  uint32_t sp = 0; // sequence pointer
	  uint32_t* cigar = bam_get_cigar(rec);
	  for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	    if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	      rp += bam_cigar_oplen(cigar[i]);
	      sp += bam_cigar_oplen(cigar[i]);
	    } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	      if (bam_cigar_oplen(cigar[i]) > c.minCnvSize) {
		if ((rp >= 0) && (rp < hdr->target_len[refIndex]) && (splitCovRight[rp] < maxCoverage - 1)) ++splitCovRight[rp];
	      }
	      rp += bam_cigar_oplen(cigar[i]);
	      if (bam_cigar_oplen(cigar[i]) > c.minCnvSize) {
		if ((rp >= 0) && (rp < hdr->target_len[refIndex]) && (splitCovLeft[rp] < maxCoverage - 1)) ++splitCovLeft[rp];
	      }
	    } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	      sp += bam_cigar_oplen(cigar[i]);
	    } else if ((bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP)) {
	      bool scleft = false;
	      if (sp == 0) scleft = true;
	      sp += bam_cigar_oplen(cigar[i]);
	      // Min clipping length of 10bp
	      if (bam_cigar_oplen(cigar[i]) > 10) {
		if ((rp >= 0) && (rp < hdr->target_len[refIndex])) {
		  if (scleft) {
		    if (splitCovLeft[rp] < maxCoverage - 1) ++splitCovLeft[rp];
		  } else {
		    if (splitCovRight[rp] < maxCoverage - 1) ++splitCovRight[rp];
		  }
		}
	      }
	    } else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	      rp += bam_cigar_oplen(cigar[i]);
	    } else std::cerr << "Warning: Unknown Cigar operation!" << std::endl;
	  }
	  
	  // Exclude secondary/supplementary alignments
	  if (rec->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
	  if ((rec->core.flag & BAM_FPAIRED) && ((rec->core.flag & BAM_FMUNMAP) || (rec->core.tid != rec->core.mtid))) continue;

	  int32_t midPoint = rec->core.pos + halfAlignmentLength(rec);
	  if (rec->core.flag & BAM_FPAIRED) {
	    // Clean-up the read store for identical alignment positions
	    if (rec->core.pos > lastAlignedPos) {
	      lastAlignedPosReads.clear();
	      lastAlignedPos = rec->core.pos;
	    }
	    
	    if ((rec->core.pos < rec->core.mpos) || ((rec->core.pos == rec->core.mpos) && (lastAlignedPosReads.find(hash_string(bam_get_qname(rec))) == lastAlignedPosReads.end()))) {
	      // First read
	      lastAlignedPosReads.insert(hash_string(bam_get_qname(rec)));
	      std::size_t hv = hash_pair(rec);
	      mateMap[hv] = true;
	      continue;
	    } else {
	      // Second read
	      std::size_t hv = hash_pair_mate(rec);
	      if ((mateMap.find(hv) == mateMap.end()) || (!mateMap[hv])) continue; // Mate discarded
	      mateMap[hv] = false;
	    }
	    
	    // update midpoint
	    int32_t isize = (rec->core.pos + alignmentLength(rec)) - rec->core.mpos;
	    if ((li.minNormalISize < isize) && (isize < li.maxNormalISize)) midPoint = rec->core.mpos + (int32_t) (isize/2);
	  }
	  
	  // Count fragment
	  if ((midPoint >= 0) && (midPoint < (int32_t) hdr->target_len[refIndex]) && (cov[midPoint] < maxCoverage - 1)) ++cov[midPoint];
	}
	// Clean-up
	if (seq != NULL) free(seq);
	if (ref != NULL) free(ref);
	bam_destroy1(rec);
	hts_itr_destroy(iter);

	// Collect split-read breakpoints
	_collectSplitBp(splitCovLeft, splitCovRight, splitBp, c.minCnvSize);
      }

      // Output BAF
      for(uint32_t i = 0; i < cvar[refIndex].size(); ++i) {
	if (gvar[refIndex][i].ref + gvar[refIndex][i].alt > 0) {
	  double baf_target = (double) gvar[refIndex][i].alt /	(double) (gvar[refIndex][i].ref + gvar[refIndex][i].alt);
	  double baf_control = 0.5;
	  if (c.hasControlFile) {
	    if (cvar[refIndex][i].ref + cvar[refIndex][i].alt > 0) {
	      baf_control = (double) cvar[refIndex][i].alt / (double) (cvar[refIndex][i].ref + cvar[refIndex][i].alt);
	    }
	  }
	  if (cvar[refIndex][i].pos != gvar[refIndex][i].pos) std::cerr << "Warning: Variant calling positions differ!" << std::endl;
	  dataOutBaf << std::string(hdr->target_name[refIndex]) << "\t" << cvar[refIndex][i].pos << "\t" << baf_control << "\t" << baf_target << std::endl;
	}
      }

      // Estimate SDs
      uint32_t widx = 0;
      for(uint32_t ws = sda.wsinit; ws < 1000000; ws = ws * 2) {
	for(uint32_t start = 0; start < hdr->target_len[refIndex]; start = start + ws) {
	  if (start + ws < hdr->target_len[refIndex]) {
	    double covsum = 0;
	    double expcov = 0;
	    for(uint32_t pos = start; pos < start + ws; ++pos) {
	      if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
		covsum += cov[pos];
		expcov += gcbias[gcContent[pos]].coverage;
	      }
	    }
	    if (expcov > 0) {
	      double cn = c.ploidy * covsum / expcov;
	      sda.cnSUM[widx] += std::abs(cn - c.ploidy) * std::abs(cn - c.ploidy);
	      ++sda.cnCount[widx];
	    }
	  }
	}
	++widx;
      }
      
      // Call & genotype CNVs
      std::vector<CNV> cnvCalls;
      callCNVs(c, gcbound, gcContent, uniqContent, gcbias, cov, hdr, refIndex, cnvCalls);
      breakpointRefinement(splitBp, cnvCalls);
      std::sort(cnvCalls.begin(), cnvCalls.end(), SortCNVs<CNV>());
      genotypeCNVs(c, sda, cnvCalls);
      
      // BED File (target intervals)
      if (c.hasBedFile) {
	// Adaptive Window Length
	{
	  // Merge overlapping BED entries
	  TChrIntervals citv;
	  _mergeOverlappingBedEntries(bedRegions[refIndex], citv);

	  // Tile merged intervals
	  double covsum = 0;
	  double expcov = 0;
	  double obsexp = 0;
	  uint32_t winlen = 0;
	  uint32_t start = 0;
	  bool endOfWindow = true;
	  typename TChrIntervals::iterator it = citv.begin();
	  if (it != citv.end()) start = it->first;
	  while(endOfWindow) {
	    endOfWindow = false;
	    for(it = citv.begin(); ((it != citv.end()) && (!endOfWindow)); ++it) {
	      if ((it->first < it->second) && (it->second <= hdr->target_len[refIndex])) {
		if (start >= it->second) {
		  if (start == it->second) {
		    // Special case
		    typename TChrIntervals::iterator itNext = it;
		    ++itNext;
		    if (itNext != citv.end()) start = itNext->first;
		  }
		  continue;
		}
		for(uint32_t pos = it->first; ((pos < it->second) && (!endOfWindow)); ++pos) {
		  if (pos < start) continue;
		  if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
		    covsum += cov[pos];
		    obsexp += gcbias[gcContent[pos]].obsexp;
		    expcov += gcbias[gcContent[pos]].coverage;
		    ++winlen;
		    if (winlen == c.window_size) {
		      obsexp /= (double) winlen;
		      double count = ((double) covsum / obsexp ) * (double) c.window_size / (double) winlen;
		      double cn = c.ploidy * covsum / expcov;
		      dataOutAdapt << std::string(hdr->target_name[refIndex]) << "\t" << start << "\t" << (pos + 1) << "\t" << winlen << "\t" << count << "\t" << cn;
		      double maf = mafSegment(c, start, pos + 1, cvar[refIndex], gvar[refIndex]);
		      if (maf != -1) dataOutAdapt << "\t" << maf << std::endl;
		      else dataOutAdapt << "\tNA" << std::endl;
		      // reset
		      covsum = 0;
		      expcov = 0;
		      obsexp = 0;
		      winlen = 0;
		      if (c.window_offset == c.window_size) {
			// Move on
			start = pos + 1;
			endOfWindow = true;
		      } else {
			// Rewind
			for(typename TChrIntervals::iterator sit = citv.begin(); ((sit != citv.end()) && (!endOfWindow)); ++sit) {
			  if ((sit->first < sit->second) && (sit->second <= hdr->target_len[refIndex])) {
			    if (start >= sit->second) continue;
			    for(uint32_t k = sit->first; ((k < sit->second) && (!endOfWindow)); ++k) {
			      if (k < start) continue;
			      if ((gcContent[k] > gcbound.first) && (gcContent[k] < gcbound.second) && (uniqContent[k] >= c.fragmentUnique * c.meanisize)) {
				++winlen;
				if (winlen == c.window_offset) {
				  start = k + 1;
				  winlen = 0;
				  endOfWindow = true;
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
	// Fixed Window Length
	{
	  for(typename TChrIntervals::iterator it = bedRegions[refIndex].begin(); it != bedRegions[refIndex].end(); ++it) {
	    if ((it->first < it->second) && (it->second <= hdr->target_len[refIndex])) {
	      double covsum = 0;
	      double expcov = 0;
	      double obsexp = 0;
	      uint32_t winlen = 0;
	      for(uint32_t pos = it->first; pos < it->second; ++pos) {
		if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
		  covsum += cov[pos];
		  obsexp += gcbias[gcContent[pos]].obsexp;
		  expcov += gcbias[gcContent[pos]].coverage;
		  ++winlen;
		}
	      }
	      if (winlen >= c.fracWindow * (it->second - it->first)) {
		obsexp /= (double) winlen;
		double count = ((double) covsum / obsexp ) * (double) (it->second - it->first) / (double) winlen;
		double cn = c.ploidy * covsum / expcov;
		dataOutFixed << std::string(hdr->target_name[refIndex]) << "\t" << it->first << "\t" << it->second << "\t" << winlen << "\t" << count << "\t" << cn;
		double maf = mafSegment(c, it->first, it->second, cvar[refIndex], gvar[refIndex]);
		if (maf != -1) dataOutFixed << "\t" << maf << std::endl;
		else dataOutFixed << "\tNA" << std::endl;
	      } else {
		dataOutFixed << std::string(hdr->target_name[refIndex]) << "\t" << it->first << "\t" << it->second << "\tNA\tNA\tNA\tNA" << std::endl;
	      }
	    }
	  }
	}
      } else {
	// Use adaptive windows (genomic tiling)
	{
	  double covsum = 0;
	  double expcov = 0;
	  double obsexp = 0;
	  uint32_t winlen = 0;
	  uint32_t start = 0;
	  uint32_t pos = 0;
	  while(pos < hdr->target_len[refIndex]) {
	    if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	      covsum += cov[pos];
	      obsexp += gcbias[gcContent[pos]].obsexp;
	      expcov += gcbias[gcContent[pos]].coverage;
	      ++winlen;
	      if (winlen == c.window_size) {
		obsexp /= (double) winlen;
		double count = ((double) covsum / obsexp ) * (double) c.window_size / (double) winlen;
		double cn = c.ploidy * covsum / expcov;
		dataOutAdapt << std::string(hdr->target_name[refIndex]) << "\t" << start << "\t" << (pos + 1) << "\t" << winlen << "\t" << count << "\t" << cn;
		double maf = mafSegment(c, start, pos + 1, cvar[refIndex], gvar[refIndex]);
		if (maf != -1) dataOutAdapt << "\t" << maf << std::endl;
		else dataOutAdapt << "\tNA" << std::endl;
		// reset
		covsum = 0;
		expcov = 0;
		obsexp = 0;
		winlen = 0;
		if (c.window_offset == c.window_size) {
		  // Move on
		  start = pos + 1;
		} else {
		  // Rewind
		  for(uint32_t k = start; k < hdr->target_len[refIndex]; ++k) {
		    if ((gcContent[k] > gcbound.first) && (gcContent[k] < gcbound.second) && (uniqContent[k] >= c.fragmentUnique * c.meanisize)) {
		      ++winlen;
		      if (winlen == c.window_offset) {
			start = k + 1;
			pos = k;
			winlen = 0;
			break;
		      }
		    }
		  }
		}
	      }
	    }
	    ++pos;
	  }
	}
	{
	  // Fixed windows (genomic tiling)
	  for(uint32_t start = 0; start < hdr->target_len[refIndex]; start = start + c.window_offset) {
	    if (start + c.window_size < hdr->target_len[refIndex]) {
	      double covsum = 0;
	      double expcov = 0;
	      double obsexp = 0;
	      uint32_t winlen = 0;
	      for(uint32_t pos = start; pos < start + c.window_size; ++pos) {
		if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
		  covsum += cov[pos];
		  obsexp += gcbias[gcContent[pos]].obsexp;
		  expcov += gcbias[gcContent[pos]].coverage;
		  ++winlen;
		}
	      }
	      if (winlen >= c.fracWindow * c.window_size) {
		obsexp /= (double) winlen;
		double count = ((double) covsum / obsexp ) * (double) c.window_size / (double) winlen;
		double cn = c.ploidy * covsum / expcov;
		dataOutFixed << std::string(hdr->target_name[refIndex]) << "\t" << start << "\t" << (start + c.window_size) << "\t" << winlen << "\t" << count << "\t" << cn;
		double maf = mafSegment(c, start, start + c.window_size, cvar[refIndex], gvar[refIndex]);
		if (maf != -1) dataOutFixed << "\t" << maf << std::endl;
		else dataOutFixed << "\tNA" << std::endl;
	      }
	    }
	  }
	}
      }
    }
	  
    // clean-up
    fai_destroy(faiRef);
    fai_destroy(faiMap);
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
    dataOutFixed.pop();
    dataOutAdapt.pop();
    dataOutBaf.pop();
    return 0;
  }

  
  int countReads(int argc, char **argv) {
    CountDNAConfig c;

    //cov.q10.d3.p0.0005.f0.8.e0.95.k0.5.gz
    //cov.q10.d3.p0.0005.f0.8.e0.97.k0.5.gz

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome file")
      ("sample,s", boost::program_options::value<std::string>(&c.sampleName)->default_value("NA12878"), "sample name")
      ("quality,q", boost::program_options::value<uint16_t>(&c.minQual)->default_value(10), "min. mapping quality")
      ("mappability,m", boost::program_options::value<boost::filesystem::path>(&c.mapFile), "input mappability map")
      ("minsize,z", boost::program_options::value<uint32_t>(&c.minCnvSize)->default_value(500), "min. CNV size")
      ("ploidy,y", boost::program_options::value<uint16_t>(&c.ploidy)->default_value(2), "baseline ploidy")
      ("fragment,e", boost::program_options::value<float>(&c.fragmentUnique)->default_value(0.97), "min. fragment uniqueness [0,1]")
      ("outprefix,o", boost::program_options::value<std::string>(&c.outprefix)->default_value("outprefix"), "output file prefix")
      ;

    boost::program_options::options_description window("Window options");
    window.add_options()
      ("window-size,i", boost::program_options::value<uint32_t>(&c.window_size)->default_value(10000), "window size")
      ("window-offset,j", boost::program_options::value<uint32_t>(&c.window_offset)->default_value(10000), "window offset")
      ("bed-intervals,b", boost::program_options::value<boost::filesystem::path>(&c.bedFile), "input BED file")
      ("fraction-window,k", boost::program_options::value<float>(&c.fracWindow)->default_value(0.25), "min. callable window fraction [0,1]")
      ;

    boost::program_options::options_description gcopt("GC options");
    gcopt.add_options()
      ("scan-window,w", boost::program_options::value<uint32_t>(&c.scanWindow)->default_value(10000), "scanning window size")
      ("fraction-unique,f", boost::program_options::value<float>(&c.uniqueToTotalCovRatio)->default_value(0.8), "uniqueness filter for scan windows [0,1]")
      ("scan-regions,r", boost::program_options::value<boost::filesystem::path>(&c.scanFile), "scanning regions in BED format")
      ("mad-cutoff,d", boost::program_options::value<uint16_t>(&c.mad)->default_value(3), "median + 3 * mad count cutoff")
      ("percentile,p", boost::program_options::value<float>(&c.exclgc)->default_value(0.0005), "excl. extreme GC fraction")
      ("no-window-selection,n", "no scan window selection")
      ;
    
    boost::program_options::options_description vcopt("Variant options");
    vcopt.add_options()
      ("basequality,a", boost::program_options::value<uint16_t>(&c.minBaseQual)->default_value(10), "min. base quality")
      ("snps,x", boost::program_options::value<uint16_t>(&c.minSnps)->default_value(3), "min. #SNPs per segment")
      ("coverage,c", boost::program_options::value<uint16_t>(&c.minCoverage)->default_value(10), "min. SNP coverage")
      ("vcffile,v", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "input VCF file")
      ("controlfile,l", boost::program_options::value<boost::filesystem::path>(&c.controlFile), "control file")
      ;

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("statsfile,t", boost::program_options::value<boost::filesystem::path>(&c.statsFile), "stats output file")
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.bamFile), "input bam file")
      ;

    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);

    // Set the visibility
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(window).add(gcopt).add(vcopt).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(window).add(gcopt).add(vcopt);

    // Parse command-line
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);

    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome")) || (!vm.count("mappability"))) {
      std::cout << std::endl;
      std::cout << "Usage: coral " << argv[0] << " [OPTIONS] -g hg19.fa -m hg19.map <aligned.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    }

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "coral ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;

    // Stats file
    if (vm.count("statsfile")) c.hasStatsFile = true;
    else c.hasStatsFile = false;

    // Control file
    if (vm.count("controlfile")) c.hasControlFile = true;
    else c.hasControlFile = false;
    
    // BED intervals
    if (vm.count("bed-intervals")) c.hasBedFile = true;
    else c.hasBedFile = false;

    // Scan regions
    if (vm.count("scan-regions")) c.hasScanFile = true;
    else c.hasScanFile = false;

    // Scan window selection
    if (vm.count("no-window-selection")) c.noScanWindowSelection = true;
    else c.noScanWindowSelection = false;

    // Check window size
    if (c.window_offset > c.window_size) c.window_offset = c.window_size;
    if (c.window_size == 0) c.window_size = 1;
    if (c.window_offset == 0) c.window_offset = 1;
    
    // Open stats file
    boost::iostreams::filtering_ostream statsOut;
    if (c.hasStatsFile) {
      statsOut.push(boost::iostreams::gzip_compressor());
      statsOut.push(boost::iostreams::file_sink(c.statsFile.string().c_str(), std::ios_base::out | std::ios_base::binary));
    }

    // Library info
    LibraryInfo li;

    // Check bam file
    if (!(boost::filesystem::exists(c.bamFile) && boost::filesystem::is_regular_file(c.bamFile) && boost::filesystem::file_size(c.bamFile))) {
      std::cerr << "Alignment file is missing: " << c.bamFile.string() << std::endl;
      return 1;
    } else {
      samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
      if (samfile == NULL) {
	std::cerr << "Fail to open file " << c.bamFile.string() << std::endl;
	return 1;
      }
      hts_idx_t* idx = sam_index_load(samfile, c.bamFile.string().c_str());
      if (idx == NULL) {
	if (bam_index_build(c.bamFile.string().c_str(), 0) != 0) {
	  std::cerr << "Fail to open index for " << c.bamFile.string() << std::endl;
	  return 1;
	}
      }
      bam_hdr_t* hdr = sam_hdr_read(samfile);
      if (hdr == NULL) {
	std::cerr << "Fail to open header for " << c.bamFile.string() << std::endl;
	return 1;
      }
      c.nchr = hdr->n_targets;
      c.minChrLen = setMinChrLen(hdr, 0.95);

      // Check matching chromosome names
      faidx_t* faiRef = fai_load(c.genome.string().c_str());
      faidx_t* faiMap = fai_load(c.mapFile.string().c_str());
      uint32_t mapFound = 0;
      uint32_t refFound = 0;
      for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) {
	std::string tname(hdr->target_name[refIndex]);
	if (faidx_has_seq(faiMap, tname.c_str())) ++mapFound;
	if (faidx_has_seq(faiRef, tname.c_str())) ++refFound;
	else {
	  std::cerr << "Warning: BAM chromosome " << tname << " not present in reference genome!" << std::endl;
	}
      }
      fai_destroy(faiRef);
      fai_destroy(faiMap);
      if (!mapFound) {
	std::cerr << "Mappability map chromosome naming disagrees with BAM file!" << std::endl;
	return 1;
      }
      if (!refFound) {
	std::cerr << "Reference genome chromosome naming disagrees with BAM file!" << std::endl;
	return 1;
      }
      
      // Estimate insert size
      getLibraryParams(c, li);
      // Fix single-end libraries
      if (!li.median) {
	li.median = 250;
	li.mad = 15;
	li.minNormalISize = 0;
	li.maxNormalISize = 400;
      }
      c.meanisize = ((int32_t) (li.median / 2)) * 2 + 1;
      if (c.hasStatsFile) {
	statsOut << "LP\t" << li.rs << ',' << li.median << ',' << li.mad << ',' << li.minNormalISize << ',' << li.maxNormalISize << std::endl;
      }

      // Clean-up
      bam_hdr_destroy(hdr);
      hts_idx_destroy(idx);
      sam_close(samfile);
    }

    // B-allele frequency
    typedef std::vector<BiallelicSupport> TVariantSupport;
    typedef std::vector<TVariantSupport> TGenomicVariants;
    TGenomicVariants gvar(c.nchr, TVariantSupport());
    TGenomicVariants cvar(c.nchr, TVariantSupport());
    if (vm.count("vcffile")) {
      if (c.hasControlFile) {
	baf(c, true, cvar);  // Het. germline variants
	baf(c, false, cvar, gvar);
      } else {
	// Tumor-only
	baf(c, false, gvar);
      }
    }

    // Scan genomic windows
    typedef std::vector<ScanWindow> TWindowCounts;
    typedef std::vector<TWindowCounts> TGenomicWindowCounts;
    TGenomicWindowCounts scanCounts(c.nchr, TWindowCounts());
    scan(c, li, scanCounts);
    
    // Select stable windows
    selectWindows(c, scanCounts);
    if (c.hasStatsFile) {
      samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
      bam_hdr_t* hdr = sam_hdr_read(samfile);
      statsOut << "SW\tchrom\tstart\tend\tselected\tcoverage\tuniqcov" <<  std::endl;
      for(uint32_t refIndex = 0; refIndex < (uint32_t) hdr->n_targets; ++refIndex) {
	for(uint32_t i = 0; i < scanCounts[refIndex].size(); ++i) {
	  statsOut << "SW\t" <<  hdr->target_name[refIndex] << '\t' << scanCounts[refIndex][i].start << '\t' << scanCounts[refIndex][i].end << '\t' << scanCounts[refIndex][i].select << '\t' << scanCounts[refIndex][i].cov << '\t' << scanCounts[refIndex][i].uniqcov << std::endl;
	}
      }
      bam_hdr_destroy(hdr);
      sam_close(samfile);
    }
    
    // GC bias estimation
    typedef std::pair<uint32_t, uint32_t> TGCBound;
    TGCBound gcbound;
    std::vector<GcBias> gcbias(c.meanisize + 1, GcBias());
    gcBias(c, scanCounts, li, gcbias, gcbound);
    if (c.hasStatsFile) {
      statsOut << "GC\tgcsum\tsample\treference\tpercentileSample\tpercentileReference\tfractionSample\tfractionReference\tobsexp\tmeancoverage" << std::endl;
      for(uint32_t i = 0; i < gcbias.size(); ++i) statsOut << "GC\t" << i << "\t" << gcbias[i].sample << "\t" << gcbias[i].reference << "\t" << gcbias[i].percentileSample << "\t" << gcbias[i].percentileReference << "\t" << gcbias[i].fractionSample << "\t" << gcbias[i].fractionReference << "\t" << gcbias[i].obsexp << "\t" << gcbias[i].coverage << std::endl;
      statsOut << "BoundsGC\t" << gcbound.first << "," << gcbound.second << std::endl;
      statsOut.pop();
    }
    
    // Count reads
    int32_t ret = bamCount(c, li, gcbias, gcbound, cvar, gvar);
    if (ret != 0) return ret;

    // Segment
    /*
    SegmentConfig segc;
    segc.k = 300;
    segc.epsilon = 1e-9;
    segc.dpthreshold = 0.5;
    segc.outfile = boost::filesystem::path(c.outprefix + ".segment.gz");
    segc.signal = boost::filesystem::path(c.outprefix + ".adaptive.cov.gz");
    segmentCovBaf(segc);
    */
    
    // Done
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Done." << std::endl;

    return 0;
  }

  
}

#endif
