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

#ifndef GCBIAS_H
#define GCBIAS_H

#include <boost/unordered_map.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <boost/dynamic_bitset.hpp>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "util.h"

namespace coralns
{
  struct GcBias {
    int32_t sample;
    int32_t reference;
    double fractionSample;
    double fractionReference;
    double percentileSample;
    double percentileReference;
    double obsexp;

    GcBias() : sample(0), reference(0), fractionSample(0), fractionReference(0), percentileSample(0), percentileReference(0), obsexp(0) {}
  };

  template<typename TConfig>
  inline std::pair<uint32_t, uint32_t>
  gcBound(TConfig const& c, std::vector<GcBias>& gcbias) {
    uint32_t lowerBound = 0;
    uint32_t upperBound = gcbias.size();
    for(uint32_t i = 0; i < gcbias.size(); ++i) {
      if ((gcbias[i].percentileSample < c.exclgc) || (gcbias[i].percentileReference < c.exclgc)) lowerBound = i;
      if ((gcbias[i].percentileSample + c.exclgc > 1) || (gcbias[i].percentileReference + c.exclgc > 1)) {
	if (i < upperBound) upperBound = i;
      }
    }
    if (lowerBound >= upperBound) upperBound = lowerBound + 1;
    /*
    // Adjust total
    uint64_t totalSampleCount = 0;
    uint64_t totalReferenceCount = 0;
    for(uint32_t i = lowerBound + 1; i < upperBound; ++i) {
      totalSampleCount += gcbias[i].sample;
      totalReferenceCount += gcbias[i].reference;
    }    
    // Re-estimate observed/expected
    for(uint32_t i = lowerBound + 1; i < upperBound; ++i) {
      gcbias[i].fractionSample = (double) gcbias[i].sample / (double) totalSampleCount;
      gcbias[i].fractionReference = (double) gcbias[i].reference / (double) totalReferenceCount;
      gcbias[i].obsexp = 1;
      if (gcbias[i].fractionReference > 0) gcbias[i].obsexp = gcbias[i].fractionSample / gcbias[i].fractionReference;
    }
    */
    
    return std::make_pair(lowerBound, upperBound);
  }

  template<typename TConfig>
  inline void
  gcBias(TConfig const& c, std::vector< std::vector<uint32_t> > const& scanCounts, LibraryInfo const& li, std::vector<GcBias>& gcbias) {

    typedef std::pair<uint32_t, uint32_t> TCountBounds;
    TCountBounds cb = estCountBounds(scanCounts);
    
    // Load bam file
    samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.bamFile.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Parse bam (contig by contig)
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Estimate GC bias" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );

    faidx_t* faiMap = fai_load(c.mapFile.string().c_str());
    faidx_t* faiRef = fai_load(c.genome.string().c_str());
    for (int refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      ++show_progress;
      if (scanCounts[refIndex].empty()) continue;

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

      // Mappability
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet uniq(hdr->target_len[refIndex], false);
      for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	if (seq[i] == 'C') uniq[i] = 1;
      }

      // Get GC- and N-content
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet gcref(hdr->target_len[refIndex], false);
      for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	if ((ref[i] == 'c') || (ref[i] == 'C') || (ref[i] == 'g') || (ref[i] == 'G')) gcref[i] = 1;
      }

      // Reference GC
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
	//std::string refslice = boost::to_upper_copy(std::string(ref + pos - halfwin, ref + pos + halfwin + 1));
	//std::cerr << refslice << ',' << c.meanisize << ',' << gcsum << ',' << usum << std::endl;
	if (usum == c.meanisize) {
	  // Valid bin?
	  uint32_t bin = pos / c.scanWindow;
	  if (bin < scanCounts[refIndex].size()) {
	    if ((scanCounts[refIndex][bin] > cb.first) && (scanCounts[refIndex][bin] < cb.second)) ++gcbias[gcsum].reference;
	  }
	}
      }

      // Mate map
      typedef boost::unordered_map<std::size_t, bool> TMateMap;
      TMateMap mateMap;
      
      // Parse BAM
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      int32_t lastAlignedPos = 0;
      std::set<std::size_t> lastAlignedPosReads;
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	if ((rec->core.flag & BAM_FPAIRED) && ((rec->core.flag & BAM_FMUNMAP) || (rec->core.tid != rec->core.mtid))) continue;
	if ((rec->core.flag & BAM_FPAIRED) && (getLayout(rec->core) != 2)) continue;
	if (rec->core.qual < c.minQual) continue;

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
	    mateMap[hv]= true;
	    continue;
	  } else {
	    // Second read
	    std::size_t hv = hash_pair_mate(rec);
	    if ((mateMap.find(hv) == mateMap.end()) || (!mateMap[hv])) continue; // Mate discarded
	    mateMap[hv] = false;
	  }
	
	  // Insert size filter
	  int32_t isize = (rec->core.pos + alignmentLength(rec)) - rec->core.mpos;
	  if ((isize <= li.minNormalISize) || (isize >= li.maxNormalISize)) continue;

	  // Count fragment mid-points
	  int32_t midPoint = rec->core.mpos + (int32_t) (isize/2);
	  int32_t fragstart = midPoint - halfwin;
	  int32_t fragend = midPoint + halfwin + 1;
	  if ((fragstart >= 0) && (fragend < (int32_t) hdr->target_len[refIndex])) {
	    int32_t usum = 0;
	    int32_t gcsum = 0;
	    for(int32_t i = fragstart; i < fragend; ++i) {
	      if (uniq[i]) ++usum;
	      if (gcref[i]) ++gcsum;
	    }
	    //std::string refslice = boost::to_upper_copy(std::string(ref + fragstart, ref + fragend));
	    //std::cerr << refslice << ',' << c.meanisize << ',' << gcsum << ',' << usum << std::endl;
	    if (usum == c.meanisize) {
	      // Valid bin?
	      uint32_t bin = midPoint / c.scanWindow;
	      if (bin < scanCounts[refIndex].size()) {
		if ((scanCounts[refIndex][bin] > cb.first) && (scanCounts[refIndex][bin] < cb.second)) ++gcbias[gcsum].sample;
	      }
	    }
	  }
	}
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);
      if (seq != NULL) free(seq);
      if (ref != NULL) free(ref);
    }

    // Determine percentiles
    uint64_t totalSampleCount = 0;
    uint64_t totalReferenceCount = 0;
    for(uint32_t i = 0; i < gcbias.size(); ++i) {
      totalSampleCount += gcbias[i].sample;
      totalReferenceCount += gcbias[i].reference;
    }
    uint64_t cumSample = 0;
    uint64_t cumReference = 0;
    for(uint32_t i = 0; i < gcbias.size(); ++i) {
      cumSample += gcbias[i].sample;
      cumReference += gcbias[i].reference;
      gcbias[i].fractionSample = (double) gcbias[i].sample / (double) totalSampleCount;
      gcbias[i].fractionReference = (double) gcbias[i].reference / (double) totalReferenceCount;
      gcbias[i].percentileSample = (double) cumSample / (double) totalSampleCount;
      gcbias[i].percentileReference = (double) cumReference / (double) totalReferenceCount;
      gcbias[i].obsexp = 1;
      if (gcbias[i].fractionReference > 0) gcbias[i].obsexp = gcbias[i].fractionSample / gcbias[i].fractionReference;
    }
      
    // Output GC-bias
    //for(uint32_t i = 0; i < gcbias.size(); ++i) std::cerr << i << "\t(" << gcbias[i].sample << "," << gcbias[i].reference << ")\t(" << gcbias[i].percentileSample << "," << gcbias[i].percentileReference << ")\t(" << gcbias[i].fractionSample << "," << gcbias[i].fractionReference << ")\t" << gcbias[i].obsexp << std::endl;

    fai_destroy(faiRef);
    fai_destroy(faiMap);
    hts_idx_destroy(idx);
    sam_close(samfile);
    bam_hdr_destroy(hdr);
  }

}

#endif
