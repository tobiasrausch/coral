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

#ifndef SCAN_H
#define SCAN_H

#include <limits>

#include <boost/icl/split_interval_map.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/progress.hpp>

#include <htslib/sam.h>

#include "version.h"
#include "util.h"


namespace coralns
{

  template<typename TConfig>
  inline int32_t
  scan(TConfig const& c, std::vector< std::vector<uint32_t> >& scanCounts) {

    // Load bam file
    samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.bamFile.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Parse BAM file
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Scanning Windows" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );

    // Iterate chromosomes
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      ++show_progress;
      if (chrNoData(c, refIndex, idx)) continue;
      // Exclude small chromosomes
      if (hdr->target_len[refIndex] < c.minChrLen) continue;
      // Exclude sex chromosomes
      if ((std::string(hdr->target_name[refIndex]) == "chrX") || (std::string(hdr->target_name[refIndex]) == "chrY") || (std::string(hdr->target_name[refIndex]) == "X") || (std::string(hdr->target_name[refIndex]) == "Y")) continue;

      // Bins on this chromosome
      uint32_t allbins = hdr->target_len[refIndex] / c.scanWindow;
      scanCounts[refIndex].resize(allbins, 0);
      
      // Mate map
      typedef boost::unordered_map<std::size_t, bool> TMateMap;
      TMateMap mateMap;

      // Count reads
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      int32_t lastAlignedPos = 0;
      std::set<std::size_t> lastAlignedPosReads;
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	if ((rec->core.flag & BAM_FPAIRED) && ((rec->core.flag & BAM_FMUNMAP) || (rec->core.tid != rec->core.mtid))) continue;
	if (rec->core.qual < c.minQual) continue;

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

	  // Isize not accurate so use a best guess
	  if (rec->core.flag & BAM_FREVERSE) midPoint = rec->core.pos + alignmentLength(rec) - (c.meanisize / 2);
	  else midPoint = rec->core.pos + (c.meanisize / 2);
	}

	// Count fragment
	if ((midPoint >= 0) && (midPoint < (int32_t) hdr->target_len[refIndex])) {
	  uint32_t bin = midPoint / c.scanWindow;
	  if (bin < allbins) ++scanCounts[refIndex][bin];
	}

	/*
	// Is this a good position to sample GC?
	if ((li.minNormalISize < isize) && (isize < li.maxNormalISize)) {
	  if ((rec->core.flag & BAM_FPAIRED) && (getLayout(rec->core) == 2)) {
	    // Check alignment quality
	    if ((!c.alignmentQ) || (getPercentIdentity(rec, ref) >= c.alignmentQ)) {
	      if (uniqContent[midPoint] == c.meanisize) {
		// Valid bin?
		uint32_t bin = midPoint / c.scanWindow;
		if (bin < scanCounts[refIndex].size()) {
		  if ((scanCounts[refIndex][bin] > cb.first) && (scanCounts[refIndex][bin] < cb.second)) {
		    if (gcpos[midPoint] < maxCoverage - 1) ++gcpos[midPoint];
		  }
		}
	      }
	    }
	  }
	}
	*/
      }
      // Clean-up
      bam_destroy1(rec);
      hts_itr_destroy(iter);
      mateMap.clear();
    }
	  
    // clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
    
    return 0;
  }

  
}

#endif
