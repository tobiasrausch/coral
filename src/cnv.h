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

#ifndef CNV_H
#define CNV_H

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

  template<typename TCoverage>
  inline void
  _collectSplitBp(TCoverage& splitCovLeft, TCoverage& splitCovRight, std::vector<uint32_t>& splitBp, uint32_t const bandwidth) {
    typedef typename TCoverage::value_type TValue;

    // Identify peaks
    TValue defaultMaxVal = 3; // Require at least 3 split-reads
    TValue maxValLeft = defaultMaxVal;
    TValue maxValRight = defaultMaxVal;
    uint32_t maxIdxLeft = 0;
    uint32_t maxIdxRight = 0;

    // Fwd scan
    for(uint32_t i = 0; i < splitCovLeft.size(); ++i) {
      if (i - maxIdxLeft > bandwidth) {
	maxValLeft = defaultMaxVal;
	maxIdxLeft = i;
      }
      if (i - maxIdxRight > bandwidth) {
	maxValRight = defaultMaxVal;
	maxIdxRight = i;
      }
      if (splitCovLeft[i] > maxValLeft) {
	maxValLeft = splitCovLeft[i];
	maxIdxLeft = i;
      } else splitCovLeft[i] = 0;
      if (splitCovRight[i] > maxValRight) {
	maxValRight = splitCovRight[i];
	maxIdxRight = i;
      } else splitCovRight[i] = 0;
    }

    // Rev scan
    maxValLeft = defaultMaxVal;
    maxValRight = defaultMaxVal;
    maxIdxLeft = splitCovLeft.size();
    maxIdxRight = splitCovRight.size();

    // Rev scan
    for(uint32_t i = splitCovLeft.size(); i > 0; --i) {
      if (maxIdxLeft - i > bandwidth) {
	maxValLeft = defaultMaxVal;
	maxIdxLeft = i;
      }
      if (maxIdxRight - i > bandwidth) {
	maxValRight = defaultMaxVal;
	maxIdxRight = i;
      }
      if (splitCovLeft[i-1] > maxValLeft) {
	maxValLeft = splitCovLeft[i-1];
	maxIdxLeft = i;
      } else splitCovLeft[i-1] = 0;
      if (splitCovRight[i-1] > maxValRight) {
	maxValRight = splitCovRight[i-1];
	maxIdxRight = i;
      } else splitCovRight[i-1] = 0;
    }

    // Erase breakpoints in close proximity (small deletions or insertions)
    uint32_t leftSum = 0;
    uint32_t rightSum = 0;
    bool erasePeak = false;
    uint32_t eraseIdx = 0;
    for(uint32_t i = 0; i < splitCovLeft.size(); ++i) {
      if (i < bandwidth) {
	leftSum += splitCovLeft[i];
	rightSum += splitCovRight[i];
      } else {
	if ((leftSum > defaultMaxVal) && (rightSum > defaultMaxVal)) {
	  erasePeak = true;
	  eraseIdx = i;
	}
	leftSum -= splitCovLeft[i - bandwidth];
	rightSum -= splitCovRight[i - bandwidth];
	leftSum += splitCovLeft[i];
	rightSum += splitCovRight[i];	
	if (erasePeak) {
	  splitCovLeft[i - bandwidth] = 0;
	  splitCovRight[i - bandwidth] = 0;
	  if (i - eraseIdx > bandwidth) erasePeak = false;
	}
      }
    }

    // Insert breakpoints
    for(uint32_t i = 0; i < splitCovLeft.size(); ++i) {
      if ((splitCovLeft[i] > 0) || (splitCovRight[i] > 0)) {
	splitBp.push_back(i);
      }
    }
  }


  template<typename TConfig, typename TGcBias, typename TCoverage, typename TGenomicVariants>
  inline void
  callCNVs(TConfig const& c, std::vector<uint32_t> const& splitBp, std::pair<uint32_t, uint32_t> const& gcbound, std::vector<uint16_t> const& gcContent, std::vector<uint16_t> const& uniqContent, TGcBias const& gcbias, TCoverage const& cov, TGenomicVariants const& cvar, TGenomicVariants const& gvar, int32_t const refIndex) {
    for(uint32_t i = 1; i < splitBp.size(); ++i) {
      uint32_t bpLeft = splitBp[i-1];
      uint32_t bpRight = splitBp[i];
      if ((bpLeft < bpRight) && (bpRight - bpLeft >= c.minCnvSize) && (bpRight <= cov.size())) {
	double covsum = 0;
	double expcov = 0;
	double obsexp = 0;
	uint32_t winlen = 0;
	for(uint32_t pos = bpLeft; pos < bpRight; ++pos) {
	  if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	    covsum += cov[pos];
	    obsexp += gcbias[gcContent[pos]].obsexp;
	    expcov += gcbias[gcContent[pos]].coverage;
	    ++winlen;
	  }
	}
	if (winlen >= c.fracWindow * (bpRight - bpLeft)) {
	  obsexp /= (double) winlen;
	  double count = ((double) covsum / obsexp ) * (double) (bpRight - bpLeft) / (double) winlen;
	  double cn = 2 * covsum / expcov;
	  std::cerr << bpLeft << "\t" << bpRight << "\t" << winlen << "\t" << count << "\t" << cn;
	  double maf = mafSegment(c, bpLeft, bpRight, cvar[refIndex], gvar[refIndex]);
	  if (maf != -1) std::cerr << "\t" << maf << std::endl;
	  else std::cerr << "\tNA" << std::endl;
	}
      }
    }
  }
  
  /*
  	// Find #mappable pos
	uint32_t mappable = 0;
	for(uint32_t pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
	  if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	    ++mappable;
	  }
	}
	typedef std::vector<TCount> TCoverage;
	TCoverage mapcov(mappable, 0);
	std::vector<uint16_t> mapGcContent(mappable, 0);
	{
	  uint32_t mapIdx = 0;
	  for(uint32_t pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
	    if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	      mapcov[mapIdx] = cov[pos];
	      mapGcContent[mapIdx] = gcContent[pos];
	      ++mapIdx;
	    }
	  }
	}

	// Find breakpoints
	typedef boost::dynamic_bitset<> TBitSet;
	TBitSet bp(mappable, false);

	// Iterate a couple CNV sizes to find breakpoint regions
	uint32_t boundarySize = 10000;
	for(uint32_t cnvSize = c.minCnvSize; cnvSize <= 1000000; cnvSize = (uint32_t) (cnvSize * 1.3)) {
	  double covsumLeft = 0;
	  double expcovLeft = 0;
	  double covsumMiddle = 0;
	  double expcovMiddle = 0;
	  double covsumRight = 0;
	  double expcovRight = 0;
	  for(int32_t pos = boundarySize; pos < (int32_t) mappable - boundarySize - cnvSize; ++pos) {
	    if (pos == boundarySize) {
	      for(int32_t i = 0; i<boundarySize; ++i) {
		covsumLeft += mapcov[i];
		expcovLeft += gcbias[mapGcContent[i]].coverage;
	      }
	      for(int32_t i = boundarySize; i < boundarySize + cnvSize; ++i) {
		covsumMiddle += mapcov[i];
		expcovMiddle += gcbias[mapGcContent[i]].coverage;
	      }
	      for(int32_t i = boundarySize + cnvSize; i < boundarySize + cnvSize + boundarySize; ++i) {
		covsumRight += mapcov[i];
		expcovRight += gcbias[mapGcContent[i]].coverage;
	      }
	    } else {
	      covsumLeft -= mapcov[pos - boundarySize - 1];
	      expcovLeft -= gcbias[mapGcContent[pos - boundarySize - 1]].coverage;
	      covsumMiddle -= mapcov[pos - 1];
	      expcovMiddle -= gcbias[mapGcContent[pos - 1]].coverage;
	      covsumRight -= mapcov[pos + cnvSize - 1];
	      expcovRight -= gcbias[mapGcContent[pos + cnvSize - 1]].coverage;
	      covsumLeft += mapcov[pos];
	      expcovLeft += gcbias[mapGcContent[pos]].coverage;
	      covsumMiddle += mapcov[pos + cnvSize];
	      expcovMiddle += gcbias[mapGcContent[pos + cnvSize]].coverage;
	      covsumRight += mapcov[pos + cnvSize + boundarySize];
	      expcovRight += gcbias[mapGcContent[pos + cnvSize + boundarySize]].coverage;
	    }
	    double cnLeft = 2 * covsumLeft / expcovLeft;
	    double cnMiddle = 2 * covsumMiddle / expcovMiddle;
	    double cnRight = 2 * covsumRight / expcovRight;
	    if ((std::abs(cnLeft - cnRight) < 0.2) && (std::abs(cnLeft - cnMiddle) > 0.8) && (std::abs(cnRight - cnMiddle) > 0.8) && ((cnMiddle < 1.5) || (cnMiddle > 2.5))) {
	      bp[pos + cnvSize / 2] = true;
	    }
	  }
	}

	// Get start and end of CNVs
	int32_t svStart = -1;
	int32_t svEnd = -1;
	for(uint32_t i = 0; i < mappable; ++i) {
	  if (bp[i]) {
	    if (svStart != -1) svEnd = i;
	    else {
	      svStart = i;
	      svEnd = i;
	    }
	  } else {
	    if (svStart != -1) {
	      if ((svStart < svEnd) && (svEnd - svStart > 50)) {
		// Start from core region
		int32_t offset = (int32_t) (0.2 * (svEnd - svStart));
		svStart += offset;
		svEnd -= offset;
		if (svStart < svEnd) {
		  double covsum = 0;
		  double expcov = 0;
		  for(int32_t pos = svStart; pos <= svEnd; ++pos) {
		    covsum += mapcov[pos];
		    expcov += gcbias[mapGcContent[pos]].coverage;
		  }
		  double cn = 2 * covsum / expcov;

		  // Use x-drop extension
		  double runningCn = cn;
		  for(int32_t extension = 0; std::abs(runningCn - cn) < 0.2; ++extension) {
		    double runningCnLeft = -1;
		    double runningCnRight = -1;
		    if (svStart > 0) {
		      double covsumLeft = covsum + mapcov[svStart - 1];
		      double expcovLeft = expcov + gcbias[mapGcContent[svStart - 1]].coverage;
		      runningCnLeft =  2 * covsumLeft / expcovLeft;
		    }
		    if (svEnd + 1 < mappable) {
		      double covsumRight = covsum + mapcov[svEnd + 1];
		      double expcovRight = expcov + gcbias[mapGcContent[svEnd + 1]].coverage;
		      runningCnRight = 2 * covsumRight / expcovRight;
		    }
		    if ((runningCnLeft != -1) && (runningCnRight != -1)) {
		      if (std::abs(runningCnLeft - cn) < (runningCnRight - cn)) {
			--svStart;
			covsum += mapcov[svStart];
			expcov += gcbias[mapGcContent[svStart]].coverage;
		      } else {
			++svEnd;
			covsum += mapcov[svEnd];
			expcov += gcbias[mapGcContent[svEnd]].coverage;
		      }
		    } else if (runningCnLeft != -1) {
		      --svStart;
		      covsum += mapcov[svStart];
		      expcov += gcbias[mapGcContent[svStart]].coverage;
		    } else if (runningCnRight != -1) {
		      ++svEnd;
		      covsum += mapcov[svEnd];
		      expcov += gcbias[mapGcContent[svEnd]].coverage;
		    } else break;
		    runningCn = 2 * covsum / expcov;
		  }

		  // Remap position
		  if ((svEnd - svStart) >= c.minCnvSize) {
		    uint32_t mapIdx = 0;
		    int32_t svStartHg = 0;
		    int32_t svEndHg = 0;
		    for(uint32_t pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
		      if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
			if (mapIdx == svStart) svStartHg = pos;
			if (mapIdx == svEnd) svEndHg = pos;
			++mapIdx;
		      }
		    }
		    std::cerr << svStartHg << '\t' << svEndHg << '\t' << cn << std::endl;
		  }
		}
	      }
	      svStart = -1;
	      svEnd = -1;
	    }
	  }
	}
  */

}

#endif
