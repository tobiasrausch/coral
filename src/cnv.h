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

  struct CNV {
    int32_t chr;
    int32_t start;
    int32_t end;
    int32_t ciposlow;
    int32_t ciposhigh;
    int32_t ciendlow;
    int32_t ciendhigh;
    float cn;
    float rdsupport;
    float penalty;
    float mappable;
    
  CNV(int32_t const c, int32_t const s, int32_t const e, int32_t const cil, int32_t const cih, int32_t const cel, int32_t ceh, float const estcn, float const sp, float const pty, float const mp) : chr(c), start(s), end(e), ciposlow(cil), ciposhigh(cih), ciendlow(cel), ciendhigh(ceh), cn(estcn), rdsupport(sp), penalty(pty), mappable(mp) {}
  };

    struct PeakInfo {
      uint32_t pos;
      uint32_t cilow;
      uint32_t cihigh;
      float peakHeight;
      float penalty;

    PeakInfo(uint32_t const p, uint32_t const cl, uint32_t const ch, float const ph) : pos(p), cilow(cl), cihigh(ch), peakHeight(ph), penalty(0) {}
    };


  
  template<typename TCoverage, typename TSplits>
  inline void
  _collectSplitBp(TCoverage& splitCovLeft, TCoverage& splitCovRight, TSplits& splitBp, uint32_t const bandwidth) {
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
	uint32_t smax = std::max(splitCovLeft[i], splitCovRight[i]);
	splitBp.push_back(std::make_pair(i, smax));
      }
    }

    // Sort breakpoints
    std::sort(splitBp.begin(), splitBp.end());
  }


  template<typename TConfig, typename TGcBias, typename TCoverage>
  inline void
    callCNVs(TConfig const& c, std::pair<uint32_t, uint32_t> const& gcbound, std::vector<uint16_t> const& gcContent, std::vector<uint16_t> const& uniqContent, TGcBias const& gcbias, TCoverage const& cov, bam_hdr_t const* hdr, int32_t const refIndex, std::vector<CNV>& cnvCalls) {

    uint32_t rdThreshold = 5;
    
    // Shrink to mappable pos
    int32_t mappable = 0;
    for(uint32_t pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
      if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) ++mappable;
    }
    TCoverage mapcov(mappable, 0);
    std::vector<uint16_t> mapGcContent(mappable, 0);
    uint32_t mapIdx = 0;
    for(uint32_t pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
      if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	mapcov[mapIdx] = cov[pos];
	mapGcContent[mapIdx] = gcContent[pos];
	++mapIdx;
      }
    }


    // Debug
    //std::map<uint32_t, uint32_t> debugMap;
    //{
    //uint32_t mapIdx = 0;
    //for(uint32_t pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
    //	if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	  //if ((pos > 23557000) && (pos < 23562000)) debugMap.insert(std::make_pair(mapIdx, pos));
    //	  if ((pos > 47290575) && (pos < 47310982)) debugMap.insert(std::make_pair(mapIdx, pos));	
    //	  ++mapIdx;
    //	}
    //}
    //}
    
    // Iterate a couple CNV sizes to find breakpoint regions
    std::vector<float> totaldiff(mappable, 0);
    for(int32_t cnvSize = c.minCnvSize / 2; cnvSize <= 100000; cnvSize = (int32_t) (cnvSize * 1.3)) {
      double covsumLeft = 0;
      double expcovLeft = 0;
      double covsumRight = 0;
      double expcovRight = 0;
      int32_t spacer = 100;
      for(int32_t pos = 0; pos < mappable - cnvSize - spacer; ++pos) {
	if (pos < cnvSize) {
	  covsumLeft += mapcov[pos];
	  expcovLeft += gcbias[mapGcContent[pos]].coverage;
	  covsumRight += mapcov[pos + spacer + cnvSize];
	  expcovRight += gcbias[mapGcContent[pos + spacer + cnvSize]].coverage;
	} else {
	  covsumLeft -= mapcov[pos - cnvSize];
	  expcovLeft -= gcbias[mapGcContent[pos - cnvSize]].coverage;
	  covsumRight -= mapcov[pos + spacer];
	  expcovRight -= gcbias[mapGcContent[pos + spacer]].coverage;
	  covsumLeft += mapcov[pos];
	  expcovLeft += gcbias[mapGcContent[pos]].coverage;
	  covsumRight += mapcov[pos + spacer + cnvSize];
	  expcovRight += gcbias[mapGcContent[pos + spacer + cnvSize]].coverage;
	  if ((expcovLeft > 0) && (expcovRight > 0)) {
	    double cnLeft = c.ploidy * covsumLeft / expcovLeft;
	    double cnRight = c.ploidy * covsumRight / expcovRight;
	    double cnShift = std::abs(cnLeft - cnRight);
	    double offset = std::abs(boost::math::lround(cnShift) - cnShift) + std::abs(boost::math::lround(cnLeft) - cnLeft) + std::abs(boost::math::lround(cnRight) - cnRight);
	    if ((boost::math::lround(cnShift) >= 1) && (offset < 1)) {
	      double scaling = boost::math::lround(cnShift) - offset;
	      totaldiff[pos + spacer / 2] += boost::math::lround(cnShift) * scaling;
	    }
	  }
	}
      }
    }
    // Erase small peaks
    for(int32_t i = 0; i < mappable; ++i) {
      if (totaldiff[i] < ((double) rdThreshold / 2.0)) totaldiff[i] = 0;
    }

    // Debug
    //for(int32_t pos = 0; pos < mappable; ++pos) {
    //if (debugMap.find(pos) != debugMap.end()) std::cerr << pos << ',' << debugMap[pos] << ',' << totaldiff[pos] << ',' << mapcov[pos] << std::endl;
    //}

    // Find peaks
    typedef std::pair<uint32_t, uint32_t> TPeakPair;
    std::vector<TPeakPair> peaks;
    std::vector<PeakInfo> peakCoords;
    {
      // Find peak boundaries
      std::vector<uint32_t> peakMetrics;
      bool inPeak = false;
      uint32_t startPos = 0;
      uint32_t peakMax = 0;
      uint32_t peakIdx = 0;
      for(uint32_t i = 1; i < totaldiff.size(); ++i) {
	if (!inPeak) {
	  if (totaldiff[i] <= 0) startPos = i;
	  else {
	    inPeak = true;
	    peakMax = totaldiff[i];
	    peakIdx = i;
	  }
	} else {
	  if (totaldiff[i] <= 0) {
	    inPeak = false;
	    peakMetrics.push_back(startPos);
	    peakMetrics.push_back(peakIdx);
	    peakMetrics.push_back(i);
	    startPos = i;
	  } else {
	    if (totaldiff[i] > peakMax) {
	      peakMax = totaldiff[i];
	      peakIdx = i;
	    }
	  }
	}
      }
      
      // Debug
      //for(uint32_t i = 0; i < peakMetrics.size(); ++i ) {
      //if (debugMap.find(peakMetrics[i]) != debugMap.end()) std::cerr << "Peak:" << i << ',' << debugMap[peakMetrics[i]] << ',' << totaldiff[peakMetrics[i]] << std::endl;
      //}
      
      // Store peak position and height
      for(uint32_t i = 0; i < peakMetrics.size(); i += 3) {
	peakCoords.push_back(PeakInfo(peakMetrics[i+1], peakMetrics[i], peakMetrics[i+2], totaldiff[peakMetrics[i+1]]));
      }

      // Estimate CN in-between peaks
      std::vector<bool> peakUsed(peakCoords.size(), false);
      for(uint32_t peakThreshold = 20; peakThreshold >= rdThreshold; --peakThreshold) {
	uint32_t lastPeak = 0;
	for(uint32_t i = 1; i < peakCoords.size(); ++i) {
	  if (peakCoords[i].peakHeight >= peakThreshold) {
	    if (((peakCoords[i].pos - peakCoords[lastPeak].pos) >= c.minCnvSize) && (!peakUsed[lastPeak]) && (!peakUsed[i])) {
	      double covsum = 0;
	      double expcov = 0;
	      for(uint32_t k = peakCoords[lastPeak].pos; k < peakCoords[i].pos; ++k) {
		covsum += mapcov[k];
		expcov += gcbias[mapGcContent[k]].coverage;
	      }
	      double cn = c.ploidy * covsum / expcov;
	      if (boost::math::lround(cn) != c.ploidy) {
		int32_t sw = 5000;
		int32_t st = std::max(0, (int32_t) (peakCoords[lastPeak].pos) - sw);
		double covsumLeft = 0;
		double expcovLeft = 0;
		for(uint32_t k = st; k < peakCoords[lastPeak].pos; ++k) {
		  covsumLeft += mapcov[k];
		  expcovLeft += gcbias[mapGcContent[k]].coverage;
		}
		double cnLeft = c.ploidy * covsumLeft / expcovLeft;
		if ((boost::math::lround(cnLeft) != boost::math::lround(cn)) && (std::abs(cnLeft - cn) > 0.5)) {
		  int32_t ed = std::min((int32_t) hdr->target_len[refIndex], (int32_t) (peakCoords[i].pos) + sw);
		  double covsumRight = 0;
		  double expcovRight = 0;
		  for(int32_t k = peakCoords[i].pos; k < ed; ++k) {
		    covsumRight += mapcov[k];
		    expcovRight += gcbias[mapGcContent[k]].coverage;
		  }
		  double cnRight = c.ploidy * covsumRight / expcovRight;
		  if ((boost::math::lround(cnRight) != boost::math::lround(cn)) && (std::abs(cnRight -cn) > 0.5)) {
		    peakUsed[lastPeak] = true;
		    peakUsed[i] = true;
		    peaks.push_back(std::make_pair(lastPeak, i));

		    // Shift peaks to optimum CN-shift, left boundary
		    // Iterate to the right
		    uint32_t mi = peakCoords[lastPeak].pos;
		    double cnShift = std::abs(cnLeft - cn);
		    double offset = std::abs(boost::math::lround(cnShift) - cnShift) + std::abs(boost::math::lround(cn) - cn) + std::abs(boost::math::lround(cnLeft) - cnLeft);
		    double bestOffset = offset;
		    int32_t bestIdx = mi;
		    while ((offset < bestOffset + 0.1) && (mi < peakCoords[lastPeak].cihigh)) {
		      covsum -= mapcov[mi];
		      expcov -= gcbias[mapGcContent[mi]].coverage;
		      covsumLeft += mapcov[mi];
		      expcovLeft += gcbias[mapGcContent[mi]].coverage;
		      ++mi;
		      cn = c.ploidy * covsum / expcov;
		      cnLeft = c.ploidy * covsumLeft / expcovLeft;
		      cnShift = std::abs(cnLeft - cn);
		      offset = std::abs(boost::math::lround(cnShift) - cnShift) + std::abs(boost::math::lround(cn) - cn) + std::abs(boost::math::lround(cnLeft) - cnLeft);
		      if (bestOffset > offset) {
			bestOffset = offset;
			bestIdx = mi;
		      }
		    }
		    peakCoords[lastPeak].pos = bestIdx;
		    peakCoords[lastPeak].penalty = bestOffset;
		    // Update cn estimate
		    covsum = 0;
		    expcov = 0;
		    for(uint32_t k = peakCoords[lastPeak].pos; k < peakCoords[i].pos; ++k) {
		      covsum += mapcov[k];
		      expcov += gcbias[mapGcContent[k]].coverage;
		    }
		    cn = c.ploidy * covsum / expcov;
		    // Iterate to the left
		    st = std::max(0, (int32_t) (peakCoords[lastPeak].pos) - sw);
		    mi = peakCoords[lastPeak].pos;
		    covsumLeft = 0;
		    expcovLeft = 0;
		    for(uint32_t k = st; k < mi; ++k) {
		      covsumLeft += mapcov[k];
		      expcovLeft += gcbias[mapGcContent[k]].coverage;
		    }
		    cnLeft = c.ploidy * covsumLeft / expcovLeft;
		    cnShift = std::abs(cnLeft - cn);
		    offset = std::abs(boost::math::lround(cnShift) - cnShift) + std::abs(boost::math::lround(cn) - cn) + std::abs(boost::math::lround(cnLeft) - cnLeft);
		    bestOffset = offset;
		    bestIdx = mi;
		    while ((offset < bestOffset + 0.1) && (mi > peakCoords[lastPeak].cilow)) {
		      --mi;
		      covsum += mapcov[mi];
		      expcov += gcbias[mapGcContent[mi]].coverage;
		      covsumLeft -= mapcov[mi];
		      expcovLeft -= gcbias[mapGcContent[mi]].coverage;
		      cn = c.ploidy * covsum / expcov;
		      cnLeft = c.ploidy * covsumLeft / expcovLeft;
		      cnShift = std::abs(cnLeft - cn);
		      offset = std::abs(boost::math::lround(cnShift) - cnShift) + std::abs(boost::math::lround(cn) - cn) + std::abs(boost::math::lround(cnLeft) - cnLeft);
		      if (bestOffset > offset) {
			bestOffset = offset;
			bestIdx = mi;
		      }
		    }
		    peakCoords[lastPeak].pos = bestIdx;
		    peakCoords[lastPeak].penalty = bestOffset;
		    // Update cn estimate
		    covsum = 0;
		    expcov = 0;
		    for(uint32_t k = peakCoords[lastPeak].pos; k < peakCoords[i].pos; ++k) {
		      covsum += mapcov[k];
		      expcov += gcbias[mapGcContent[k]].coverage;
		    }
		    cn = c.ploidy * covsum / expcov;
		    // Adjust right peak
		    // Iterate to the left
		    mi = peakCoords[i].pos;
		    cnShift = std::abs(cnRight - cn);
		    offset = std::abs(boost::math::lround(cnShift) - cnShift) + std::abs(boost::math::lround(cn) - cn) + std::abs(boost::math::lround(cnRight) - cnRight);
		    bestOffset = offset;
		    bestIdx = mi;
		    while ((offset < bestOffset + 0.1) && (mi > peakCoords[i].cilow)) {
		      --mi;
		      covsum -= mapcov[mi];
		      expcov -= gcbias[mapGcContent[mi]].coverage;
		      covsumRight += mapcov[mi];
		      expcovRight += gcbias[mapGcContent[mi]].coverage;
		      cn = c.ploidy * covsum / expcov;
		      cnRight = c.ploidy * covsumRight / expcovRight;
		      cnShift = std::abs(cnRight - cn);
		      offset = std::abs(boost::math::lround(cnShift) - cnShift) + std::abs(boost::math::lround(cn) - cn) + std::abs(boost::math::lround(cnRight) - cnRight);
		      if (bestOffset > offset) {
			bestOffset = offset;
			bestIdx = mi;
		      }
		    }
		    peakCoords[i].pos = bestIdx;
		    peakCoords[i].penalty = bestOffset;
		    // Update cn estimate
		    covsum = 0;
		    expcov = 0;
		    for(uint32_t k = peakCoords[lastPeak].pos; k < peakCoords[i].pos; ++k) {
		      covsum += mapcov[k];
		      expcov += gcbias[mapGcContent[k]].coverage;
		    }
		    cn = c.ploidy * covsum / expcov;
		    // Adjust right peak
		    // Iterate to the right
		    ed = std::min((int32_t) hdr->target_len[refIndex], (int32_t) (peakCoords[i].pos) + sw);
		    mi = peakCoords[i].pos;
		    covsumRight = 0;
		    expcovRight = 0;
		    for(int32_t k = mi; k < ed; ++k) {
		      covsumRight += mapcov[k];
		      expcovRight += gcbias[mapGcContent[k]].coverage;
		    }
		    cnRight = c.ploidy * covsumRight / expcovRight;
		    cnShift = std::abs(cnRight - cn);
		    offset = std::abs(boost::math::lround(cnShift) - cnShift) + std::abs(boost::math::lround(cn) - cn) + std::abs(boost::math::lround(cnRight) - cnRight);
		    bestOffset = offset;
		    bestIdx = mi;
		    while ((offset < bestOffset + 0.1) && (mi < peakCoords[i].cihigh)) {
		      covsum += mapcov[mi];
		      expcov += gcbias[mapGcContent[mi]].coverage;
		      covsumRight -= mapcov[mi];
		      expcovRight -= gcbias[mapGcContent[mi]].coverage;
		      ++mi;
		      cn = c.ploidy * covsum / expcov;
		      cnRight = c.ploidy * covsumRight / expcovRight;
		      cnShift = std::abs(cnRight - cn);
		      offset = std::abs(boost::math::lround(cnShift) - cnShift) + std::abs(boost::math::lround(cn) - cn) + std::abs(boost::math::lround(cnRight) - cnRight);
		      //std::cerr << offset << ',' << bestOffset << ',' << mi << ',' << bestIdx << ',' << cnShift << ',' << cn << ',' << cnRight << ',' << mapcov[mi] << std::endl;
		      if (bestOffset > offset) {
			bestOffset = offset;
			bestIdx = mi;
		      }
		    }
		    peakCoords[i].pos = bestIdx;
		    peakCoords[i].penalty = bestOffset;
		    // Adjust confidence interval
		    // Update cn estimate
		    covsum = 0;
		    expcov = 0;
		    for(uint32_t k = peakCoords[lastPeak].pos; k < peakCoords[i].pos; ++k) {
		      covsum += mapcov[k];
		      expcov += gcbias[mapGcContent[k]].coverage;
		    }
		    cn = c.ploidy * covsum / expcov;
		    // Iterate to the left
		    st = std::max(0, (int32_t) (peakCoords[lastPeak].pos) - sw);
		    mi = peakCoords[lastPeak].pos;
		    covsumLeft = 0;
		    expcovLeft = 0;
		    for(uint32_t k = st; k < mi; ++k) {
		      covsumLeft += mapcov[k];
		      expcovLeft += gcbias[mapGcContent[k]].coverage;
		    }
		    cnLeft = c.ploidy * covsumLeft / expcovLeft;
		    cnShift = std::abs(cnLeft - cn);
		    offset = std::abs(boost::math::lround(cnShift) - cnShift) + std::abs(boost::math::lround(cn) - cn) + std::abs(boost::math::lround(cnLeft) - cnLeft);
		    bestOffset = offset;
		    while ((offset < bestOffset + 0.1) && (mi > peakCoords[lastPeak].cilow)) {
		      --mi;
		      covsum += mapcov[mi];
		      expcov += gcbias[mapGcContent[mi]].coverage;
		      covsumLeft -= mapcov[mi];
		      expcovLeft -= gcbias[mapGcContent[mi]].coverage;
		      cn = c.ploidy * covsum / expcov;
		      cnLeft = c.ploidy * covsumLeft / expcovLeft;
		      cnShift = std::abs(cnLeft - cn);
		      offset = std::abs(boost::math::lround(cnShift) - cnShift) + std::abs(boost::math::lround(cn) - cn) + std::abs(boost::math::lround(cnLeft) - cnLeft);
		    }
		    peakCoords[lastPeak].cilow = mi;
		    // Update cn estimate
		    covsum = 0;
		    expcov = 0;
		    for(uint32_t k = peakCoords[lastPeak].pos; k < peakCoords[i].pos; ++k) {
		      covsum += mapcov[k];
		      expcov += gcbias[mapGcContent[k]].coverage;
		    }
		    cn = c.ploidy * covsum / expcov;
		    // Iterate to the right
		    st = std::max(0, (int32_t) (peakCoords[lastPeak].pos) - sw);
		    mi = peakCoords[lastPeak].pos;
		    covsumLeft = 0;
		    expcovLeft = 0;
		    for(uint32_t k = st; k < mi; ++k) {
		      covsumLeft += mapcov[k];
		      expcovLeft += gcbias[mapGcContent[k]].coverage;
		    }
		    cnLeft = c.ploidy * covsumLeft / expcovLeft;
		    cnShift = std::abs(cnLeft - cn);
		    offset = std::abs(boost::math::lround(cnShift) - cnShift) + std::abs(boost::math::lround(cn) - cn) + std::abs(boost::math::lround(cnLeft) - cnLeft);
		    bestOffset = offset;
		    while ((offset < bestOffset + 0.1) && (mi < peakCoords[lastPeak].cihigh)) {
		      covsum -= mapcov[mi];
		      expcov -= gcbias[mapGcContent[mi]].coverage;
		      covsumLeft += mapcov[mi];
		      expcovLeft += gcbias[mapGcContent[mi]].coverage;
		      ++mi;
		      cn = c.ploidy * covsum / expcov;
		      cnLeft = c.ploidy * covsumLeft / expcovLeft;
		      cnShift = std::abs(cnLeft - cn);
		      offset = std::abs(boost::math::lround(cnShift) - cnShift) + std::abs(boost::math::lround(cn) - cn) + std::abs(boost::math::lround(cnLeft) - cnLeft);
		    }
		    peakCoords[lastPeak].cihigh = mi;
		    // Update cn estimate
		    covsum = 0;
		    expcov = 0;
		    for(uint32_t k = peakCoords[lastPeak].pos; k < peakCoords[i].pos; ++k) {
		      covsum += mapcov[k];
		      expcov += gcbias[mapGcContent[k]].coverage;
		    }
		    cn = c.ploidy * covsum / expcov;
		    // Adjust right confidence interval
		    // Iterate to the left
		    ed = std::min((int32_t) hdr->target_len[refIndex], (int32_t) (peakCoords[i].pos) + sw);
		    mi = peakCoords[i].pos;
		    covsumRight = 0;
		    expcovRight = 0;
		    for(int32_t k = mi; k < ed; ++k) {
		      covsumRight += mapcov[k];
		      expcovRight += gcbias[mapGcContent[k]].coverage;
		    }
		    cnRight = c.ploidy * covsumRight / expcovRight;
		    cnShift = std::abs(cnRight - cn);
		    offset = std::abs(boost::math::lround(cnShift) - cnShift) + std::abs(boost::math::lround(cn) - cn) + std::abs(boost::math::lround(cnRight) - cnRight);
		    bestOffset = offset;
		    while ((offset < bestOffset + 0.1) && (mi > peakCoords[i].cilow)) {
		      --mi;
		      covsum -= mapcov[mi];
		      expcov -= gcbias[mapGcContent[mi]].coverage;
		      covsumRight += mapcov[mi];
		      expcovRight += gcbias[mapGcContent[mi]].coverage;
		      cn = c.ploidy * covsum / expcov;
		      cnRight = c.ploidy * covsumRight / expcovRight;
		      cnShift = std::abs(cnRight - cn);
		      offset = std::abs(boost::math::lround(cnShift) - cnShift) + std::abs(boost::math::lround(cn) - cn) + std::abs(boost::math::lround(cnRight) - cnRight);
		    }
		    peakCoords[i].cilow = mi;
		    // Update cn estimate
		    covsum = 0;
		    expcov = 0;
		    for(uint32_t k = peakCoords[lastPeak].pos; k < peakCoords[i].pos; ++k) {
		      covsum += mapcov[k];
		      expcov += gcbias[mapGcContent[k]].coverage;
		    }
		    cn = c.ploidy * covsum / expcov;
		    // Iterate to the right
		    ed = std::min((int32_t) hdr->target_len[refIndex], (int32_t) (peakCoords[i].pos) + sw);
		    mi = peakCoords[i].pos;
		    covsumRight = 0;
		    expcovRight = 0;
		    for(int32_t k = mi; k < ed; ++k) {
		      covsumRight += mapcov[k];
		      expcovRight += gcbias[mapGcContent[k]].coverage;
		    }
		    cnRight = c.ploidy * covsumRight / expcovRight;
		    cnShift = std::abs(cnRight - cn);
		    offset = std::abs(boost::math::lround(cnShift) - cnShift) + std::abs(boost::math::lround(cn) - cn) + std::abs(boost::math::lround(cnRight) - cnRight);
		    bestOffset = offset;
		    while ((offset < bestOffset + 0.1) && (mi < peakCoords[i].cihigh)) {
		      covsum += mapcov[mi];
		      expcov += gcbias[mapGcContent[mi]].coverage;
		      covsumRight -= mapcov[mi];
		      expcovRight -= gcbias[mapGcContent[mi]].coverage;
		      ++mi;
		      cn = c.ploidy * covsum / expcov;
		      cnRight = c.ploidy * covsumRight / expcovRight;
		      cnShift = std::abs(cnRight - cn);
		      offset = std::abs(boost::math::lround(cnShift) - cnShift) + std::abs(boost::math::lround(cn) - cn) + std::abs(boost::math::lround(cnRight) - cnRight);
		    }
		    peakCoords[i].cihigh = mi;
		  }
		}
	      }
	    }
	    lastPeak = i;
	  }
	}
      }
    }

    // Translate back from mappable space to genomic space
    std::map<uint32_t, uint32_t> mapToGenomic;
    {      
      // Collect peak coordinates
      std::vector<uint32_t> coords;
      for(uint32_t i = 0; i < peaks.size(); ++i) {
	coords.push_back(peakCoords[peaks[i].first].cilow);
	coords.push_back(peakCoords[peaks[i].first].pos);
	coords.push_back(peakCoords[peaks[i].first].cihigh);
	coords.push_back(peakCoords[peaks[i].second].cilow);
	coords.push_back(peakCoords[peaks[i].second].pos);
	coords.push_back(peakCoords[peaks[i].second].cihigh);
      }
      std::sort(coords.begin(), coords.end());
      uint32_t pointer = 0;
      mapIdx = 0;
      for(uint32_t pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
	if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	  while ((pointer < coords.size()) && (coords[pointer] == mapIdx)) {
	    mapToGenomic.insert(std::make_pair(coords[pointer], pos));
	    ++pointer;
	  }
	  ++mapIdx;
	}
      }
    }

    // Translate peaks into genomic space
    for(uint32_t i = 0; i < peaks.size(); ++i) {
      double rdsupport = std::min(peakCoords[peaks[i].first].peakHeight, peakCoords[peaks[i].second].peakHeight);
      double penalty = std::max(peakCoords[peaks[i].first].penalty, peakCoords[peaks[i].second].penalty);
      double covsum = 0;
      double expcov = 0;
      for(uint32_t k = peakCoords[peaks[i].first].pos; k < peakCoords[peaks[i].second].pos; ++k) {
	covsum += mapcov[k];
	expcov += gcbias[mapGcContent[k]].coverage;
      }
      double cn = c.ploidy * covsum / expcov;
      int32_t svSize = mapToGenomic[peakCoords[peaks[i].second].pos] - mapToGenomic[peakCoords[peaks[i].first].pos];
      int32_t mapSize = peakCoords[peaks[i].second].pos - peakCoords[peaks[i].first].pos;
      double mappableFrac = (double) mapSize / (double) svSize;
      if (svSize >= (int32_t) c.minCnvSize) {
	cnvCalls.push_back(CNV(refIndex, mapToGenomic[peakCoords[peaks[i].first].pos], mapToGenomic[peakCoords[peaks[i].second].pos], mapToGenomic[peakCoords[peaks[i].first].cilow], mapToGenomic[peakCoords[peaks[i].first].cihigh], mapToGenomic[peakCoords[peaks[i].second].cilow], mapToGenomic[peakCoords[peaks[i].second].cihigh], cn, rdsupport, penalty, mappableFrac));
      }
    }
  }
  

}

#endif
