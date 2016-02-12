/*
============================================================================
Strand-Seq Utility Functions
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

#ifndef STRANDUTIL_H
#define STRANDUTIL_H

namespace streq
{

  inline uint32_t alignmentLength(bam1_t const* rec) {
    uint32_t* cigar = bam_get_cigar(rec);
    uint32_t alen = 0;
    for (std::size_t i = 0; i < rec->core.n_cigar; ++i)
      if (bam_cigar_op(cigar[i]) == BAM_CMATCH) alen+=bam_cigar_oplen(cigar[i]);
    return alen;
  }

  inline uint32_t halfAlignmentLength(bam1_t const* rec) {
    return (alignmentLength(rec) / 2);
  }

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


  template<typename TCounter>
  inline uint32_t
  _getReadSupportPercentile(TCounter const& watsonCount, TCounter const& crickCount, uint32_t perc) {
    TCounter support(watsonCount.size(), 0);
    typename TCounter::iterator itSupport = support.begin();
    typename TCounter::const_iterator itCrick =crickCount.begin();
    for(typename TCounter::const_iterator itWatson = watsonCount.begin(); itWatson != watsonCount.end(); ++itWatson, ++itCrick, ++itSupport)
      *itSupport = *itWatson + *itCrick;
    std::sort(support.begin(), support.end());
    return support[(int32_t) ((perc * watsonCount.size()) / 100)];
  }

  template<typename TWRatioVector>
  inline std::pair<float, float>
  _biModalMinima(TWRatioVector& wRatio) {
    std::sort(wRatio.begin(), wRatio.end());
    uint32_t binCount = (std::pow(wRatio.size(), 0.1/0.3) * (wRatio[wRatio.size()-1] - wRatio[0])) / (wRatio[(int) (3*wRatio.size()/4)] - wRatio[(int) (wRatio.size()/4)]);
    std::vector<uint32_t> histogram(binCount + 1, 0);
    for(typename TWRatioVector::iterator itW = wRatio.begin(); itW != wRatio.end(); ++itW) ++histogram[(uint32_t)((*itW)*binCount)];
    uint32_t smallestBinVal = histogram[0];
    float crickCut = 0;
    for(std::size_t i = 0; i<binCount/2; ++i) {
      if (histogram[i] < smallestBinVal) {
	smallestBinVal = histogram[i];
	crickCut = (float) i / (float) binCount;
      }
    }
    smallestBinVal = histogram[binCount];
    float watsonCut = 1;
    for(std::size_t i = binCount; i>binCount/2; --i) {
      if (histogram[i] < smallestBinVal) {
	smallestBinVal = histogram[i];
	watsonCut = (float) i / (float) binCount;
      }
    }
    return std::make_pair(crickCut, watsonCut);
  }

  template<typename TFractionVector>
  inline std::pair<float, float>
  _getWWStrandBounds(TFractionVector& watfraction) {
    float watmedian = 0;
    _getMedian(watfraction.begin(), watfraction.end(), watmedian);
    float watmad = 0;
    _getMAD(watfraction.begin(), watfraction.end(), watmedian, watmad);
    return std::make_pair(watmedian, watmad);
  }

}

#endif
