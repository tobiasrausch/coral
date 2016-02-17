/*
============================================================================
Strand-Seq Watershed Segmentation
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

#ifndef STRANDWATERSHED_H
#define STRANDWATERSHED_H

#include <assert.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>

namespace streq {

  struct WaterComp {
    int32_t leftInd;
    int32_t rightInd;
    int32_t minInd;

    WaterComp(int32_t ind) : leftInd(ind), rightInd(ind), minInd(ind) {}
  };

  template<typename TValue>
  struct Peak {
    int32_t minInd;
    int32_t maxInd;
    TValue shift;

    Peak<TValue>() : minInd(0), maxInd(0), shift(0) {}
    Peak<TValue>(int32_t m, int32_t a, TValue s) : minInd(m), maxInd(a), shift(s) {}

    bool operator<(const Peak<TValue>& other) const {
      return ((shift < other.shift) || ((shift == other.shift) && (minInd < other.minInd)));
    }
  };


  // Returns (global minimum, index) and all paired extrema
  template<typename TValue, typename TPeaks>
  inline std::pair<TValue, int32_t>
  watershed(std::vector<TValue> const& data, TPeaks& peaks) {
    typedef std::vector<TValue> TValueVector;
    typedef typename TPeaks::value_type TPeak;

    typedef std::pair<TValue, int32_t> TDataIndex;
    typedef std::vector<TDataIndex> TDataIndexVec;
    TDataIndexVec sdata; 
    std::vector<int32_t> colors(data.size(), -1);
    std::vector<WaterComp> comp;
    uint32_t totalComp = 0;
    std::size_t ind = 0;

    for(typename TValueVector::const_iterator itV = data.begin(); itV != data.end(); ++itV, ++ind) sdata.push_back(std::make_pair(*itV, ind));
    std::sort(sdata.begin(), sdata.end());
    
    for (typename TDataIndexVec::iterator itDI = sdata.begin(); itDI != sdata.end(); ++itDI) {
      int ind = itDI->second;
      if (ind == 0) {
	if (colors[ind + 1] == -1) {
	  comp.push_back(WaterComp(ind)); 
	  colors[ind] = totalComp++;
	} else {
	  colors[ind] = colors[ind+1];
	  comp[colors[ind]].leftInd = ind;
	}
      } else if (ind == data.size() - 1) {
	if (colors[ind - 1] == -1) {
	  comp.push_back(WaterComp(ind)); 
	  colors[ind] = totalComp++;
	} else {
	  colors[ind] = colors[ind-1];
	  comp[colors[ind]].rightInd = ind;
	}
      } else {
	if (colors[ind-1] == -1 && colors[ind+1] == -1) {
	  // local min
	  comp.push_back(WaterComp(ind)); 
	  colors[ind] = totalComp++;
	} else if (colors[ind-1] != -1 && colors[ind+1] == -1) {
	  colors[ind] = colors[ind-1];
	  comp[colors[ind]].rightInd = ind;
	} else if (colors[ind-1] == -1 && colors[ind+1] != -1) {
	  colors[ind] = colors[ind+1];
	  comp[colors[ind]].leftInd = ind;
	} else if (colors[ind-1] != -1 && colors[ind+1] != -1) {
	  // local max
	  int firstInd = (data[comp[colors[ind+1]].minInd] < data[comp[colors[ind-1]].minInd]) ? comp[colors[ind-1]].minInd : comp[colors[ind+1]].minInd;
	  int secondInd = ind;
	  if (data[secondInd] > data[firstInd]) peaks.push_back(TPeak(firstInd, secondInd, data[secondInd] - data[firstInd]));
	  else if (data[firstInd] > data[secondInd]) peaks.push_back(TPeak(secondInd, firstInd, data[firstInd] - data[secondInd]));
	  else if (firstInd < secondInd) peaks.push_back(TPeak(firstInd, secondInd, data[secondInd] - data[firstInd]));
	  else peaks.push_back(TPeak(secondInd, firstInd, data[firstInd] - data[secondInd]));
	  
	  // merge components
	  int32_t survivorComp = colors[ind-1];
	  int32_t destroyedComp = colors[ind+1];
	  if ((data[comp[survivorComp].minInd] > data[comp[destroyedComp].minInd]) || ((data[comp[survivorComp].minInd] == data[comp[destroyedComp].minInd]) && (survivorComp > destroyedComp))) {
	    survivorComp = colors[ind+1];
	    destroyedComp = colors[ind-1];
	  }
	  // change color of the edges
	  colors[comp[destroyedComp].rightInd] = survivorComp;
	  colors[comp[destroyedComp].leftInd] = survivorComp;
    
	  // update edge indexes
	  if (comp[survivorComp].minInd > comp[destroyedComp].minInd) comp[survivorComp].leftInd = comp[destroyedComp].leftInd;
	  else comp[survivorComp].rightInd = comp[destroyedComp].rightInd;

	  // color ind
	  colors[ind] = colors[ind-1];
	}
      }
    }
    std::sort(peaks.begin(), peaks.end());
    return std::make_pair(data[comp[0].minInd], comp[0].minInd);
  }

	

}
#endif
