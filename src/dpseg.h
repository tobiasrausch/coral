#ifndef DPSEG_H
#define DPSEG_H

namespace coralns
{

    template<typename TConfig, typename TSignalMatrix>
    inline void
    dpseg(TConfig const& c, TSignalMatrix const& sm, Recap& res) {
      typedef Recap::TPrecision TPrecision;
      typedef Recap::TPrecVector TPrecVector;
      
      uint32_t nrow = sm.shape()[0];
      uint32_t ncol = sm.shape()[1];      
      uint32_t kmax = (uint32_t) std::min((double) res.jump.size(), std::floor((double) nrow / 4.0));

      std::vector<int32_t> b(res.jump.begin(), res.jump.end());
      b.push_back(-1);
      b.push_back(nrow-1);
      std::sort(b.begin(), b.end());
      uint32_t k = b.size() - 1;

      TSignalMatrix csSM;
      colcumsum(sm, csSM);
      TSignalMatrix smalls(boost::extents[nrow+1][ncol]);
      TPrecVector smallv(nrow+1, 0);
      for(uint32_t j = 0; j < ncol; ++j) smalls[0][j] = 0;
      for(uint32_t i = 1; i < nrow + 1; ++i) {
	for(uint32_t j = 0; j < ncol; ++j) {
	  smalls[i][j] = csSM[i-1][j];
	  smallv[i] += sm[i-1][j]*sm[i-1][j];
	}
	smallv[i] += smallv[i-1];
      }

      TSignalMatrix jmat(boost::extents[k][k]);
      for(uint32_t i = 0; i < k; ++i) {
	for(uint32_t j = 0; j < k; ++j) {
	  int32_t istart = b[i] + 1;
	  int32_t iend = b[j+1];
	  TPrecision cs = 0;
	  for(uint32_t col = 0; col < ncol; ++col) cs += (smalls[iend + 1][col] - smalls[istart][col]) * (smalls[iend + 1][col] - smalls[istart][col]);
	  cs /= (TPrecision) (iend - istart + 1);
	  jmat[i][j] = smallv[iend + 1] - smallv[istart] - cs;
	}
      }

      TSignalMatrix vmat(boost::extents[kmax+1][k]);
      typedef boost::multi_array<uint32_t, 2> TIndexMatrix;
      TIndexMatrix jump(boost::extents[kmax][k]);
      for(uint32_t j = 0; j < k; ++j) vmat[0][j] = jmat[0][j];
      for(uint32_t ki = 0; ki < kmax; ++ki) {
	for(uint32_t j = ki + 1; j < k; ++j) {
	  TPrecision val = vmat[ki][ki] + jmat[ki + 1][j];
	  uint32_t ind = ki;
	  for(uint32_t itrange = ki; itrange<=j-1; ++itrange) {
	    if (vmat[ki][itrange] + jmat[itrange + 1][j] < val) {
	      val = vmat[ki][itrange] + jmat[itrange + 1][j];
	      ind = itrange;
	    }
	  }
	  vmat[ki+1][j] = val;
	  jump[ki][j] = ind;
	}
      }

      typedef std::vector<uint32_t> TIndexVector;
      typedef std::vector<TIndexVector> TIndexVectorMap;
      TIndexVectorMap resjump(kmax, TIndexVector());      
      for(uint32_t ki = 0; ki < kmax; ++ki) {
	resjump[ki].resize(ki+1, 0);
	resjump[ki][ki] = jump[ki][k-1];
	for(int32_t i = (int32_t)(ki) - 1; i >= 0; --i) resjump[ki][i] = jump[i][resjump[ki][i+1]];
      }

      for(uint32_t ki = 0; ki < kmax; ++ki)
	for(uint32_t i = 0; i < resjump[ki].size(); ++i) resjump[ki][i] = b[resjump[ki][i]+1];

      std::vector<TPrecision> jvec(kmax+1, 0);
      for(uint32_t ki = 0; ki < kmax + 1; ++ki) jvec[ki] = std::log(vmat[ki][k-1]);
      std::vector<TPrecision> jtild(kmax+1, 0);
      for(uint32_t ki = 0; ki < kmax + 1; ++ki) jtild[ki] = (jvec[kmax] - jvec[ki])/(jvec[kmax] - jvec[0]) * kmax + 1;
      std::vector<TPrecision> diffjtild(kmax, 0);
      for(uint32_t ki = 0; ki < kmax; ++ki) diffjtild[ki] = jtild[ki+1] - jtild[ki];
      int32_t bestIdx = -1;
      for(uint32_t ki = 0; ki < kmax - 1; ++ki)
	if (diffjtild[ki+1] - diffjtild[ki] > c.dpthreshold) bestIdx = ki + 1;
      if (bestIdx == -1) bestIdx = 0;
      res.kbest = bestIdx;
      res.kbestjump = resjump[res.kbest];
    }

}

#endif


