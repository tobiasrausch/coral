#ifndef MATRIX_H
#define MATRIX_H

namespace coralns
{

  template<typename TSignalMatrix, typename TValue>
  inline void
  colmean(TSignalMatrix const& sm, std::vector<TValue>& cm) {
    std::size_t nrow = sm.shape()[0];
    std::size_t ncol = sm.shape()[1];
    cm.resize(ncol, 0);
    for(uint32_t i = 0; i < nrow; ++i)
      for(uint32_t j = 0; j < ncol; ++j)
	cm[j] += sm[i][j];
    for(uint32_t j = 0; j < ncol; ++j) cm[j] /= nrow;
  }
  
  template<typename TSignalMatrix>
  inline void
  colcumsum(TSignalMatrix const& sm, TSignalMatrix& cs) {
    std::size_t nrow = sm.shape()[0];
    std::size_t ncol = sm.shape()[1];
    cs.resize(boost::extents[nrow][ncol]);
    if (nrow > 0) {
      for(uint32_t j = 0; j < ncol; ++j) cs[0][j] = sm[0][j];
      for(uint32_t i = 1; i < nrow; ++i)
	for(uint32_t j = 0; j < ncol; ++j)
	  cs[i][j] = cs[i-1][j] + sm[i][j];
    }
  }
  
  template<typename TSignalMatrix, typename TValue>
  inline void
  leftmultiplybyXt(TSignalMatrix const& sm, std::vector<TValue> const& weights, TSignalMatrix& cHat) {
    std::size_t nrow = sm.shape()[0];
    std::size_t ncol = sm.shape()[1];
    
    TSignalMatrix colcs;
    colcumsum(sm, colcs);
    cHat.resize(boost::extents[nrow - 1][ncol]);
    for(uint32_t i = 0; i < nrow - 1; ++i) 
      for(uint32_t j = 0; j < ncol; ++j) cHat[i][j] = ((i+1) * colcs[nrow-1][j] / nrow - colcs[i][j]) * weights[i];
  }
  
  template<typename TSignalMatrix>
  inline void
  flipud(TSignalMatrix const& s, TSignalMatrix& sud) {
    uint32_t nrow = s.shape()[0];
    uint32_t ncol = s.shape()[1];
    sud.resize(boost::extents[nrow][ncol]);
    for(uint32_t i = 0; i < nrow; ++i)
      for(uint32_t j = 0; j < ncol; ++j)
	sud[nrow-i-1][j] = s[i][j];
  }
  
  template<typename TSignalMatrix, typename TValue>
  inline void
  multiplyXtXbysparse(std::vector<uint32_t> const& Amatrix, TSignalMatrix const& smallw, std::vector<TValue> const& weights, TSignalMatrix& r) {
    uint32_t nrow = smallw.shape()[0];
    uint32_t ncol = smallw.shape()[1];
    TSignalMatrix tmpVal(boost::extents[nrow][ncol]);
    for(uint32_t i = 0; i < Amatrix.size(); ++i)
      for(uint32_t j = 0; j < ncol; ++j)
	tmpVal[i][j] = smallw[i][j] * weights[Amatrix[i]];
    
    uint32_t n = weights.size();
    TSignalMatrix s(boost::extents[n-1][ncol]);
    for(uint32_t i = 0; i < (n-1); ++i)
      for(uint32_t j = 0; j < ncol; ++j)
	s[i][j] = 0;
    for(uint32_t i = 0; i < Amatrix.size(); ++i)
      for(uint32_t j = 0; j < ncol; ++j)
	s[Amatrix[i]][j] = tmpVal[i][j];
    
    
    TSignalMatrix sud;
    flipud(s, sud);
    TSignalMatrix cssud;
    colcumsum(sud, cssud);
    flipud(cssud, s);

    std::vector<TValue> u(ncol, 0);
    for(uint32_t j = 0; j < ncol; ++j)
      for(uint32_t i = 0; i < Amatrix.size(); ++i)
	u[j] += tmpVal[i][j] * (TValue) (Amatrix[i] + 1);
    
    for(uint32_t i = 0; i < (n-1); ++i)
      for(uint32_t j = 0; j < ncol; ++j)
	s[i][j] = s[i][j] - u[j]/(TValue)n;
    colcumsum(s, r);
    
    for(uint32_t i = 0; i < (n-1); ++i)
      for(uint32_t j = 0; j < ncol; ++j)
	r[i][j] *= weights[i];
  }
  
  
  template<typename TSignalMatrix, typename TValue>
  inline void
  leftmultiplybyinvXAtXA(std::vector<uint32_t> const& Amatrix, TSignalMatrix const& cHat, std::vector<TValue> const& weights, TSignalMatrix& r) {
    std::vector<uint32_t> diff(Amatrix.size() + 1);
    int32_t last = -1;
    for(uint32_t i = 0; i < Amatrix.size(); ++i) {
      diff[i] = (uint32_t) (Amatrix[i] - last);
      last = Amatrix[i];
    }
    diff[Amatrix.size()] = (uint32_t) ((weights.size() - 1) - last);
    
    uint32_t ncol = cHat.shape()[1];
    TSignalMatrix cHatSub(boost::extents[Amatrix.size()][ncol]);
    TSignalMatrix delta(boost::extents[Amatrix.size() + 1][ncol]);
    for(uint32_t i = 0; i < Amatrix.size(); ++i) {
      for(uint32_t j = 0; j < ncol; ++j) {
	cHatSub[i][j] = cHat[Amatrix[i]][j] / weights[Amatrix[i]];
	if (!i) delta[i][j] = cHatSub[i][j] / diff[i];
	else delta[i][j] = (cHatSub[i][j] - cHatSub[i-1][j]) / diff[i];
      }
    }
    for(uint32_t j = 0; j < ncol; ++j) delta[Amatrix.size()][j] = (0 - cHatSub[Amatrix.size()-1][j]) / diff[Amatrix.size()];
    
    r.resize(boost::extents[Amatrix.size()][ncol]);
    for(uint32_t i = 1; i < delta.shape()[0]; ++i) {
      for(uint32_t j = 0; j < delta.shape()[1]; ++j) {
	r[i-1][j] = (-1 * (delta[i][j] - delta[i-1][j])) / weights[Amatrix[i-1]];
      }
    }
  }
  
}

#endif


