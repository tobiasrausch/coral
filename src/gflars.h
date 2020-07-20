#ifndef GFLARS_H
#define GFLARS_H

namespace coralns
{

  struct Recap {
    typedef double TPrecision;
    typedef std::vector<TPrecision> TPrecVector;
    typedef std::vector<uint32_t> TIndexVector;
    typedef boost::multi_array<TPrecision, 2> TSignalMatrix;
    typedef boost::unordered_map<uint32_t, TSignalMatrix> TSignalMap;
    
    uint32_t kbest;
    TIndexVector kbestjump;
    TPrecVector lambda;
    TPrecVector meansignal;
    TIndexVector jump;
    TSignalMap value;
    TSignalMatrix smooth;
    TSignalMatrix updown;
  };
  

  template<typename TSignalMatrix, typename TValue>
  inline void
  defaultweights(TSignalMatrix const& sm, std::vector<TValue>& weights) {
    TValue n = (TValue) sm.shape()[0];
    weights.resize(n);
    for(uint32_t i = 1; i <= n; ++i) weights[i-1] = (TValue) std::sqrt((TValue) n / (TValue) (i * (n - i)));
  }

  template<typename TConfig, typename TSignalMatrix>
  inline void
  gflars(TConfig const& c, TSignalMatrix const& sm, Recap& res) {
    typedef Recap::TPrecision TPrecision;
    typedef Recap::TPrecVector TPrecVector;
    
    res.lambda.resize(c.k,0);
    res.jump.resize(c.k,0);
    
    // Get weights
    TPrecVector w;
    defaultweights(sm, w);
    colmean(sm, res.meansignal);
    
    TSignalMatrix cHat;
    leftmultiplybyXt(sm, w, cHat);
    
    for(uint32_t iter = 0; iter < c.k; ++iter) {
      std::size_t nrow = cHat.shape()[0];
      std::size_t ncol = cHat.shape()[1];
      TPrecVector cHatSquareNorm(nrow, 0);
      TPrecision bigcHat = 0;
      uint32_t besti = 0;
      for(uint32_t i = 0; i < nrow; ++i) {
	for(uint32_t j = 0; j < ncol; ++j) cHatSquareNorm[i] += cHat[i][j] * cHat[i][j];
	if (cHatSquareNorm[i] > bigcHat) {
	  bigcHat = cHatSquareNorm[i];
	  besti = i;
	}
	//std::cout << cHatSquareNorm[i] << std::endl;
      }
      //std::cout << bigcHat << ',' << besti << std::endl;
      
      if (!iter) res.jump[0] = besti;
      typedef std::pair<uint32_t, uint32_t> TValIdxPair;
      std::vector<TValIdxPair> sub;
      for(uint32_t i = 0; i<= iter; ++i) sub.push_back(std::make_pair(res.jump[i], i));
      std::sort(sub.begin(), sub.end());
      std::vector<uint32_t> Amatrix(iter+1);
      std::vector<uint32_t> Imatrix(iter+1);
      for(uint32_t i = 0; i<= iter; ++i) {
	Amatrix[i] = sub[i].first;
	Imatrix[i] = sub[i].second;
      }
      
      TSignalMatrix smallw; 
      leftmultiplybyinvXAtXA(Amatrix, cHat, w, smallw);
      // Debug
      //for(uint32_t i = 0; i < smallw.shape()[0]; ++i) {
      //for(uint32_t j = 0; j < smallw.shape()[1]; ++j) std::cout << smallw[i][j] << ',';
      //std::cout << std::endl;  }
      
      TSignalMatrix smalla;
      multiplyXtXbysparse(Amatrix, smallw, w, smalla);
      
      TPrecVector a1(smalla.shape()[0], 0);
      for(uint32_t i = 0; i < smalla.shape()[0]; ++i) {
	a1[i] = bigcHat;
	for(uint32_t j = 0; j < smalla.shape()[1]; ++j) a1[i] -= smalla[i][j] * smalla[i][j];
	//std::cout << a1[i] << std::endl;
      }
      TPrecVector a2(cHat.shape()[0], 0);
      for(uint32_t i = 0; i < cHat.shape()[0]; ++i) {
	a2[i] = bigcHat;
	for(uint32_t j = 0; j < cHat.shape()[1]; ++j) a2[i] -= smalla[i][j] * cHat[i][j];
	//std::cout << a2[i] << std::endl;
      }
      TPrecVector a3(nrow, 0);
      for(uint32_t i = 0; i < a3.size(); ++i) {
	a3[i] = bigcHat - cHatSquareNorm[i];
	//std::cout << a3[i] << std::endl;
      }

      TPrecVector gammaTemp(2*nrow, 0);
      for(uint32_t i = 0; i < a1.size(); ++i) {
	if (a1[i] > c.epsilon) {
	  gammaTemp[i] = (a2[i] + std::sqrt(a2[i]*a2[i] - a1[i]*a3[i]))/a1[i];
	  gammaTemp[i+nrow-1] = (a2[i] - std::sqrt(a2[i]*a2[i] - a1[i]*a3[i]))/a1[i];
	}
      }
      for(uint32_t i = 0; i < a1.size(); ++i) {
	if ((a1[i] <= c.epsilon) && (a2[i] > c.epsilon)) {
	  gammaTemp[i] = a3[i] / (2*a2[i]);
	  gammaTemp[i+nrow-1] = a3[i] / (2*a2[i]);
	}
      }
      TPrecision maxg = gammaTemp[0] + 1;
      for(uint32_t i = 0; i < gammaTemp.size(); ++i)
	if (gammaTemp[i] + 1 > maxg) maxg = gammaTemp[i] + 1;
      for(uint32_t i = 0; i < a1.size(); ++i) {
	if ((a1[i] <= c.epsilon) && (a2[i] <= c.epsilon)) {
	  gammaTemp[i] = maxg;
	  gammaTemp[i + nrow] = maxg; 
	}
      }
      for(uint32_t i = 0; i < Amatrix.size(); ++i) {
	gammaTemp[Amatrix[i]] = maxg;
	gammaTemp[Amatrix[i] + nrow - 1] = maxg;
      }
      for(uint32_t i = 0; i < gammaTemp.size(); ++i) {
	if (gammaTemp[i] <= 0) gammaTemp[i] = maxg;
	else {
	  std::complex<TPrecision> z(gammaTemp[i]);
	  if (z.imag() != 0) gammaTemp[i] = maxg;
	}
      }
      TPrecision gamma = gammaTemp[0];
      uint32_t nexttoadd = 0;
      for(uint32_t i = 0; i < gammaTemp.size(); ++i) {
	if (gammaTemp[i] < gamma) {
	  gamma = gammaTemp[i];
	  nexttoadd = i;
	}
      }
      //std::cout << gamma << ',' << nexttoadd << ',' << std::endl;
      
      TSignalMatrix resvalue(boost::extents[iter+1][ncol]);
      for(uint32_t i = 0; i < Imatrix.size(); ++i)
	for(uint32_t j = 0; j < ncol; ++j) resvalue[i][j] = gamma * smallw[Imatrix[i]][j];
      if (iter > 0) {
	for(uint32_t i = 0; i < iter; ++i)
	  for(uint32_t j = 0; j < ncol; ++j)
	    resvalue[i][j] += res.value[iter-1][i][j];
      }
      res.value.insert(std::make_pair(iter, resvalue));
      res.lambda[iter] = std::sqrt(bigcHat);
      if (iter + 1 < c.k) {
	res.jump[iter+1] = 1 + ((nexttoadd-1) % (nrow - 1));
	//res.jump[iter+1] = ( nexttoadd % (nrow - 1));
	for(uint32_t i = 0; i < nrow; ++i)
	  for(uint32_t j = 0; j < ncol; ++j)
	    cHat[i][j] -= gamma * smalla[i][j];
      }
    }
  }

}

#endif


