/*
============================================================================
Strand-Seq Haplotyping
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

#ifndef STRANDHAPLO_H
#define STRANDHAPLO_H

namespace streq
{

  template<typename TConfig, typename TFileCounts, typename THetSnpSet>
  inline void
  _extractHetSNP(TConfig const& c, TFileCounts const& fCount, int32_t refIndex, THetSnpSet& hetSNP) {
    typedef typename TFileCounts::value_type TGenomicCounts;
    typedef typename TGenomicCounts::value_type TCountVector;

    // Sum-up ref and alt counts for the given chromosome
    TCountVector sumCounts(fCount[0][refIndex].size(), RACount());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      typename TCountVector::iterator itSC = sumCounts.begin();
      for(typename TCountVector::const_iterator itC = fCount[file_c][refIndex].begin(); itC != fCount[file_c][refIndex].end(); ++itC, ++itSC) {
	itSC->watsonRef += itC->watsonRef;
	itSC->watsonAlt += itC->watsonAlt;
	itSC->crickRef += itC->crickRef;
	itSC->crickAlt += itC->crickAlt;
	itSC->falseBase += itC->falseBase;
      }
    }
    
    for(typename TCountVector::const_iterator itSC = sumCounts.begin(); itSC != sumCounts.end(); ++itSC) {
      // No nucleotide mis-call seen
      if (itSC->falseBase == 0) {
	// REF and ALT observed
	if (((itSC->watsonRef + itSC->crickRef) > 0) && ((itSC->watsonAlt + itSC->crickAlt) > 0)) {
	  // It's an informative SNP about strands
	  if (((itSC->watsonRef + itSC->watsonAlt) > 0) && ((itSC->crickRef + itSC->crickAlt) > 0)) {
	    hetSNP.insert((std::size_t) (itSC - sumCounts.begin()));
	  }
	}
      }
    }
  }


  template<typename TConfig, typename TGenomicSnps, typename TFileCounts, typename TFileWRatio, typename TFileWindows>
  inline void
    _processVariationData(TConfig const& c, bam_hdr_t const* hdr, TGenomicSnps const& snps, TFileCounts const& fCount, TFileWRatio const& fWR, TFileWindows& wcWindows, TFileWindows& wcFlipWindows) {

    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Process variation data" << std::endl;
    boost::progress_display sprog( hdr->n_targets);
    for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
      ++sprog;
      if (hdr->target_len[refIndex] < c.window) continue;

      // Find informative het. SNPs
      typedef std::set<std::size_t> THetSnp;
      THetSnp hetSNP;
      _extractHetSNP(c, fCount, refIndex, hetSNP);

      // Build Haplotypes for wcWindows
      std::vector<uint32_t> support(c.files.size(), 0);
      typedef std::vector<bool> TFlipVector;
      typedef std::vector<TFlipVector> TFlipMatrix;
      TFlipMatrix fm;
      fm.resize(c.files.size(), TFlipVector());
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) fm[file_c].resize(c.files.size(), false);
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	if (wcWindows[file_c][refIndex].empty()) continue;
	for(unsigned int file_d = file_c + 1; file_d < c.files.size(); ++file_d) {
	  if (wcWindows[file_d][refIndex].empty()) continue;
	  uint32_t diffCount = 0;
	  uint32_t agreeCount = 0;
	  // Iterate informative SNPs
	  for(typename THetSnp::const_iterator itHS = hetSNP.begin(); itHS != hetSNP.end(); ++itHS) {
	    uint32_t curBin = snps[refIndex][*itHS].pos / c.window;
	    if ((fWR[file_c][refIndex][curBin] == -1) || (fWR[file_d][refIndex][curBin] == -1)) continue; // Blacklisted bin
	    if ((wcWindows[file_c][refIndex].find(curBin) != wcWindows[file_c][refIndex].end()) && (wcWindows[file_d][refIndex].find(curBin) != wcWindows[file_d][refIndex].end())) {
	      bool cWAllele = false;
	      bool watsonSuccess = _getWatsonAllele(fCount[file_c][refIndex][*itHS], cWAllele);
	      bool cCAllele = false;
	      bool crickSuccess = _getCrickAllele(fCount[file_c][refIndex][*itHS], cCAllele);
	      // At least one allele must be called
	      if ((watsonSuccess) && (!crickSuccess)) cCAllele = !cWAllele;
	      else if ((!watsonSuccess) && (crickSuccess)) cWAllele = !cCAllele;
	      else if ((!watsonSuccess) && (!crickSuccess)) continue;  // Not covered
	      else if (cWAllele == cCAllele) continue; // Incorrect genotyping
	      bool dWAllele = false;
	      watsonSuccess = _getWatsonAllele(fCount[file_d][refIndex][*itHS], dWAllele);
	      bool dCAllele = false;
	      crickSuccess = _getCrickAllele(fCount[file_d][refIndex][*itHS], dCAllele);
	      // At least one allele must be called
	      if ((watsonSuccess) && (!crickSuccess)) dCAllele = !dWAllele;
	      else if ((!watsonSuccess) && (crickSuccess)) dWAllele = !dCAllele;
	      else if ((!watsonSuccess) && (!crickSuccess)) continue;  // Not covered
	      else if (dWAllele == dCAllele) continue; // Incorrect genotyping

	      // Same haplotype?
	      if (cWAllele == dWAllele) ++agreeCount;
	      else ++diffCount;

	      // Debug code
	      //std::cerr << file_c << ',' << refIndex << ',' << snps[refIndex][*itHS].pos << ',' <<  snps[refIndex][*itHS].ref << ',' <<  snps[refIndex][*itHS].alt << ',' << fCount[file_c][refIndex][*itHS].watsonRef << ',' << fCount[file_c][refIndex][*itHS].watsonAlt << ',' << fCount[file_c][refIndex][*itHS].crickRef << ',' << fCount[file_c][refIndex][*itHS].crickAlt << ',' << std::endl;
	      //std::cerr << file_d << ',' << refIndex << ',' << snps[refIndex][*itHS].pos << ',' <<  snps[refIndex][*itHS].ref << ',' <<  snps[refIndex][*itHS].alt << ',' << fCount[file_d][refIndex][*itHS].watsonRef << ',' << fCount[file_d][refIndex][*itHS].watsonAlt << ',' << fCount[file_d][refIndex][*itHS].crickRef << ',' << fCount[file_d][refIndex][*itHS].crickAlt << ',' << std::endl;
	    }
	  }
	  if (diffCount>agreeCount) {
	    support[file_c] += diffCount;
	    support[file_d] += diffCount;
	    fm[file_c][file_d] = true;
	    fm[file_d][file_c] = true;
	  } else {
	    support[file_c] += agreeCount;
	    support[file_d] += agreeCount;
	  }
	}
      }
      // Find best covered sample (most informative het. SNPs)
      uint32_t bestSupport = 0;
      uint32_t bestIndex = 0;
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	if (support[file_c] > bestSupport) {
	  bestSupport = support[file_c];
	  bestIndex = file_c;
	}
      }
      if (bestSupport > 0) {
	for(unsigned int file_d = 0; file_d < c.files.size(); ++file_d) {
	  if (file_d == bestIndex) continue; // No flip required
	  if (fm[bestIndex][file_d]) {
	    // Flip required
	    if (wcWindows[file_d][refIndex].empty()) continue;
	    uint32_t bins = hdr->target_len[refIndex] / c.window + 1;
	    for(uint32_t bin = 0; bin < bins; ++bin) {
	      if (wcWindows[file_d][refIndex].find(bin) != wcWindows[file_d][refIndex].end()) {
		wcFlipWindows[file_d][refIndex].insert(bin);
		wcWindows[file_d][refIndex].erase(bin);
	      }
	    }
	  }
	}
      }
    }
  }


}

#endif
