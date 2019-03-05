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

#ifndef VARIANTS_H
#define VARIANTS_H

#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <htslib/sam.h>


namespace coralns
{

  struct BiallelicVariant {
    uint32_t pos;
    char ref;
    char alt;

    explicit BiallelicVariant(uint32_t const p) : pos(p), ref('N'), alt('N') {}
    BiallelicVariant(uint32_t const p, char const r, char const a) : pos(p), ref(r), alt(a) {}
  };

  struct BiallelicSupport {
    uint32_t pos;
    uint16_t ref;
    uint16_t alt;
    explicit BiallelicSupport(uint32_t const p) : pos(p), ref(0), alt(0) {}
    BiallelicSupport(uint32_t const p, uint16_t const r, uint16_t const a) : pos(p), ref(r), alt(a) {}
  };
  
  template<typename TRecord>
  struct SortVariants : public std::binary_function<TRecord, TRecord, bool> {
    inline bool operator()(TRecord const& s1, TRecord const& s2) const {
      return s1.pos < s2.pos;
    }
  };

  template<typename TVarSupport>
  inline double
  mafSegment(uint32_t const s, uint32_t const e, int32_t const minSnps, TVarSupport const& vs) {
    uint32_t minMarkers = std::max(minSnps, 1);
    if ((!vs.empty()) && (s < e)) {
      typename TVarSupport::const_iterator vIt = std::lower_bound(vs.begin(), vs.end(), BiallelicSupport(s), SortVariants<BiallelicSupport>());
      typename TVarSupport::const_iterator vItEnd = std::upper_bound(vs.begin(), vs.end(), BiallelicSupport(e), SortVariants<BiallelicSupport>());
      if ((vIt != vs.end()) && (vIt != vItEnd)) {
	std::vector<double> mafvec;
	for(;vIt != vItEnd; ++vIt) {
	  if ((vIt->pos >= s) && (vIt->pos < e)) {
	    double maf = (double) vIt->alt / (double) (vIt->alt + vIt->ref);
	    if (maf > 0.5) maf = 1.0 - maf;
	    mafvec.push_back(maf);
	  }
	}
	if (mafvec.size() >= minMarkers) {
	  std::sort(mafvec.begin(), mafvec.end());
	  return mafvec[mafvec.size() / 2];
	}
      }
    }
    return -1;
  }

  template<typename TVariants>
  inline bool
  _loadVariants(htsFile* ifile, hts_idx_t* bcfidx, bcf_hdr_t* hdr, std::string const& sample, std::string const& chrom, TVariants& pV) {
    typedef typename TVariants::value_type TVariant;

    // Get sample index (if present)
    int32_t sampleIndex = -1;
    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i)
      if (hdr->samples[i] == sample) sampleIndex = i;
        
    // Genotypes
    int ngt = 0;
    int32_t* gt = NULL;

    // Collect het. bi-allelic variants for this chromosome
    int32_t chrid = bcf_hdr_name2id(hdr, chrom.c_str());
    int32_t lastpos = -1;
    if (chrid < 0) return false;
    hts_itr_t* itervcf = bcf_itr_querys(bcfidx, hdr, chrom.c_str());
    if (itervcf != NULL) {
      bcf1_t* rec = bcf_init1();
      while (bcf_itr_next(ifile, itervcf, rec) >= 0) {
	// Only bi-allelic variants
	if (rec->n_allele == 2) {
	  bcf_unpack(rec, BCF_UN_ALL);
	  bcf_get_genotypes(hdr, rec, &gt, &ngt);
	  bool includeVar = false;
	  if (sampleIndex != -1) {
	    if ((bcf_gt_allele(gt[sampleIndex*2]) != -1) && (bcf_gt_allele(gt[sampleIndex*2 + 1]) != -1) && (!bcf_gt_is_missing(gt[sampleIndex*2])) && (!bcf_gt_is_missing(gt[sampleIndex*2 + 1]))) {
	      int gt_type = bcf_gt_allele(gt[sampleIndex*2]) + bcf_gt_allele(gt[sampleIndex*2 + 1]);
	      // Only het. variants
	      if (gt_type == 1) includeVar = true;
	    }
	  } else includeVar = true;
	  if (includeVar) {
	    if (rec->pos != lastpos) {
	      // Only one variant per position
	      std::string ref = boost::to_upper_copy(std::string(rec->d.allele[0]));
	      std::string alt = boost::to_upper_copy(std::string(rec->d.allele[1]));
	      // Only SNPs
	      if ((ref.size()==1) && (alt.size()==1) && (ref[0] != alt[0])) {
		pV.push_back(TVariant(rec->pos, ref[0], alt[0]));
		lastpos = rec->pos;
	      }
	    }
	  }
	}
      }
      bcf_destroy(rec);
      hts_itr_destroy(itervcf);
    }
    if (gt != NULL) free(gt);
    return true;
  }

}

#endif
