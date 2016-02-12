/*
============================================================================
Strand-Seq SNP Functions
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

#ifndef STRANDSNP_H
#define STRANDSNP_H

namespace streq
{

  struct Snp {
    uint32_t pos;
    char ref;
    char alt;
    
    Snp(uint32_t p, char r, char a) : pos(p), ref(r), alt(a) {}
  };

  struct RACount {
    uint16_t watsonRef;
    uint16_t watsonAlt;
    uint16_t crickRef;
    uint16_t crickAlt;
    uint16_t falseBase;
  
    RACount() : watsonRef(0), watsonAlt(0), crickRef(0), crickAlt(0), falseBase(0) {}
    RACount(uint16_t wr, uint16_t wa, uint16_t cr, uint16_t ca, uint16_t f) : watsonRef(wr), watsonAlt(wa), crickRef(cr), crickAlt(ca), falseBase(f) {}
  };


  template<typename TRecord>
  struct SortSnps : public std::binary_function<TRecord, TRecord, bool> {
      inline bool operator()(TRecord const& s1, TRecord const& s2) const {
	return s1.pos < s2.pos;
      }
    };

  template<typename TRefAltCount>
  inline bool
  _getWatsonAllele(TRefAltCount const& ra, bool& allele) {
    if (ra.watsonRef + ra.watsonAlt > 0) {
      if (ra.watsonRef > 2*ra.watsonAlt) {
	allele = false;
	return true;
      } else if (2 * ra.watsonRef < ra.watsonAlt) {
	allele = true;
	return true;
      }
    }
    return false;
  }

  template<typename TRefAltCount>
  inline bool
  _getCrickAllele(TRefAltCount const& ra, bool& allele) {
    if (ra.crickRef + ra.crickAlt > 0) {
      if (ra.crickRef > 2*ra.crickAlt) {
	allele = false;
	return true;
      } else if (2 * ra.crickRef < ra.crickAlt) {
	allele = true;
	return true;
      }
    }
    return false;
  }


  template<typename TSnpVector, typename TCountVector>
  inline void
  _refAltCount(bam1_t const* rec, TSnpVector const& snps, TCountVector& chrCounts) {
    typedef typename TSnpVector::value_type TSnp;
  
    // Annotate SNPs
    typename TSnpVector::const_iterator iSnp = std::lower_bound(snps.begin(), snps.end(), TSnp(rec->core.pos, 'A', 'A'), SortSnps<TSnp>());
    typename TSnpVector::const_iterator iSnpEnd = std::upper_bound(snps.begin(), snps.end(), TSnp(rec->core.pos + alignmentLength(rec), 'A', 'A'), SortSnps<TSnp>());
    if (iSnp != iSnpEnd) {
      std::string sequence;
      sequence.resize(rec->core.l_qseq);
      uint8_t* seqptr = bam_get_seq(rec);
      for (std::size_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
      uint32_t* cigar = bam_get_cigar(rec);
      for(;iSnp != iSnpEnd; ++iSnp) {
	int32_t gp = rec->core.pos; // Genomic position
	int32_t sp = 0; // Sequence position
	bool foundChar = false;
	for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	  if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) sp += bam_cigar_oplen(cigar[i]);
	  else if (bam_cigar_op(cigar[i]) == BAM_CINS) sp += bam_cigar_oplen(cigar[i]);
	  else if (bam_cigar_op(cigar[i]) == BAM_CDEL) gp += bam_cigar_oplen(cigar[i]);
	  else if (bam_cigar_op(cigar[i]) == BAM_CMATCH) {
	    if (gp + bam_cigar_oplen(cigar[i]) < iSnp->pos) {
	      gp += bam_cigar_oplen(cigar[i]);
	      sp += bam_cigar_oplen(cigar[i]);
	    } else {
	      for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]); ++k, ++sp, ++gp) {
		if (gp == iSnp->pos) {
		  foundChar = true;
		  break;
		}
	      }
	      if (foundChar) break;
	    }
	  }
	}
	if (foundChar) {
	  if (sequence[sp] == iSnp->ref) { 
	    if (rec->core.flag & BAM_FREAD1) 
	      if (rec->core.flag & BAM_FREVERSE) ++chrCounts[(uint32_t) (iSnp - snps.begin())].crickRef;
	      else ++chrCounts[(uint32_t) (iSnp - snps.begin())].watsonRef;
	    else
	      if (rec->core.flag & BAM_FREVERSE) ++chrCounts[(uint32_t) (iSnp - snps.begin())].watsonRef;
	      else ++chrCounts[(uint32_t) (iSnp - snps.begin())].crickRef;
	  } else if (sequence[sp] == iSnp->alt) {
	    if (rec->core.flag & BAM_FREAD1) 
	      if (rec->core.flag & BAM_FREVERSE) ++chrCounts[(uint32_t) (iSnp - snps.begin())].crickAlt;
	      else ++chrCounts[(uint32_t) (iSnp - snps.begin())].watsonAlt;
	    else
	      if (rec->core.flag & BAM_FREVERSE) ++chrCounts[(uint32_t) (iSnp - snps.begin())].watsonAlt;
	      else ++chrCounts[(uint32_t) (iSnp - snps.begin())].crickAlt;
	  } else ++chrCounts[(uint32_t) (iSnp - snps.begin())].falseBase;
	}
      }
    }
  }

  template<typename TConfig, typename TGenomicSnps>
  inline int32_t 
  _loadVariationData(TConfig const& c, bam_hdr_t const* hdr, TGenomicSnps& snps) {
    typedef typename TGenomicSnps::value_type TSnpVector;
    typedef typename TSnpVector::value_type TSnp;

    // Open VCF file
    htsFile* ifile = bcf_open(c.variation.string().c_str(), "r");
    if (ifile == NULL) {
      std::cerr << "SNP VCF files is missing " << c.variation.string() << std::endl;
      return 1;
    }
    hts_idx_t* bcfidx = bcf_index_load(c.variation.string().c_str());
    if (bcfidx == NULL) {
      std::cerr << "SNP VCF index file is missing " << c.variation.string() << std::endl;
      return 1;
    }
    bcf_hdr_t* vcfh = bcf_hdr_read(ifile);
    
    // Load SNPs
    snps.clear();
    snps.resize(hdr->n_targets);
    for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
      if (hdr->target_len[refIndex] < c.window) continue;
      std::string chrName(hdr->target_name[refIndex]);
      uint32_t chrid = bcf_hdr_name2id(vcfh, chrName.c_str());
      hts_itr_t* itervcf = bcf_itr_queryi(bcfidx, chrid, 0, hdr->target_len[refIndex]);
      bcf1_t* var = bcf_init();
      while (bcf_itr_next(ifile, itervcf, var) >= 0) {
	bcf_unpack(var, BCF_UN_STR);
	std::vector<std::string> alleles;
	for(std::size_t i = 0; i<var->n_allele; ++i) alleles.push_back(std::string(var->d.allele[i]));
	// Only bi-allelic SNPs
	if ((alleles.size() == 2) && (alleles[0].size() == 1) && (alleles[1].size() == 1)) snps[refIndex].push_back(TSnp(var->pos, alleles[0][0], alleles[1][0]));
      }
      bcf_destroy(var);
      hts_itr_destroy(itervcf);
      
      // Sort Snps by position
      std::sort(snps[refIndex].begin(), snps[refIndex].end(), SortSnps<TSnp>());
    }
    // Close VCF
    bcf_hdr_destroy(vcfh);
    hts_idx_destroy(bcfidx);
    bcf_close(ifile);
    
    return 0;
  }



}

#endif
