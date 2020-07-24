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

#ifndef BAF_H
#define BAF_H

#include <iostream>
#include <vector>
#include <fstream>

#include <boost/math/distributions/binomial.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "util.h"
#include "variants.h"

namespace coralns {


  inline double
  binomTest(uint32_t x, uint32_t n, double p) {
    boost::math::binomial binomialdist(n, p);
    double cutoff = pdf(binomialdist, x);
    double pval = 0.0;
    for(uint32_t k = 0; k <= n; ++k) {
      double p = pdf(binomialdist, k);
      if (p <= cutoff) pval +=p;
    }
    return pval;
  }
  
  template<typename TConfig, typename TGenomicVariants>
  inline double
  baf(TConfig const& c, bool const runControl, TGenomicVariants& cvar, TGenomicVariants& gvar) {
    double maf = 0;
    uint32_t mafCount = 0;
    
    // Open BAM file
    samFile* samfile;
    if (runControl) samfile = sam_open(c.controlFile.string().c_str(), "r");
    else samfile = sam_open(c.bamFile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx;
    if (runControl) idx = sam_index_load(samfile, c.controlFile.string().c_str());
    else idx = sam_index_load(samfile, c.bamFile.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Load bcf file
    htsFile* ibcffile = bcf_open(c.vcffile.string().c_str(), "r");
    hts_idx_t* bcfidx = bcf_index_load(c.vcffile.string().c_str());
    bcf_hdr_t* bcfhdr = bcf_hdr_read(ibcffile);

    // Parse BAM
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Variant Calling" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );
    
    // Iterate all chromosomes
    faidx_t* faiRef = fai_load(c.genome.string().c_str());
    for (int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      ++show_progress;
      if (chrNoData(c, refIndex, idx)) continue;

      // Load het. markers
      typedef std::vector<BiallelicVariant> TVariants;
      TVariants pv;
      if (!_loadVariants(ibcffile, bcfidx, bcfhdr, c.sampleName, hdr->target_name[refIndex], pv)) continue;
      if (pv.empty()) continue;

      // Sort variants
      std::sort(pv.begin(), pv.end(), SortVariants<BiallelicVariant>());

      // Check presence in reference
      std::string tname(hdr->target_name[refIndex]);
      int32_t seqlen = faidx_seq_len(faiRef, tname.c_str());
      if (seqlen == - 1) continue;
      else seqlen = -1;
      char* ref = faidx_fetch_seq(faiRef, tname.c_str(), 0, faidx_seq_len(faiRef, tname.c_str()), &seqlen);

      // Count REF and ALT support
      typedef std::vector<uint16_t> TAlleleSupport;
      int32_t maxCoverage = std::numeric_limits<uint16_t>::max();
      TAlleleSupport refS(pv.size(), 0);
      TAlleleSupport altS(pv.size(), 0);
      {
	hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	while (sam_itr_next(samfile, iter, rec) >= 0) {
	  if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	  if ((rec->core.flag & BAM_FPAIRED) && ((rec->core.flag & BAM_FMUNMAP) || (rec->core.tid != rec->core.mtid))) continue;
	  if (rec->core.qual < c.minQual) continue;
	  
	  // Fetch contained variants
	  typename TVariants::const_iterator vIt = std::lower_bound(pv.begin(), pv.end(), BiallelicVariant(rec->core.pos), SortVariants<BiallelicVariant>());
	  typename TVariants::const_iterator vItEnd = std::upper_bound(pv.begin(), pv.end(), BiallelicVariant(lastAlignedPosition(rec)), SortVariants<BiallelicVariant>());
	  if (vIt != vItEnd) {
	    // Get read sequence
	    std::string sequence;
	    sequence.resize(rec->core.l_qseq);
	    uint8_t* seqptr = bam_get_seq(rec);
	    for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

	    // Get base qualities
	    typedef std::vector<uint8_t> TQuality;
	    TQuality quality;
	    quality.resize(rec->core.l_qseq);
	    uint8_t* qualptr = bam_get_qual(rec);
	    for (int i = 0; i < rec->core.l_qseq; ++i) quality[i] = qualptr[i];
	    
	    // Parse CIGAR
	    uint32_t* cigar = bam_get_cigar(rec);
	    for(;vIt != vItEnd; ++vIt) {
	      uint32_t gp = rec->core.pos; // Genomic position
	      uint32_t sp = 0; // Sequence position
	      for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
		if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) sp += bam_cigar_oplen(cigar[i]);
		else if (bam_cigar_op(cigar[i]) == BAM_CINS) sp += bam_cigar_oplen(cigar[i]);
		else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) sp += bam_cigar_oplen(cigar[i]);
		else if (bam_cigar_op(cigar[i]) == BAM_CDEL) gp += bam_cigar_oplen(cigar[i]);
		else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) gp += bam_cigar_oplen(cigar[i]);
		else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
		  //Nop
		} else if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
		  if (gp + (int32_t) bam_cigar_oplen(cigar[i]) < vIt->pos) {
		    gp += bam_cigar_oplen(cigar[i]);
		    sp += bam_cigar_oplen(cigar[i]);
		  } else {
		    for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]); ++k, ++sp, ++gp) {
		      if (gp == vIt->pos) {
			if (quality[sp] >= c.minBaseQual) {
			  // Check REF and ALT alleles
			  std::string rnuc = boost::to_upper_copy(std::string(1, ref[gp]));
			  if ((gp < hdr->target_len[refIndex]) && (vIt->ref == rnuc[0]) && (sp < sequence.size())) {
			    if ((sequence[sp] == vIt->alt) && (altS[vIt-pv.begin()] < maxCoverage)) ++altS[vIt-pv.begin()];
			    else if ((sequence[sp] == vIt->ref) && (refS[vIt-pv.begin()] < maxCoverage)) ++refS[vIt-pv.begin()];
			  }
			}
		      }
		    }
		  }
		}
		else std::cerr << "Warning: Unknown Cigar options!" << std::endl;
	      }
	    }
	  }
	}
	bam_destroy1(rec);
	hts_itr_destroy(iter);
      }
      // Clean-up      
      if (ref != NULL) free(ref);

      // Filter for het. variants
      for (uint32_t i = 0; i<pv.size(); ++i) {
	if (cvar.empty()) {
	  // Control or tumor-only
	  // Require min. coverage
	  uint32_t totalcov = refS[i] + altS[i];
	  if (totalcov > c.minCoverage) {
	    double minHet = 0.1;
	    if (c.hasControlFile) minHet = 0.3; // More stringent because control available
	    uint32_t minSupport = std::max((int32_t) std::ceil(minHet * (double) totalcov), 2);
	    if ((altS[i] >= minSupport) && (refS[i] >= minSupport)) {
	      double vaf = (double) altS[i] / (double) (altS[i] + refS[i]);
	      if (vaf < 0.5) maf += vaf;
	      else maf += (1 - vaf);
	      gvar[refIndex].push_back(BiallelicSupport(pv[i].pos, refS[i], altS[i]));
	      ++mafCount;
	    }
	  }
	} else {
	  // Matched control available
	  // Make sure variant is present in control het. variant set
	  typename TGenomicVariants::value_type::const_iterator vIt = std::lower_bound(cvar[refIndex].begin(), cvar[refIndex].end(), BiallelicSupport(pv[i].pos), SortVariants<BiallelicSupport>());
	  if ((vIt != cvar[refIndex].end()) && (vIt->pos == pv[i].pos)) {
	    double vaf = (double) altS[i] / (double) (altS[i] + refS[i]);
	    if (vaf < 0.5) maf += vaf;
	    else maf += (1 - vaf);
	    gvar[refIndex].push_back(BiallelicSupport(pv[i].pos, refS[i], altS[i]));
	    ++mafCount;
	  }
	}
      }

      // Sort by position
      std::sort(gvar[refIndex].begin(), gvar[refIndex].end(), SortVariants<BiallelicSupport>());
    }

    // Close BCF
    bcf_hdr_destroy(bcfhdr);
    hts_idx_destroy(bcfidx);
    bcf_close(ibcffile);
    
    // Clean-up
    fai_destroy(faiRef);
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);

    return (double) maf /  (double) mafCount;
  }


  template<typename TConfig, typename TGenomicVariants>
  inline double
  baf(TConfig const& c, bool const runControl, TGenomicVariants& gvar) {
    TGenomicVariants cvar;
    return baf(c, runControl, cvar, gvar);
  }
  
}

#endif
