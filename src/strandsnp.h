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

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/tokenizer.hpp>


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


  inline std::string
  _parseBamHeader(bam_hdr_t const* hdr) {
    typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
    boost::char_separator<char> sep("\n");
    std::string header(hdr->text, 0, hdr->l_text);
    Tokenizer tokens(header, sep);
    for(Tokenizer::iterator tokIter = tokens.begin(); tokIter!=tokens.end(); ++tokIter) {
      if ((tokIter->at(0) == '@') && (tokIter->at(1) == 'R') && (tokIter->at(2) == 'G')) {
	boost::char_separator<char> sp("\t");
	Tokenizer fields(*tokIter, sp);
	for(Tokenizer::iterator itF = fields.begin(); itF!=fields.end(); ++itF) {
	  if ((itF->at(0) == 'S') && (itF->at(1) == 'M') && (itF->at(2) == ':')) {
	    std::string smtag(*itF);
	    return smtag.substr(3);
	  }
	}
      }
    }
    return "UnknownSample";
  }


  template<typename TConfig, typename TGenomicSnps, typename TFileCounts>
  inline void
  outputVCF(TConfig const& c, TGenomicSnps const& snps, TFileCounts const& fCount) {
    typedef typename TFileCounts::value_type TGenomicCounts;
    
    // Open wc bam file (not sorted and no index yet!!!)
    samFile* samfile = sam_open(c.wc.string().c_str(), "r");
    bam_hdr_t* bamhd = sam_hdr_read(samfile);
    std::string smtag(_parseBamHeader(bamhd));

    // Info
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Writing phased VCF" << std::endl;
    boost::progress_display show_progress( 2*bamhd->n_targets );

    // Output all structural variants
    htsFile *fp = hts_open(c.outvcf.string().c_str(), "wg");
    bcf_hdr_t *hdr = bcf_hdr_init("w");
    
    // Print vcf header
    boost::gregorian::date today = now.date();
    std::string datestr("##fileDate=");
    datestr += boost::gregorian::to_iso_string(today);
    bcf_hdr_append(hdr, datestr.c_str());
    bcf_hdr_append(hdr, "##FILTER=<ID=LowQual,Description=\"Read support below 3 or mapping quality below 20.\">");
    bcf_hdr_append(hdr, "##INFO=<ID=H1DP,Number=1,Type=Integer,Description=\"Number of high-quality reads in H1\">");
    bcf_hdr_append(hdr, "##INFO=<ID=H2DP,Number=1,Type=Integer,Description=\"Number of high-quality reads in H2\">");
    bcf_hdr_append(hdr, "##INFO=<ID=H1DP4,Number=4,Type=Integer,Description=\"Number of high-quality ref-fwd, ref-reverse, alt-fwd and alt-reverse bases for H1\">");
    bcf_hdr_append(hdr, "##INFO=<ID=H2DP4,Number=4,Type=Integer,Description=\"Number of high-quality ref-fwd, ref-reverse, alt-fwd and alt-reverse bases for H2\">");
    bcf_hdr_append(hdr, "##INFO=<ID=PSMETHOD,Number=1,Type=String,Description=\"Phasing method\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    // Add reference
    for (int i = 0; i<bamhd->n_targets; ++i) {
      std::string refname("##contig=<ID=");
      refname += std::string(bamhd->target_name[i]) + ",length=" + boost::lexical_cast<std::string>(bamhd->target_len[i]) + ">";
      bcf_hdr_append(hdr, refname.c_str());
    }
    // Add sample
    uint32_t numSample = 1;
    std::string sampleName(smtag);
    bcf_hdr_add_sample(hdr, sampleName.c_str());
    bcf_hdr_add_sample(hdr, NULL);
    bcf_hdr_write(fp, hdr);

    // Count by haplotype
    TGenomicCounts snpsH1(bamhd->n_targets);
    TGenomicCounts snpsH2(bamhd->n_targets);
    for (int refIndex = 0; refIndex<bamhd->n_targets; ++refIndex) {
      if (bamhd->target_len[refIndex] < c.window) continue;
      snpsH1[refIndex].resize(snps[refIndex].size(), RACount());
      snpsH2[refIndex].resize(snps[refIndex].size(), RACount());
    }
    bam1_t* r = bam_init1();
    while (sam_read1(samfile, bamhd, r) >= 0) {
      if (r->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
      if ((r->core.qual < c.minMapQual) || (r->core.tid<0)) continue;

      uint8_t *psptr = bam_aux_get(r, "PS");
      if (psptr) {
	// Only consider reads belonging to a phased block
	char* ps = (char*) (psptr + 1);
	std::string rPS = std::string(ps);
	
	// Get the haplotype
	uint8_t* hpptr = bam_aux_get(r, "HP");
	if (hpptr) {
	  int hap = bam_aux2i(hpptr);
	  if (hap == 1) _refAltCount(r, snps[r->core.tid], snpsH1[r->core.tid]);
	  else if (hap == 2) _refAltCount(r, snps[r->core.tid], snpsH2[r->core.tid]);
	}
      }
    }
    bam_destroy1(r);
    show_progress += bamhd->n_targets;

    for (int refIndex = 0; refIndex<bamhd->n_targets; ++refIndex) {
      ++show_progress;
      if (bamhd->target_len[refIndex] < c.window) continue;

      // Find informative het. SNPs
      typedef std::set<std::size_t> THetSnp;
      THetSnp hetSNP;
      _extractHetSNP(c, fCount, refIndex, hetSNP);

      // Iterate all SNPs
      int32_t *gts = (int*) malloc(numSample * 2 * sizeof(int));
      bcf1_t *rec = bcf_init();
      for(typename THetSnp::const_iterator itHS = hetSNP.begin(); itHS != hetSNP.end(); ++itHS) {
	bool h1Allele = false;
	bool h1Success = _getWatsonAllele(snpsH1[refIndex][*itHS], h1Allele);
	bool h2Allele = false;
	bool h2Success = _getCrickAllele(snpsH2[refIndex][*itHS], h2Allele);
	if ((h1Success) && (h2Success) && (h1Allele != h2Allele)) {
	  // Output main vcf fields
	  rec->rid = bcf_hdr_name2id(hdr, bamhd->target_name[refIndex]);
	  rec->pos = snps[refIndex][*itHS].pos;
	  rec->qual = 0;
	  std::string id(".");
	  bcf_update_id(hdr, rec, id.c_str());
	  std::string alleles;
	  alleles.resize(3);
	  alleles[0] = snps[refIndex][*itHS].ref;
	  alleles[1] = ',';
	  alleles[2] = snps[refIndex][*itHS].alt;
	  bcf_update_alleles_str(hdr, rec, alleles.c_str());
	  int32_t tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");
	  bcf_update_filter(hdr, rec, &tmpi, 1);
	  
	  // Add INFO fields
	  std::string psmethod("StrandSeq_v0.0.1");
	  bcf_update_info_string(hdr,rec, "PSMETHOD", psmethod.c_str());
	  int32_t h1dp4[4];
	  h1dp4[0] = snpsH1[refIndex][*itHS].watsonRef;
	  h1dp4[1] = snpsH1[refIndex][*itHS].crickRef;
	  h1dp4[2] = snpsH1[refIndex][*itHS].watsonAlt;
	  h1dp4[3] = snpsH1[refIndex][*itHS].crickAlt;
	  int32_t h1dp = h1dp4[0] + h1dp4[1] + h1dp4[2] + h1dp4[3];
	  bcf_update_info_int32(hdr, rec, "H1DP", &h1dp, 1);
	  bcf_update_info_int32(hdr, rec, "H1DP4", h1dp4, 4);
	  int32_t h2dp4[4];
	  h2dp4[0] = snpsH2[refIndex][*itHS].watsonRef;
	  h2dp4[1] = snpsH2[refIndex][*itHS].crickRef;
	  h2dp4[2] = snpsH2[refIndex][*itHS].watsonAlt;
	  h2dp4[3] = snpsH2[refIndex][*itHS].crickAlt;
	  int32_t h2dp = h2dp4[0] + h2dp4[1] + h2dp4[2] + h2dp4[3];
	  bcf_update_info_int32(hdr, rec, "H2DP", &h2dp, 1);
	  bcf_update_info_int32(hdr, rec, "H2DP4", h2dp4, 4);
	  
	  // Add genotypes
	  //bcf_gt_unphased(0);
	  if (h1Allele) gts[0] = bcf_gt_phased(1);
	  else gts[0] = bcf_gt_phased(0);
	  if (h2Allele) gts[1] = bcf_gt_phased(1);
	  else gts[1] = bcf_gt_phased(0);
	  bcf_update_genotypes(hdr, rec, gts, numSample * 2);
	  bcf_write1(fp, hdr, rec);
	  bcf_clear1(rec);
	}
      }

      // Clean-up
      free(gts);
      bcf_destroy1(rec);
    }
    
    // Close BAM file
    bam_hdr_destroy(bamhd);
    sam_close(samfile);
    
    // Close VCF file
    bcf_hdr_destroy(hdr);
    hts_close(fp);
  }



}

#endif
