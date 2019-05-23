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

#ifndef GMM_H
#define GMM_H

#include <boost/math/distributions/normal.hpp>
#include <htslib/sam.h>
#include <htslib/vcf.h>


namespace coralns
{

  #define SMALLEST_GL -1000
  #define MAX_CN 9
  
   // Convert string to char*
  struct cstyle_str {
    const char* operator ()(const std::string& s) {
      return s.c_str();
    }
  };

  struct SDAggregator {
    uint32_t wsinit;
    std::vector<double> cnSUM;
    std::vector<uint32_t> cnCount;

    SDAggregator(uint32_t wsi) : wsinit(wsi) {
      for(uint32_t ws = wsi; ws < 1000000; ws = ws * 2) {
	cnSUM.push_back(0);
        cnCount.push_back(0);
      }
    }
  };
  
  template<typename TValue>
  inline double
  getSD(SDAggregator const& sda, TValue const svsize) {
    double defaultSD = 0.5;
    uint32_t widx = 0;
    for(uint32_t ws = sda.wsinit; ws < 1000000; ws = ws * 2) {
      if ((TValue) (2 * ws) > svsize) {
	if (sda.cnCount[widx] > 1) {
	  double sd = sda.cnSUM[widx] / (double) (sda.cnCount[widx] - 1);
	  if (sd > defaultSD) return defaultSD;
	  else return sd;
	} else return defaultSD;
      }
      ++widx;
    }
    return defaultSD;
  }

  template<typename TPrecision>
    struct BoLog {
      typedef TPrecision value_type;
      
      std::vector<TPrecision> phred2prob;

      BoLog() {
	for(int i = 0; i <= boost::math::round(-10 * SMALLEST_GL); ++i) phred2prob.push_back(std::pow(TPrecision(10), -(TPrecision(i)/TPrecision(10))));
      }
    };


  template<typename TBoLog>
  inline void
  _computeGLs(TBoLog const& bl, SDAggregator const& sda, CNV const& cnv, float* gls, int32_t* gqval, int32_t* gts) {
    typedef typename TBoLog::value_type FLP;
    FLP gl[MAX_CN+1];
    
    // Compute CN probabilities
    for(uint32_t geno=0; geno<=MAX_CN; ++geno) gl[geno]=0;
    for(uint32_t geno=0; geno<=MAX_CN; ++geno) {
      double sd = getSD(sda, cnv.end - cnv.start);
      boost::math::normal gaus(geno, sd);
      gl[geno] = boost::math::pdf(gaus, cnv.cn);
    }

    // Normalize probabilities
    FLP totalGeno = 0;
    for(uint32_t geno=0; geno<=MAX_CN; ++geno) totalGeno += gl[geno];
    for(uint32_t geno=0; geno<=MAX_CN; ++geno) {
      gl[geno] = std::log10(gl[geno] / totalGeno);
    }

    // Re-scale to best GL
    unsigned int glBest=0;
    FLP glBestVal=gl[glBest];
    for(unsigned int geno=1; geno<=MAX_CN; ++geno) {
      if (gl[geno] > glBestVal) {
	glBestVal=gl[geno];
	glBest = geno;
      }
    }
    // Rescale by best genotype
    for(unsigned int geno=0; geno<=MAX_CN; ++geno) {
      gl[geno] -= glBestVal;
      // Cap at smallest GL
      gl[geno] = (gl[geno] > SMALLEST_GL) ? gl[geno] : SMALLEST_GL;
    }

    // Phred-scaled genotype likelihoods
    uint32_t pl[MAX_CN+1];
    uint32_t plSum = 0;
    FLP denominator = 0;
    for(unsigned int geno=0; geno<=MAX_CN; ++geno) {
      pl[geno] = (uint32_t) boost::math::round(-10 * gl[geno]);
      plSum += pl[geno];
      denominator += bl.phred2prob[pl[geno]];
    }
    if (plSum > 0) {
      FLP likelihood = (FLP) std::log10(1-1/denominator);
      likelihood = (likelihood > SMALLEST_GL) ? likelihood : SMALLEST_GL;
      gqval[0] = (int32_t) boost::math::round(-10 * likelihood);
    } else gqval[0] = 0;
    gts[0] = bcf_gt_missing;
    gts[1] = bcf_gt_missing;
    for(unsigned int geno=0; geno<=MAX_CN; ++geno) gls[geno] = (float) gl[geno];
  }

  
  template<typename TConfig>
  inline void
  genotypeCNVs(TConfig const& c, SDAggregator const& sda, std::vector<CNV> const& cnvCalls) {
    // BoLog class
    BoLog<double> bl;
    
    // Open one bam file header
    samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* bamhd = sam_hdr_read(samfile);
  
    // Output all structural variants
    std::string filename = c.outprefix + ".bcf";
    htsFile *fp = hts_open(filename.c_str(), "wb");
    bcf_hdr_t *hdr = bcf_hdr_init("w");

    // Print vcf header
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    boost::gregorian::date today = now.date();
    std::string datestr("##fileDate=");
    datestr += boost::gregorian::to_iso_string(today);
    bcf_hdr_append(hdr, datestr.c_str());
    bcf_hdr_append(hdr, "##ALT=<ID=CNV,Description=\"Copy-number variant\">");
    bcf_hdr_append(hdr, "##FILTER=<ID=LowQual,Description=\"Low quality CNV call.\">");
    bcf_hdr_append(hdr, "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END\">");
    bcf_hdr_append(hdr, "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS\">");
    bcf_hdr_append(hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the CNV\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SRL,Number=1,Type=Integer,Description=\"Split-read support left breakpoint\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SRR,Number=1,Type=Integer,Description=\"Split-read support right breakpoint\">");
    bcf_hdr_append(hdr, "##INFO=<ID=PTY,Number=1,Type=Float,Description=\"CN penalty of CNV and boundary regions\">");
    bcf_hdr_append(hdr, "##INFO=<ID=MPB,Number=1,Type=Float,Description=\"Mappable fraction\">");
    bcf_hdr_append(hdr, "##INFO=<ID=RDS,Number=1,Type=Float,Description=\"Read-depth increase/decrease support\">");
    bcf_hdr_append(hdr, "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">");
    bcf_hdr_append(hdr, "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect SV\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=CNL,Number=G,Type=Float,Description=\"Log10-scaled copy-number likelihoods for CN0, CN1, CN2, ..., CN9\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Per-sample genotype filter\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Read-depth based copy-number estimate for autosomal sites\">");
    // Add reference
    std::string refloc("##reference=");
    refloc += c.genome.string();
    bcf_hdr_append(hdr, refloc.c_str());
    for (int i = 0; i<bamhd->n_targets; ++i) {
      std::string refname("##contig=<ID=");
      refname += std::string(bamhd->target_name[i]) + ",length=" + boost::lexical_cast<std::string>(bamhd->target_len[i]) + ">";
      bcf_hdr_append(hdr, refname.c_str());
    }
    // Add samples
    bcf_hdr_add_sample(hdr, c.sampleName.c_str());
    bcf_hdr_add_sample(hdr, NULL);
    bcf_hdr_write(fp, hdr);
    
    if (!cnvCalls.empty()) {
      // Genotype arrays
      int32_t *gts = (int*) malloc(bcf_hdr_nsamples(hdr) * 2 * sizeof(int));
      float *gls = (float*) malloc(bcf_hdr_nsamples(hdr) * (MAX_CN + 1) * sizeof(float));
      int32_t *cnest = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
      int32_t *gqval = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
      std::vector<std::string> ftarr;
      ftarr.resize(bcf_hdr_nsamples(hdr));

      // Iterate all structural variants
      bcf1_t *rec = bcf_init();
      for(uint32_t i = 0; i < cnvCalls.size(); ++i) {
	// Output main vcf fields
	int32_t tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");
	//tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "LowQual");
	rec->rid = bcf_hdr_name2id(hdr, bamhd->target_name[cnvCalls[i].chr]);
	int32_t svStartPos = cnvCalls[i].start;
	if (svStartPos < 1) svStartPos = 1;
	int32_t svEndPos = cnvCalls[i].end;
	if (svEndPos < 1) svEndPos = 1;
	if (svEndPos >= (int32_t) bamhd->target_len[cnvCalls[i].chr]) svEndPos = bamhd->target_len[cnvCalls[i].chr] - 1;
	rec->pos = svStartPos;
	std::string id("CNV");
	std::string padNumber = boost::lexical_cast<std::string>(i + 1);
	padNumber.insert(padNumber.begin(), 8 - padNumber.length(), '0');
	id += padNumber;
	bcf_update_id(hdr, rec, id.c_str());
	std::string alleles = "N,<CNV>";
	bcf_update_alleles_str(hdr, rec, alleles.c_str());
	bcf_update_filter(hdr, rec, &tmpi, 1);
      
	// Add INFO fields
	if ((cnvCalls[i].srleft > 0) && (cnvCalls[i].srright > 0)) bcf_update_info_flag(hdr, rec, "PRECISE", NULL, 1);
	else bcf_update_info_flag(hdr, rec, "IMPRECISE", NULL, 1);	
	std::string svt("CNV");
	bcf_update_info_string(hdr, rec, "SVTYPE", svt.c_str());
	std::string coralVersion("EMBL.CORALv");
	coralVersion += coralVersionNumber;
	bcf_update_info_string(hdr,rec, "SVMETHOD", coralVersion.c_str());
	tmpi = svEndPos;
	bcf_update_info_int32(hdr, rec, "END", &tmpi, 1);
	tmpi = cnvCalls[i].srleft;
	bcf_update_info_int32(hdr, rec, "SRL", &tmpi, 1);
	tmpi = cnvCalls[i].srright;
	bcf_update_info_int32(hdr, rec, "SRR", &tmpi, 1);
	int32_t ciend[2];
	ciend[0] = cnvCalls[i].ciendlow - cnvCalls[i].end;
	if (ciend[0] >= 0) ciend[0] = -1;
	ciend[1] = cnvCalls[i].ciendhigh - cnvCalls[i].end;
	if (ciend[1] <= 0) ciend[1] = 1;
	int32_t cipos[2];
	cipos[0] = cnvCalls[i].ciposlow - cnvCalls[i].start;
	if (cipos[0] >= 0) cipos[0] = -1;
	cipos[1] = cnvCalls[i].ciposhigh - cnvCalls[i].start;
	if (cipos[1] <= 0) cipos[1] = 1;
	bcf_update_info_int32(hdr, rec, "CIPOS", cipos, 2);
	bcf_update_info_int32(hdr, rec, "CIEND", ciend, 2);
      	float tmpf = cnvCalls[i].penalty;
	bcf_update_info_float(hdr, rec, "PTY", &tmpf, 1);
	tmpf = cnvCalls[i].mappable;
	bcf_update_info_float(hdr, rec, "MPB", &tmpf, 1);
	tmpf = cnvCalls[i].rdsupport;
	bcf_update_info_float(hdr, rec, "RDS", &tmpf, 1);
      
	// Add genotype columns
	int32_t cnEstimate = (int32_t) boost::math::lround(cnvCalls[i].cn);
	cnest[0] = cnEstimate;

	// Compute GLs
	_computeGLs(bl, sda, cnvCalls[i], gls, gqval, gts);
	
	// Genotype filter
	if (gqval[0] < 15) ftarr[0] = "LowQual";
	else ftarr[0] = "PASS";
	rec->qual = gqval[0];
	bcf_update_genotypes(hdr, rec, gts, bcf_hdr_nsamples(hdr) * 2);
	bcf_update_format_float(hdr, rec, "CNL",  gls, bcf_hdr_nsamples(hdr) * (MAX_CN + 1));
	bcf_update_format_int32(hdr, rec, "GQ", gqval, bcf_hdr_nsamples(hdr));
	std::vector<const char*> strp(bcf_hdr_nsamples(hdr));
	std::transform(ftarr.begin(), ftarr.end(), strp.begin(), cstyle_str());
	bcf_update_format_string(hdr, rec, "FT", &strp[0], bcf_hdr_nsamples(hdr));
	bcf_update_format_int32(hdr, rec, "CN", cnest, bcf_hdr_nsamples(hdr));
	bcf_write1(fp, hdr, rec);
	bcf_clear1(rec);
      }
      bcf_destroy1(rec);
    
      // Clean-up
      free(gts);
      free(gls);
      free(cnest);
      free(gqval);
    }
  
    // Close BAM file
    bam_hdr_destroy(bamhd);
    sam_close(samfile);
    
    // Close VCF file
    bcf_hdr_destroy(hdr);
    hts_close(fp);
  
    // Build index
    bcf_index_build(filename.c_str(), 14);
  }
  

}

#endif
