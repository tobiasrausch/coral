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

#include <htslib/sam.h>
#include <htslib/vcf.h>



namespace coralns
{

   // Convert string to char*
  struct cstyle_str {
    const char* operator ()(const std::string& s) {
      return s.c_str();
    }
  };


  template<typename TConfig>
  inline void
  genotypeCNVs(TConfig const& c, std::vector<CNV> const& cnvCalls) {
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
    bcf_hdr_append(hdr, "##INFO=<ID=SRQ,Number=1,Type=Float,Description=\"Split-read consensus alignment quality\">");
    bcf_hdr_append(hdr, "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">");
    bcf_hdr_append(hdr, "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect SV\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Log10-scaled genotype likelihoods for RR,RA,AA genotypes\">");
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
      float *gls = (float*) malloc(bcf_hdr_nsamples(hdr) * 3 * sizeof(float));
      int32_t *cnest = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
      int32_t *gqval = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
      std::vector<std::string> ftarr;
      ftarr.resize(bcf_hdr_nsamples(hdr));


      for(uint32_t i = 0; i < cnvCalls.size(); ++i) {
	std::cerr << bamhd->target_name[cnvCalls[i].chr] << '\t' << cnvCalls[i].start << '\t' << cnvCalls[i].end << '\t' << "CN=" << cnvCalls[i].cn << ";CIPOS=" << cnvCalls[i].ciposlow << "," << cnvCalls[i].ciposhigh << ";CIEND=" << cnvCalls[i].ciendlow << "," << cnvCalls[i].ciendhigh << ",RDSUPPORT=" << cnvCalls[i].rdsupport << ",PENALTY=" << cnvCalls[i].penalty << ",MAPPABLE=" << cnvCalls[i].mappable << std::endl;
	/*
	for(uint32_t k = 0; k < splitBp.size(); ++k) {
	  if ((cnvCalls[i].ciposlow < splitBp[k].first) && (splitBp[k].first < cnvCalls[i].ciposhigh)) {
	    std::cerr << splitBp[k].first << ',' << splitBp[k].second << std::endl;
	  }
	  if ((cnvCalls[i].ciendlow < splitBp[k].first) && (splitBp[k].first < cnvCalls[i].ciendhigh)) {
	    std::cerr << splitBp[k].first << ',' << splitBp[k].second << std::endl;
	  }
	}
	*/
      }
      
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
	//bcf_update_info_flag(hdr, rec, "PRECISE", NULL, 1);
	bcf_update_info_flag(hdr, rec, "IMPRECISE", NULL, 1);
	std::string svt("CNV");
	bcf_update_info_string(hdr, rec, "SVTYPE", svt.c_str());
	std::string coralVersion("EMBL.CORALv");
	coralVersion += coralVersionNumber;
	bcf_update_info_string(hdr,rec, "SVMETHOD", coralVersion.c_str());
	tmpi = svEndPos;
	bcf_update_info_int32(hdr, rec, "END", &tmpi, 1);
	int32_t ciend[2];
	ciend[0] = cnvCalls[i].ciendlow;
	ciend[1] = cnvCalls[i].ciendhigh;
	int32_t cipos[2];
	cipos[0] = cnvCalls[i].ciposlow;
	cipos[1] = cnvCalls[i].ciposhigh;
	bcf_update_info_int32(hdr, rec, "CIPOS", cipos, 2);
	bcf_update_info_int32(hdr, rec, "CIEND", ciend, 2);
      	float tmpf = 0;
	bcf_update_info_float(hdr, rec, "SRQ", &tmpf, 1);
      
	// Add genotype columns
	cnest[0] = 0;

	// Compute GLs
	gqval[0] = 10;
	//if (svIter->precise) _computeGLs(bl, jctCountMap[file_c][svIter->id].ref, jctCountMap[file_c][svIter->id].alt, gls, gqval, gts, file_c);
	//else _computeGLs(bl, spanCountMap[file_c][svIter->id].ref, spanCountMap[file_c][svIter->id].alt, gls, gqval, gts, file_c);
	
	// Genotype filter
	if (gqval[0] < 15) ftarr[0] = "LowQual";
	else ftarr[0] = "PASS";
	rec->qual = 0;
	bcf_update_genotypes(hdr, rec, gts, bcf_hdr_nsamples(hdr) * 2);
	bcf_update_format_float(hdr, rec, "GL",  gls, bcf_hdr_nsamples(hdr) * 3);
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
