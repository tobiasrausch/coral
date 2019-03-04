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


  struct BafConfig {
    bool tumorOnly;
    uint16_t minMapQual;
    uint16_t minBaseQual;
    boost::filesystem::path genome;
    boost::filesystem::path mapFile;
    boost::filesystem::path outfile;
    boost::filesystem::path vcffile;
    std::vector<boost::filesystem::path> files;
  };

  template<typename TConfig>
  inline int32_t
  bafRun(TConfig& c) {
  /*

    // Load bam files
    samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.bamfile.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Load bcf file
    htsFile* ibcffile = bcf_open(c.vcffile.string().c_str(), "r");
    hts_idx_t* bcfidx = bcf_index_load(c.vcffile.string().c_str());
    bcf_hdr_t* bcfhdr = bcf_hdr_read(ibcffile);
    
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Assign reads to haplotypes" << std::endl;
    boost::progress_display show_progress(hdr->n_targets);

    // Allele support file
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.as.string().c_str(), std::ios_base::out | std::ios_base::binary));
    dataOut << "chr\tpos\tid\tref\talt\tdepth\trefsupport\taltsupport\tgt\taf\tpvalue" << std::endl;
  
    // Assign reads to SNPs
    faidx_t* fai = fai_load(c.genome.string().c_str());
    for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
      std::string chrName(hdr->target_name[refIndex]);
      ++show_progress;

      // Load het. markers
      typedef std::vector<BiallelicVariant> TPhasedVariants;
      TPhasedVariants pv;
      if (!_loadVariants(ibcffile, bcfidx, bcfhdr, c.sample, chrName, pv)) continue;
      if (pv.empty()) continue;

      // Sort variants
      std::sort(pv.begin(), pv.end(), SortVariants<BiallelicVariant>());

      // Load reference
      int32_t seqlen = -1;
      char* seq = NULL;
      seq = faidx_fetch_seq(fai, chrName.c_str(), 0, hdr->target_len[refIndex], &seqlen);
      
      // Annotate REF and ALT support
      typedef std::vector<uint32_t> TAlleleSupport;
      TAlleleSupport ref(pv.size(), 0);
      TAlleleSupport alt(pv.size(), 0);
      hts_itr_t* itr = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* r = bam_init1();
      while (sam_itr_next(samfile, itr, r) >= 0) {
	if (r->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	if ((r->core.qual < c.minMapQual) || (r->core.tid<0)) continue;
	if ((r->core.flag & BAM_FPAIRED) && (r->core.flag & BAM_FMUNMAP)) continue;

	// Fetch contained variants
	TPhasedVariants::const_iterator vIt = std::lower_bound(pv.begin(), pv.end(), BiallelicVariant(r->core.pos), SortVariants<BiallelicVariant>());
	TPhasedVariants::const_iterator vItEnd = std::upper_bound(pv.begin(), pv.end(), BiallelicVariant(lastAlignedPosition(r)), SortVariants<BiallelicVariant>());
	if (vIt != vItEnd) {
	  // Get read sequence
	  std::string sequence;
	  sequence.resize(r->core.l_qseq);
	  uint8_t* seqptr = bam_get_seq(r);
	  for (int32_t i = 0; i < r->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

	  // Get base qualities
	  typedef std::vector<uint8_t> TQuality;
	  TQuality quality;
	  quality.resize(r->core.l_qseq);
	  uint8_t* qualptr = bam_get_qual(r);
	  for (int i = 0; i < r->core.l_qseq; ++i) quality[i] = qualptr[i];
	  
	  // Parse CIGAR
	  uint32_t* cigar = bam_get_cigar(r);
	  for(;vIt != vItEnd; ++vIt) {
	    int32_t gp = r->core.pos; // Genomic position
	    int32_t sp = 0; // Sequence position
	    bool varFound = false;
	    for (std::size_t i = 0; ((i < r->core.n_cigar) && (!varFound)); ++i) {
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
		      varFound = true;
		      if (quality[sp] >= c.minBaseQual) {
			// Check REF allele
			if (vIt->ref == std::string(seq + gp, seq + gp + vIt->ref.size())) {
			  // Check ALT allele
			  if ((sp + vIt->alt.size() < sequence.size()) && (sp + vIt->ref.size() < sequence.size())) {
			    if (vIt->ref.size() == vIt->alt.size()) {
			      // SNP
			      if ((sequence.substr(sp, vIt->alt.size()) == vIt->alt) && (sequence.substr(sp, vIt->ref.size()) != vIt->ref)) {
				++alt[vIt-pv.begin()];
			      } else if ((sequence.substr(sp, vIt->alt.size()) != vIt->alt) && (sequence.substr(sp, vIt->ref.size()) == vIt->ref)) {
				++ref[vIt-pv.begin()];
			      }
			    }
			  } else if (vIt->ref.size() < vIt->alt.size()) {
			    // Insertion
			    int32_t diff = vIt->alt.size() - vIt->ref.size();
			    std::string refProbe = vIt->ref + std::string(seq + gp + vIt->ref.size(), seq + gp + vIt->ref.size() + diff);
			    if ((sequence.substr(sp, vIt->alt.size()) == vIt->alt) && (sequence.substr(sp, vIt->alt.size()) != refProbe)) {
			      ++alt[vIt-pv.begin()];
			    } else if ((sequence.substr(sp, vIt->alt.size()) != vIt->alt) && (sequence.substr(sp, vIt->alt.size()) == refProbe)) {
			      ++ref[vIt-pv.begin()];
			    }
			  } else {
			    // Deletion
			    int32_t diff = vIt->ref.size() - vIt->alt.size();
			    std::string altProbe = vIt->alt + std::string(seq + gp + vIt->ref.size(), seq + gp + vIt->ref.size() + diff);
			    if ((sequence.substr(sp, vIt->ref.size()) == altProbe) && (sequence.substr(sp, vIt->ref.size()) != vIt->ref)) {
			      ++alt[vIt-pv.begin()];
			    } else if ((sequence.substr(sp, vIt->ref.size()) != altProbe) && (sequence.substr(sp, vIt->ref.size()) == vIt->ref)) {
			      ++ref[vIt-pv.begin()];
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	      else {
		std::cerr << "Unknown Cigar options" << std::endl;
		return 1;
	      }
	    }
	  }
	}
      }
      bam_destroy1(r);
      hts_itr_destroy(itr);
      if (seqlen) free(seq);

      // Output (phased) allele support
      hts_itr_t* itervcf = bcf_itr_querys(bcfidx, bcfhdr, chrName.c_str());
      if (itervcf != NULL) {
	bcf1_t* recvcf = bcf_init1();
	for (uint32_t i = 0; i<pv.size(); ++i) {
	  // Fetch variant annotation from VCF
	  int32_t itrRet = 0;
	  do {
	    itrRet = bcf_itr_next(ibcffile, itervcf, recvcf);
	    if (itrRet >= 0) {
	      bcf_unpack(recvcf, BCF_UN_SHR);
	      std::vector<std::string> alleles;
	      for(std::size_t k = 0; k<recvcf->n_allele; ++k) alleles.push_back(std::string(recvcf->d.allele[k]));
	      if ((recvcf->pos == pv[i].pos) && (pv[i].ref == alleles[0]) && (pv[i].alt == alleles[1])) break;
	    } else {
	      std::cerr << "Error: Variant not found! " << chrName << ":" << (pv[i].pos + 1) << std::endl;
	      return 1;
	    }
	  } while (itrRet >= 0);
	  uint32_t totalcov = ref[i] + alt[i];
	  std::string hapstr = "0/1";
	  if (c.isPhased) {
	    if (pv[i].hap) hapstr = "1|0";
	    else hapstr = "0|1";
	  }
	  if (totalcov > 0) {
	    double h1af = 0;
	    double vaf = (double) alt[i] / (double) totalcov;
	    if (pv[i].hap) h1af = (double) alt[i] / (double) totalcov;
	    else h1af = (double) ref[i] / (double) totalcov;
	    double pval = binomTest(alt[i], totalcov, 0.5);
	    dataOut << chrName << "\t" << (pv[i].pos + 1) << "\t" << recvcf->d.id << "\t" << pv[i].ref << "\t" << pv[i].alt << "\t" << totalcov << "\t" << ref[i] << "\t" << alt[i] << "\t" << hapstr << "\t";
	    if (c.isPhased) dataOut << h1af << "\t";
	    else dataOut << vaf << "\t";
	    dataOut << pval << std::endl;
	  } else {
	    if (c.outputAll) {
	      // No coverage
	      dataOut << chrName << "\t" << (pv[i].pos + 1) << "\t" << recvcf->d.id << "\t" << pv[i].ref << "\t" << pv[i].alt << "\t" << totalcov << "\t" << ref[i] << "\t" << alt[i] << "\t" << hapstr << "\tNA\tNA" << std::endl;
	    }
	  }
	}
	bcf_destroy(recvcf);
	hts_itr_destroy(itervcf);
      }
    }
    fai_destroy(fai);

    // Close bam
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);

    // Close output allele file
    dataOut.pop();

    // Close BCF
    bcf_hdr_destroy(bcfhdr);
    hts_idx_destroy(bcfidx);
    bcf_close(ibcffile);
    
    // End
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
  */
    return 0;

  }

  int baf(int argc, char **argv) {
    BafConfig c;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("map-qual,m", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(10), "min. mapping quality")
      ("base-qual,b", boost::program_options::value<uint16_t>(&c.minBaseQual)->default_value(10), "min. base quality")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "input genome file")
      ("mappability,m", boost::program_options::value<boost::filesystem::path>(&c.mapFile), "input mappability map")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("baf.gz"), "BAF output file")
      ("vcffile,v", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "input BCF file")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<std::vector<boost::filesystem::path> >(&c.files), "input bam files")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    
    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome")) || (!vm.count("mappability")) || (!vm.count("vcffile"))) {
      std::cout << "Usage: coral " << argv[0] << " [OPTIONS] -g <genome.fa> -v <snps.bcf> -m <map.fa> <tumor.bam> <control.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    }

    // Check reference
    if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
      std::cerr << "Reference file is missing: " << c.genome.string() << std::endl;
      return 1;
    } else {
      faidx_t* fai = fai_load(c.genome.string().c_str());
      if (fai == NULL) {
	if (fai_build(c.genome.string().c_str()) == -1) {
	  std::cerr << "Fail to open genome fai index for " << c.genome.string() << std::endl;
	  return 1;
	} else fai = fai_load(c.genome.string().c_str());
      }
      fai_destroy(fai);
    }

    // Check mappability
    if (!(boost::filesystem::exists(c.mapFile) && boost::filesystem::is_regular_file(c.mapFile) && boost::filesystem::file_size(c.mapFile))) {
      std::cerr << "Mappability file is missing: " << c.mapFile.string() << std::endl;
      return 1;
    } else {
      faidx_t* fai = fai_load(c.mapFile.string().c_str());
      if (fai == NULL) {
	if (fai_build(c.mapFile.string().c_str()) == -1) {
	  std::cerr << "Fail to open mappability fai index for " << c.mapFile.string() << std::endl;
	  return 1;
	} else fai = fai_load(c.mapFile.string().c_str());
      }
      fai_destroy(fai);
    }
    
    // Check input BAM files
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      if (!(boost::filesystem::exists(c.files[file_c]) && boost::filesystem::is_regular_file(c.files[file_c]) && boost::filesystem::file_size(c.files[file_c]))) {
	std::cerr << "Alignment file is missing: " << c.files[file_c].string() << std::endl;
	return 1;
      }
      samFile* samfile = sam_open(c.files[file_c].string().c_str(), "r");
      if (samfile == NULL) {
	std::cerr << "Fail to open file " << c.files[file_c].string() << std::endl;
	return 1;
      }
      hts_idx_t* idx = sam_index_load(samfile, c.files[file_c].string().c_str());
      if (idx == NULL) {
	std::cerr << "Fail to open index for " << c.files[file_c].string() << std::endl;
	return 1;
      }
      bam_hdr_t* hdr = sam_hdr_read(samfile);
      if (hdr == NULL) {
	std::cerr << "Fail to open header for " << c.files[file_c].string() << std::endl;
	return 1;
      }
      // Do chromosomes match?
      faidx_t* fai = fai_load(c.genome.string().c_str());
      for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) {
	std::string tname(hdr->target_name[refIndex]);
	if (!faidx_has_seq(fai, tname.c_str())) {
	  std::cerr << "BAM file chromosome " << hdr->target_name[refIndex] << " is NOT present in your reference file " << c.genome.string() << std::endl;
	  return 1;
	}
      }
      fai_destroy(fai);

      // Clean-up
      bam_hdr_destroy(hdr);
      hts_idx_destroy(idx);
      sam_close(samfile);
    }
      
    // Check VCF/BCF file
    if (vm.count("vcffile")) {
      if (!(boost::filesystem::exists(c.vcffile) && boost::filesystem::is_regular_file(c.vcffile) && boost::filesystem::file_size(c.vcffile))) {
	std::cerr << "Input SNP VCF/BCF file is missing: " << c.vcffile.string() << std::endl;
	return 1;
      }
      htsFile* ifile = bcf_open(c.vcffile.string().c_str(), "r");
      if (ifile == NULL) {
	std::cerr << "Fail to open file " << c.vcffile.string() << std::endl;
	return 1;
      }
      hts_idx_t* bcfidx = bcf_index_load(c.vcffile.string().c_str());
      if (bcfidx == NULL) {
	std::cerr << "Fail to open index file for " << c.vcffile.string() << std::endl;
	return 1;
      }
      bcf_hdr_t* hdr = bcf_hdr_read(ifile);
      if (hdr == NULL) {
	std::cerr << "Fail to open header for " << c.vcffile.string() << std::endl;
	return 1;
      }
      bcf_hdr_destroy(hdr);
      hts_idx_destroy(bcfidx);
      bcf_close(ifile);
    }
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return bafRun(c);
  }


}

#endif
