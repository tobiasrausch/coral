/*
============================================================================
Single Cell Sequencing Analysis Methods
============================================================================
Copyright (C) 2018 Tobias Rausch

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

#ifndef COUNT_H
#define COUNT_H

#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <fstream>

#define BOOST_DISABLE_ASSERTS
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "util.h"
#include "json.h"

namespace sc
{

  struct CountConfig {
    unsigned short minMapQual;
    uint32_t window;
    std::string method;
    std::vector<std::string> sampleName;
    boost::filesystem::path outfile;
    std::vector<boost::filesystem::path> files;
  };
  
  
  struct WindowClassifier {
  float watsonCut;
    float crickCut;
    float median;
    float mad;
    
    WindowClassifier(): watsonCut(0), crickCut(0), median(0), mad(0) {}
    WindowClassifier(float w, float c, float m, float a) : watsonCut(w), crickCut(c), median(m), mad(a) {}
  };

  int count(int argc, char **argv) {

#ifdef PROFILE
    ProfilerStart("count.prof");
#endif
    
    CountConfig c;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("type,t", boost::program_options::value<std::string>(&c.method)->default_value("StrandSeq"), "single cell seq. method [StrandSeq|Malbac]")
      ("map-qual,q", boost::program_options::value<unsigned short>(&c.minMapQual)->default_value(1), "min. mapping quality")
      ("window,w", boost::program_options::value<uint32_t>(&c.window)->default_value(200000), "window length")
      ("outfile", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.json.gz"), "output file")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input bam file")
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
    if ((vm.count("help")) || (!vm.count("input-file"))) {
      std::cout << "Usage: sc " << argv[0] << " [OPTIONS] <sc1.bam> <sc2.bam> ... <scN.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    } 
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "sc ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    // Check BAM files
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      if (!(boost::filesystem::exists(c.files[file_c]) && boost::filesystem::is_regular_file(c.files[file_c]) && boost::filesystem::file_size(c.files[file_c]))) {
	std::cerr << "Input BAM file is missing: " << c.files[file_c].string() << std::endl;
	return 1;
      }
    }
    
    // Load bam files
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    typedef std::vector<bam_hdr_t*> THeader;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    THeader hdr(c.files.size());
    c.sampleName.resize(c.files.size());
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      if (samfile[file_c] == NULL) {
	std::cerr << "Fail to open file " << c.files[file_c].string() << std::endl;
	return 1;
      }
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
      if (idx[file_c] == NULL) {
	std::cerr << "Fail to open index for " << c.files[file_c].string() << std::endl;
	return 1;
      }
      hdr[file_c] = sam_hdr_read(samfile[file_c]);
      if (hdr[file_c] == NULL) {
	std::cerr << "Fail to open header for " << c.files[file_c].string() << std::endl;
	return 1;
      }
      std::string sampleName;
      if (!getSMTag(std::string(hdr[file_c]->text), c.files[file_c].stem().string(), sampleName)) {
	std::cerr << "Only one sample (@RG:SM) is allowed per input BAM file " << c.files[file_c].string() << std::endl;
	return 1;
      } else c.sampleName[file_c] = sampleName;
    }

    // Watson-Crick Counter
    typedef std::pair<uint32_t, uint32_t> TWatsonCrick;
    typedef std::vector<TWatsonCrick> TChrWC;
    typedef std::vector<TChrWC> TGenomicWC;
    TGenomicWC gWC(hdr[0]->n_targets, TChrWC());
    for (int refIndex = 0; refIndex<hdr[0]->n_targets; ++refIndex) {
      if (hdr[0]->target_len[refIndex] < c.window) continue;
      uint32_t bins = hdr[0]->target_len[refIndex] / c.window + 1;
      gWC[refIndex].resize(bins);
    }
    
    // Parse bam (contig by contig)
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "BAM file parsing" << std::endl;
    boost::progress_display show_progress( hdr[0]->n_targets );
    for (int refIndex = 0; refIndex<hdr[0]->n_targets; ++refIndex) {
      ++show_progress;
      if (!gWC[refIndex].size()) continue;
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr[0]->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	  if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	  if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
	
	  int32_t pos = rec->core.pos + halfAlignmentLength(rec);
	  if (rec->core.flag & BAM_FREAD1) 
	    if (rec->core.flag & BAM_FREVERSE) ++gWC[refIndex][(int) (pos / c.window)].second;
	    else ++gWC[refIndex][(int) (pos / c.window)].first;
	  else
	    if (rec->core.flag & BAM_FREVERSE) ++gWC[refIndex][(int) (pos / c.window)].first;
	    else ++gWC[refIndex][(int) (pos / c.window)].second;
	}
	bam_destroy1(rec);
	hts_itr_destroy(iter);
      }
    }

    // Output
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Output Single Cell Counts" << std::endl;
    scJsonOut(c, hdr, gWC);

    // End
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;

    // Close bam
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
      bam_hdr_destroy(hdr[file_c]);
    }

#ifdef PROFILE
    ProfilerStop();
#endif
    return 0;
  }

}

#endif
