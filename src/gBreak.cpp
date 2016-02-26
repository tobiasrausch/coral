/*
============================================================================
Strand-Seq Breakpoint Caller and Genotyper
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
#include <boost/dynamic_bitset.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "strandutil.h"
#include "strandsnp.h"
#include "strandhaplo.h"
#include "strandwatershed.h"
#include "strandglog.h"

using namespace streq;

struct Config {
  unsigned short minMapQual;
  uint16_t readcount;
  boost::filesystem::path wc;
  boost::filesystem::path ww;
  boost::filesystem::path outvcf;
};

template<typename TValue, typename TPosition>
inline void
_movingAverage(std::vector<TValue>& spp, TPosition const windowSize) {
  TValue movingAverage = 0;
  for(std::size_t i = 0; (i<windowSize) && (i < (TPosition) spp.size()); ++i) movingAverage += spp[i];
  for(std::size_t i = windowSize; i < spp.size() ; ++i) {
    movingAverage -= spp[i-windowSize];
    movingAverage += spp[i];
    spp[i - windowSize/2] = movingAverage / windowSize;
  }
  // Fill bounds
  for(std::size_t i = 0; i < windowSize/2; ++i) {
    spp[i] = spp[windowSize/2];
    spp[spp.size() - i - 1] = spp[spp.size() - windowSize/2 - 1];
  }
}


template<typename TConfig, typename TGenomeSegments>
inline int32_t
_findBreakpoints(TConfig const& c, TGenomeSegments& gSegments) {
  typedef typename TGenomeSegments::value_type TChrSegments;

  // Load WW bam
  samFile* samfile = sam_open(c.ww.string().c_str(), "r");
  if (samfile == NULL) {
    std::cerr << "Fail to open file " << c.ww.string() << std::endl;
    return 1;
  }
  hts_idx_t* idx = sam_index_load(samfile, c.ww.string().c_str());
  if (idx == NULL) {
    std::cerr << "Fail to open index for " << c.ww.string() << std::endl;
    return 1;
  }
  bam_hdr_t* hdr = sam_hdr_read(samfile);
  if (hdr == NULL) {
    std::cerr << "Fail to open header for " << c.ww.string() << std::endl;
    return 1;
  }

  // Set up genomic segments
  gSegments.resize(hdr->n_targets, TChrSegments());

  // Parse bam (contig by contig)
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Call Breakpoints" << std::endl;
  uint32_t genomesize = 0;
  for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) genomesize += hdr->target_len[refIndex];
  boost::progress_display show_progress( genomesize );
  for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
    show_progress += hdr->target_len[refIndex];
    // Breakpoint support
    std::vector<int32_t> midPos;
    typedef uint16_t TPeakHeight;
    typedef std::vector<TPeakHeight> TBreakpointProfile;
    TBreakpointProfile profile;

    // Running counts & position
    boost::dynamic_bitset<> pre(c.readcount);
    boost::dynamic_bitset<> post(c.readcount);
    std::vector<int32_t> pos(c.readcount);
    int32_t p = 0;
    int32_t oldPos = 0;
    uint32_t totalRC = 0;
    hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
    bam1_t* rec = bam_init1();
    while (sam_itr_next(samfile, iter, rec) >= 0) {
      if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
      if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;

      // Shift positions
      int32_t middlePos = pos[p];
      pos[p] = rec->core.pos;

      // Shift watson-crick bits
      pre[p] = post[p];
      if (rec->core.flag & BAM_FREAD1)
	if (rec->core.flag & BAM_FREVERSE) post[p] = 1;
	else post[p] = 0;
      else
	if (rec->core.flag & BAM_FREVERSE) post[p] = 0;
	else post[p] = 1;
      p = (p+1) % c.readcount;
      ++totalRC;
      
      // Get breakpoint
      if ((oldPos != middlePos) && (totalRC >= 2*c.readcount)) {
	profile.push_back(std::abs((int) pre.count() - (int) post.count()));
	midPos.push_back(middlePos);
	oldPos = middlePos;
      }
    }
    bam_destroy1(rec);
    hts_itr_destroy(iter);
    if (profile.empty()) continue;

    // Find peaks in the profile
    typedef Peak<TPeakHeight> TPeak;
    std::vector<TPeak> peaks;
    watershed(profile, peaks);
    TPeakHeight minShift = c.readcount / 3;
    TPeakHeight maxValley = c.readcount / 10;

    // Define chromosome segments
    TChrSegments chrSegments;
    chrSegments.push_back(std::make_pair(0,0));
    chrSegments.push_back(std::make_pair(profile.size() - 1, profile.size() - 1));
    for(typename std::vector<TPeak>::const_iterator it = std::lower_bound(peaks.begin(), peaks.end(), TPeak(0, 0, minShift)); it != peaks.end(); ++it) {
      int32_t leftValley = 0;
      for(int32_t i = it->maxInd; i>=0; --i)
	if (profile[i]<=maxValley) { leftValley = i; break; }
      int32_t rightValley = profile.size() - 1;
      for(int32_t i = it->maxInd; i < (int32_t) profile.size(); ++i)
	if (profile[i]<=maxValley) { rightValley = i; break; }
      // Make sure segments are non-overlapping
      bool nonover = true;
      for(typename TChrSegments::iterator itSeg = chrSegments.begin(); itSeg != chrSegments.end(); ++itSeg) {
	if ((itSeg->first <= rightValley) && (itSeg->second >= leftValley)) {
	  // Merge
	  nonover = false;
	  itSeg->first = std::min(itSeg->first, leftValley);
	  itSeg->second = std::max(itSeg->second, rightValley);
	}
      }
      if (nonover) chrSegments.push_back(std::make_pair(leftValley, rightValley));
    }
    std::sort(chrSegments.begin(), chrSegments.end());

    // Insert inner intervals
    typename TChrSegments::iterator itSegPrev = chrSegments.begin();
    if (itSegPrev != chrSegments.end()) {
      typename TChrSegments::iterator itSeg = chrSegments.begin();
      for(++itSeg; itSeg != chrSegments.end(); ++itSeg, ++itSegPrev) {
	// Parse WW file
	gSegments[refIndex].push_back(std::make_pair(midPos[itSegPrev->second], midPos[itSeg->first]));
      }
    }
  }
  // Clean-up
  bam_hdr_destroy(hdr);
  hts_idx_destroy(idx);
  sam_close(samfile);

  return 0;
}



template<typename TConfig, typename TGenomeSegments>
inline int32_t
_genotypeSegments(TConfig const& c, TGenomeSegments const& gSegments) {
  typedef typename TGenomeSegments::value_type TChrSegments;

  // Load WW & WC bam
  samFile* wwfile = sam_open(c.ww.string().c_str(), "r");
  if (wwfile == NULL) {
    std::cerr << "Fail to open file " << c.ww.string() << std::endl;
    return 1;
  }
  hts_idx_t* wwidx = sam_index_load(wwfile, c.ww.string().c_str());
  if (wwidx == NULL) {
    std::cerr << "Fail to open index for " << c.ww.string() << std::endl;
    return 1;
  }
  bam_hdr_t* hdr = sam_hdr_read(wwfile);
  if (hdr == NULL) {
    std::cerr << "Fail to open header for " << c.ww.string() << std::endl;
    return 1;
  }
  samFile* wcfile = sam_open(c.wc.string().c_str(), "r");
  if (wcfile == NULL) {
    std::cerr << "Fail to open file " << c.wc.string() << std::endl;
    return 1;
  }
  hts_idx_t* wcidx = sam_index_load(wcfile, c.wc.string().c_str());
  if (wcidx == NULL) {
    std::cerr << "Fail to open index for " << c.wc.string() << std::endl;
    return 1;
  }

  // Genotype segments
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Genotype Segments" << std::endl;
  uint32_t genomesize = 0;
  for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) genomesize += hdr->target_len[refIndex];
  boost::progress_display show_progress( genomesize );

  // Write VCF header
  htsFile *fp = hts_open(c.outvcf.string().c_str(), "wg");
  bcf_hdr_t *hdr_out = bcf_hdr_init("w");
  boost::gregorian::date today = now.date();
  std::string datestr("##fileDate=");
  datestr += boost::gregorian::to_iso_string(today);
  bcf_hdr_append(hdr_out, datestr.c_str());
  bcf_hdr_append(hdr_out, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">");
  bcf_hdr_append(hdr_out, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
  bcf_hdr_append(hdr_out, "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Log10-scaled genotype likelihoods for RR,RA,AA genotypes\">");
  bcf_hdr_append(hdr_out, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
  bcf_hdr_append(hdr_out, "##FORMAT=<ID=WT,Number=1,Type=Integer,Description=\"#Watson reads\">");
  bcf_hdr_append(hdr_out, "##FORMAT=<ID=CR,Number=1,Type=Integer,Description=\"#Crick reads\">");
  // Add reference
  for (int i = 0; i<hdr->n_targets; ++i) {
    std::string refname("##contig=<ID=");
    refname += std::string(hdr->target_name[i]) + ",length=" + boost::lexical_cast<std::string>(hdr->target_len[i]) + ">";
    bcf_hdr_append(hdr_out, refname.c_str());
  }
  // Add sample
  uint32_t numSample = 1;
  std::string smtag(_parseBamHeader(hdr));
  std::string sampleName(smtag);
  bcf_hdr_add_sample(hdr_out, sampleName.c_str());
  bcf_hdr_add_sample(hdr_out, NULL);
  bcf_hdr_write(fp, hdr_out);

  // Iterate chromosomes
  for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
    show_progress += hdr->target_len[refIndex];

    // Genotype segments
    bcf1_t* rout = bcf_init();
    int32_t* gts = (int*) malloc(numSample * 2 * sizeof(int));
    float* gls = (float*) malloc(numSample * 3 * sizeof(float));
    int32_t* gqval = (int*) malloc(numSample * sizeof(int));
    int32_t* wtcount = (int*) malloc(numSample * sizeof(int));
    int32_t* crcount = (int*) malloc(numSample * sizeof(int));
    for(typename TChrSegments::const_iterator itSeg = gSegments[refIndex].begin(); itSeg != gSegments[refIndex].end(); ++itSeg) {
      // Parse WW file
      std::vector<uint8_t> altWWQual;
      std::vector<uint8_t> refWWQual;
      hts_itr_t* iterWW = sam_itr_queryi(wwidx, refIndex, itSeg->first, itSeg->second);
      bam1_t* recWW = bam_init1();
      while (sam_itr_next(wwfile, iterWW, recWW) >= 0) {
	if (recWW->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	if ((recWW->core.qual < c.minMapQual) || (recWW->core.tid<0)) continue;

	if (recWW->core.flag & BAM_FREAD1)
	  if (recWW->core.flag & BAM_FREVERSE) altWWQual.push_back(recWW->core.qual);
	  else refWWQual.push_back(recWW->core.qual);
	else
	  if (recWW->core.flag & BAM_FREVERSE) refWWQual.push_back(recWW->core.qual);
	  else altWWQual.push_back(recWW->core.qual);
      }
      bam_destroy1(recWW);
      hts_itr_destroy(iterWW);

      // Parse WC file
      std::vector<uint8_t> altWCQual;
      std::vector<uint8_t> refWCQual;
      hts_itr_t* iterWC = sam_itr_queryi(wcidx, refIndex, itSeg->first, itSeg->second);
      bam1_t* recWC = bam_init1();
      while (sam_itr_next(wcfile, iterWC, recWC) >= 0) {
	if (recWC->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	if ((recWC->core.qual < c.minMapQual) || (recWC->core.tid<0)) continue;

	if (recWC->core.flag & BAM_FREAD1)
	  if (recWC->core.flag & BAM_FREVERSE) altWCQual.push_back(recWC->core.qual);
	  else refWCQual.push_back(recWC->core.qual);
	else
	  if (recWC->core.flag & BAM_FREVERSE) refWCQual.push_back(recWC->core.qual);
	  else altWCQual.push_back(recWC->core.qual);
      }
      bam_destroy1(recWC);
      hts_itr_destroy(iterWC);
      
      // Min. number of reads
      uint32_t wwCount = refWWQual.size() + altWWQual.size();
      uint32_t wcCount = refWCQual.size() + altWCQual.size();
      if ((wwCount) && (wcCount)) {
	int32_t gtWC = _genotype(refWCQual, altWCQual, gls, gqval, gts);
	if (gtWC != -1) {
	  int32_t gtWW = _genotype(refWWQual, altWWQual, gls, gqval, gts);
	  if (gtWW != -1) {
	    if ((gtWW) && (((gtWW == 1) && (gtWC != 1)) || ((gtWW == 2) && (gtWC==1)))) {
	      if (gtWW == 0) {
		gts[0] = bcf_gt_phased(0);
		gts[1] = bcf_gt_phased(0);
	      } else if (gtWW == 2) {
		gts[0] = bcf_gt_phased(1);
		gts[1] = bcf_gt_phased(1);
	      } else {
		if (gtWC == 0) {
		  // Watson-Watson  (inversion is on crick -> HP=2)
		  gts[0] = bcf_gt_phased(0);
		  gts[1] = bcf_gt_phased(1);
		} else {
		  gts[0] = bcf_gt_phased(1);
		  gts[1] = bcf_gt_phased(0);
		}
	      }
	      rout->rid = bcf_hdr_name2id(hdr_out, hdr->target_name[refIndex]);
	      rout->pos = itSeg->first;
	      rout->qual = gqval[0];
	      std::string id(".");
	      bcf_update_id(hdr_out, rout, id.c_str());
	      std::string alleles("N,<INV>");
	      bcf_update_alleles_str(hdr_out, rout, alleles.c_str());
	      int32_t tmpi = bcf_hdr_id2int(hdr_out, BCF_DT_ID, "PASS");
	      bcf_update_filter(hdr_out, rout, &tmpi, 1);
	      int32_t endpos = itSeg->second;
	      bcf_update_info_int32(hdr_out, rout, "END", &endpos, 1);	    
	      wtcount[0] = refWWQual.size();
	      crcount[0] = altWWQual.size();
	      bcf_update_genotypes(hdr_out, rout, gts, numSample * 2);
	      bcf_update_format_float(hdr_out, rout, "GL",  gls, numSample * 3);
	      bcf_update_format_int32(hdr_out, rout, "GQ", gqval, numSample);
	      bcf_update_format_int32(hdr_out, rout, "WT", wtcount, numSample);
	      bcf_update_format_int32(hdr_out, rout, "CR", crcount, numSample);
	      bcf_write1(fp, hdr_out, rout);
	      bcf_clear1(rout);
	    }
	  }
	}
      }
    }
    // Clean-up
    free(gts);
    free(gls);
    free(gqval);
    free(wtcount);
    free(crcount);
    bcf_destroy1(rout);
  }
  // Close VCF
  bcf_hdr_destroy(hdr_out);
  hts_close(fp);

  // Close bam
  hts_idx_destroy(wwidx);
  sam_close(wwfile);
  bam_hdr_destroy(hdr);
  hts_idx_destroy(wcidx);
  sam_close(wcfile);

  return 0;
}


int main(int argc, char **argv) {

#ifdef PROFILE
  ProfilerStart("gBreak.prof");
#endif

  Config c;

  // Parameter
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("map-qual,q", boost::program_options::value<unsigned short>(&c.minMapQual)->default_value(1), "min. mapping quality")
    ("readcount,r", boost::program_options::value<uint16_t>(&c.readcount)->default_value(250), "window length in terms of #reads (<65000)")
    ("diffstrand,d", boost::program_options::value<boost::filesystem::path>(&c.wc), "different strand bam (wc.bam)")
    ("outvcf,o", boost::program_options::value<boost::filesystem::path>(&c.outvcf)->default_value("variants.vcf.gz"), "output variants file")
    ;

  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&c.ww), "input bam file")
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
  if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("diffstrand"))) {
    std::cout << "Usage: " << argv[0] << " [OPTIONS] -d <sample.wc.bam> <sample.ww.bam>" << std::endl;
    std::cout << visible_options << "\n";
    return 1;
  } 
  
  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  // Genome segments
  typedef std::pair<int32_t, int32_t> TPeakWidth;
  typedef std::vector<TPeakWidth> TChrSegments;
  typedef std::vector<TChrSegments> TGenomeSegments;
  TGenomeSegments gSegments;
  int r = _findBreakpoints(c, gSegments);
  if (r != 0) return r;

  // Genotype segments
  r = _genotypeSegments(c, gSegments);
  if (r != 0) return r;

  // End
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;

#ifdef PROFILE
  ProfilerStop();
#endif


  return 0;
}
