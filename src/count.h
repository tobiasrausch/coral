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

#ifndef COUNT_H
#define COUNT_H

#include <limits>

#include <boost/icl/split_interval_map.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/progress.hpp>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "version.h"
#include "scan.h"
#include "util.h"


namespace coralns
{

  struct CountDNAConfig {
    uint32_t nchr;
    uint32_t meanisize;
    uint32_t window_size;
    uint32_t window_offset;
    uint32_t scanWindow;
    uint32_t minChrLen;
    uint16_t minQual;
    uint16_t mad;
    float exclgc;
    float alignmentQ;
    std::string sampleName;
    boost::filesystem::path genome;
    boost::filesystem::path mapFile;
    boost::filesystem::path bamFile;
    boost::filesystem::path outfile;
  };

  template<typename TConfig>
  inline int32_t
  bamCount(TConfig const& c, LibraryInfo const& li, std::vector<GcBias> const& gcbias, std::pair<uint32_t, uint32_t> const& gcbound) {
    
    // Load bam file
    samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.bamFile.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Parse BAM file
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Count fragments" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );

    // Open output file
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
    dataOut << "chr\tstart\tend\t" << c.sampleName << std::endl;
    
    // Iterate chromosomes
    faidx_t* faiMap = fai_load(c.mapFile.string().c_str());
    faidx_t* faiRef = fai_load(c.genome.string().c_str());
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      ++show_progress;
      if (chrNoData(c, refIndex, idx)) continue;

      // Check presence in mappability map
      std::string tname(hdr->target_name[refIndex]);
      int32_t seqlen = faidx_seq_len(faiMap, tname.c_str());
      if (seqlen == - 1) continue;
      else seqlen = -1;
      char* seq = faidx_fetch_seq(faiMap, tname.c_str(), 0, faidx_seq_len(faiMap, tname.c_str()), &seqlen);

      // Check presence in reference
      seqlen = faidx_seq_len(faiRef, tname.c_str());
      if (seqlen == - 1) continue;
      else seqlen = -1;
      char* ref = faidx_fetch_seq(faiRef, tname.c_str(), 0, faidx_seq_len(faiRef, tname.c_str()), &seqlen);

      // Get GC and Mappability
      std::vector<uint16_t> uniqContent(hdr->target_len[refIndex], 0);
      std::vector<uint16_t> gcContent(hdr->target_len[refIndex], 0);
      {
	// Mappability map
	typedef boost::dynamic_bitset<> TBitSet;
	TBitSet uniq(hdr->target_len[refIndex], false);
	for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	  if (seq[i] == 'C') uniq[i] = 1;
	}

	// GC map
	typedef boost::dynamic_bitset<> TBitSet;
	TBitSet gcref(hdr->target_len[refIndex], false);
	for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	  if ((ref[i] == 'c') || (ref[i] == 'C') || (ref[i] == 'g') || (ref[i] == 'G')) gcref[i] = 1;
	}

	// Sum across fragments
	int32_t halfwin = (int32_t) (c.meanisize / 2);
	int32_t usum = 0;
	int32_t gcsum = 0;
	for(int32_t pos = halfwin; pos < (int32_t) hdr->target_len[refIndex] - halfwin; ++pos) {
	  if (pos == halfwin) {
	    for(int32_t i = pos - halfwin; i<=pos+halfwin; ++i) {
	      usum += uniq[i];
	      gcsum += gcref[i];
	    }
	  } else {
	    usum -= uniq[pos - halfwin - 1];
	    gcsum -= gcref[pos - halfwin - 1];
	    usum += uniq[pos + halfwin];
	    gcsum += gcref[pos + halfwin];
	  }
	  gcContent[pos] = gcsum;
	  uniqContent[pos] = usum;
	}
      }
      
      // Coverage track
      typedef uint16_t TCount;
      uint32_t maxCoverage = std::numeric_limits<TCount>::max();
      typedef std::vector<TCount> TCoverage;
      TCoverage cov(hdr->target_len[refIndex], 0);

      // Mate map
      typedef boost::unordered_map<std::size_t, bool> TMateMap;
      TMateMap mateMap;
      
      // Count reads
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      int32_t lastAlignedPos = 0;
      std::set<std::size_t> lastAlignedPosReads;
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	if ((rec->core.flag & BAM_FPAIRED) && ((rec->core.flag & BAM_FMUNMAP) || (rec->core.tid != rec->core.mtid))) continue;
	if (rec->core.qual < c.minQual) continue;

	int32_t midPoint = rec->core.pos + halfAlignmentLength(rec);
	if (rec->core.flag & BAM_FPAIRED) {
	  // Clean-up the read store for identical alignment positions
	  if (rec->core.pos > lastAlignedPos) {
	    lastAlignedPosReads.clear();
	    lastAlignedPos = rec->core.pos;
	  }
	
	  if ((rec->core.pos < rec->core.mpos) || ((rec->core.pos == rec->core.mpos) && (lastAlignedPosReads.find(hash_string(bam_get_qname(rec))) == lastAlignedPosReads.end()))) {
	    // First read
	    lastAlignedPosReads.insert(hash_string(bam_get_qname(rec)));
	    std::size_t hv = hash_pair(rec);
	    mateMap[hv] = true;
	    continue;
	  } else {
	    // Second read
	    std::size_t hv = hash_pair_mate(rec);
	    if ((mateMap.find(hv) == mateMap.end()) || (!mateMap[hv])) continue; // Mate discarded
	    mateMap[hv] = false;
	  }

	  // Insert size filter
	  int32_t isize = (rec->core.pos + alignmentLength(rec)) - rec->core.mpos;
	  if ((li.minNormalISize < isize) && (isize < li.maxNormalISize)) {
	    midPoint = rec->core.mpos + (int32_t) (isize/2);
	  } else {
	    if (rec->core.flag & BAM_FREVERSE) midPoint = rec->core.pos + alignmentLength(rec) - (c.meanisize / 2);
	    else midPoint = rec->core.pos + (c.meanisize / 2);
	  }
	}

	// Count fragment
	if ((midPoint >= 0) && (midPoint < (int32_t) hdr->target_len[refIndex]) && (cov[midPoint] < maxCoverage - 1)) ++cov[midPoint];
      }
      // Clean-up
      if (seq != NULL) free(seq);
      if (ref != NULL) free(ref);
      bam_destroy1(rec);
      hts_itr_destroy(iter);
      mateMap.clear();

      for(uint32_t start = 0; start < hdr->target_len[refIndex]; start = start + c.window_offset) {
	if (start + c.window_size < hdr->target_len[refIndex]) {
	  double covsum = 0;
	  double expcov = 0;
	  double obsexp = 0;
	  int32_t winlen = 0;
	  for(uint32_t pos = start; pos < start + c.window_size; ++pos) {
	    if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] == c.meanisize)) {
	      covsum += cov[pos];
	      obsexp += gcbias[gcContent[pos]].obsexp;
	      expcov += gcbias[gcContent[pos]].coverage;
	      ++winlen;
	    }
	  }
	  if (2 * winlen > c.window_size) {
	    obsexp /= (double) winlen;
	    double count = ((double) covsum / obsexp ) * (double) c.window_size / (double) winlen;
	    double cn = 2 * covsum / expcov;
	    dataOut << std::string(hdr->target_name[refIndex]) << "\t" << start << "\t" << (start + c.window_size) << "\t" << cn << std::endl;
	  }
	}
      }
    }
	  
    // clean-up
    fai_destroy(faiRef);
    fai_destroy(faiMap);
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
    dataOut.pop();
    
    return 0;
  }

  
  int countReads(int argc, char **argv) {
    CountDNAConfig c;
    c.meanisize = 300;   // Best guess for illumina PE
    
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome file")
      ("sample,s", boost::program_options::value<std::string>(&c.sampleName)->default_value("NA12878"), "sample name")
      ("quality,q", boost::program_options::value<uint16_t>(&c.minQual)->default_value(10), "min. mapping quality")
      ("mappability,m", boost::program_options::value<boost::filesystem::path>(&c.mapFile), "input mappability map")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("cov.gz"), "coverage output file")
      ;

    boost::program_options::options_description window("Window options");
    window.add_options()
      ("window-size,i", boost::program_options::value<uint32_t>(&c.window_size)->default_value(10000), "window size")
      ("window-offset,j", boost::program_options::value<uint32_t>(&c.window_offset)->default_value(10000), "window offset")
      ;

    boost::program_options::options_description gcopt("GC options");
    gcopt.add_options()
      ("scan-window,c", boost::program_options::value<uint32_t>(&c.scanWindow)->default_value(100000), "scanning window size")
      ("mad-cutoff,d", boost::program_options::value<uint16_t>(&c.mad)->default_value(5), "median + 5 * mad count cutoff")
      ("alignmentqual,a", boost::program_options::value<float>(&c.alignmentQ)->default_value(0), "alignment fractional identity (0:deactivate)")
      ("percentile,p", boost::program_options::value<float>(&c.exclgc)->default_value(0.005), "excl. extreme GC fraction")
      ;      

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.bamFile), "input bam file")
      ;

    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);

    // Set the visibility
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(window).add(gcopt).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(window).add(gcopt);

    // Parse command-line
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);

    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file"))) {
      std::cout << std::endl;
      std::cout << "Usage: alfred " << argv[0] << " [OPTIONS] <aligned.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    }

    // Check bam file
    if (!(boost::filesystem::exists(c.bamFile) && boost::filesystem::is_regular_file(c.bamFile) && boost::filesystem::file_size(c.bamFile))) {
      std::cerr << "Alignment file is missing: " << c.bamFile.string() << std::endl;
      return 1;
    } else {
      samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
      if (samfile == NULL) {
	std::cerr << "Fail to open file " << c.bamFile.string() << std::endl;
	return 1;
      }
      hts_idx_t* idx = sam_index_load(samfile, c.bamFile.string().c_str());
      if (idx == NULL) {
	if (bam_index_build(c.bamFile.string().c_str(), 0) != 0) {
	  std::cerr << "Fail to open index for " << c.bamFile.string() << std::endl;
	  return 1;
	}
      }
      bam_hdr_t* hdr = sam_hdr_read(samfile);
      if (hdr == NULL) {
	std::cerr << "Fail to open header for " << c.bamFile.string() << std::endl;
	return 1;
      }
      c.nchr = hdr->n_targets;
      c.minChrLen = setMinChrLen(hdr, 0.95);

      // Clean-up
      bam_hdr_destroy(hdr);
      hts_idx_destroy(idx);
      sam_close(samfile);
    }

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "coral ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;

    // Scan genomic windows
    typedef std::vector<uint32_t> TWindowCounts;
    typedef std::vector<TWindowCounts> TGenomicWindowCounts;
    TGenomicWindowCounts scanCounts(c.nchr, TWindowCounts());
    scan(c, scanCounts);
    //for(uint32_t i = 0; i < c.nchr; ++i) { for(uint32_t k = 0; k < scanCounts[i].size(); ++k) { std::cerr << scanCounts[i][k] << std::endl;  }   }

    // Estimate library parameters
    LibraryInfo li;
    getLibraryParams(c, scanCounts, li);
    c.meanisize = ((int32_t) (li.median / 2)) * 2 + 1;
    //std::cerr << li.rs << ',' << li.median << ',' << li.mad << ',' << li.minNormalISize << ',' << li.maxNormalISize << std::endl;

    // GC bias estimation
    std::vector<GcBias> gcbias(c.meanisize + 1, GcBias());
    gcBias(c, scanCounts, li, gcbias);

    // Correctable GC range
    typedef std::pair<uint32_t, uint32_t> TGCBound;
    TGCBound gcbound = gcBound(c, gcbias);
    //std::cerr << gcbound.first << "," << gcbound.second << std::endl;

    // Count reads
    return bamCount(c, li, gcbias, gcbound);
  }

  
}

#endif
