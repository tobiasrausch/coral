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

#include "bed.h"
#include "version.h"
#include "scan.h"
#include "util.h"


namespace coralns
{

  struct CountDNAConfig {
    bool hasStatsFile;
    bool hasBedFile;
    bool hasScanFile;
    bool noScanWindowSelection;
    uint32_t nchr;
    uint32_t meanisize;
    uint32_t window_size;
    uint32_t window_offset;
    uint32_t scanWindow;
    uint32_t minChrLen;
    uint32_t minCnvSize;
    uint16_t minQual;
    uint16_t mad;
    float exclgc;
    float uniqueToTotalCovRatio;
    float fracWindow;
    float fragmentUnique;
    std::string sampleName;
    boost::filesystem::path genome;
    boost::filesystem::path statsFile;
    boost::filesystem::path mapFile;
    boost::filesystem::path bamFile;
    boost::filesystem::path bedFile;
    boost::filesystem::path scanFile;
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

    // BED regions
    typedef std::set<std::pair<uint32_t, uint32_t> > TChrIntervals;
    typedef std::vector<TChrIntervals> TRegionsGenome;
    TRegionsGenome bedRegions;
    if (c.hasBedFile) {
      if (!_parsePotOverlappingIntervals(c.bedFile.string(), c.hasBedFile, hdr, bedRegions)) {
	std::cerr << "Couldn't parse BED intervals. Do the chromosome names match?" << std::endl;
	return 1;
      }
    }
    
    // Parse BAM file
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Count fragments" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );

    // Open output file
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
    dataOut << "chr\tstart\tend\t" << c.sampleName << "_counts\t" << c.sampleName << "_CN" << std::endl;
    
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

	// Sum across fragment
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

      if (c.hasBedFile) {
	for(typename TChrIntervals::iterator it = bedRegions[refIndex].begin(); it != bedRegions[refIndex].end(); ++it) {
	  if ((it->first < it->second) && (it->second < hdr->target_len[refIndex])) {
	    double covsum = 0;
	    double expcov = 0;
	    double obsexp = 0;
	    uint32_t winlen = 0;
	    for(uint32_t pos = it->first; pos < it->second; ++pos) {
	      if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
		covsum += cov[pos];
		obsexp += gcbias[gcContent[pos]].obsexp;
		expcov += gcbias[gcContent[pos]].coverage;
		++winlen;
	      }
	    }
	    if (2 * winlen > (it->second - it->first)) {
	      obsexp /= (double) winlen;
	      double count = ((double) covsum / obsexp ) * (double) (it->second - it->first) / (double) winlen;
	      double cn = 2 * covsum / expcov;
	      dataOut << std::string(hdr->target_name[refIndex]) << "\t" << it->first << "\t" << it->second << "\t" << count << "\t" << cn << std::endl;
	    } else {
	      dataOut << std::string(hdr->target_name[refIndex]) << "\t" << it->first << "\t" << it->second << "\tNA\tNA" << std::endl;
	    }
	  }
	}
      } else {
	// Find #mappable pos
	uint32_t mappable = 0;
	for(uint32_t pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
	  if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	    ++mappable;
	  }
	}
	typedef std::vector<TCount> TCoverage;
	TCoverage mapcov(mappable, 0);
	std::vector<uint16_t> mapGcContent(mappable, 0);
	{
	  uint32_t mapIdx = 0;
	  for(uint32_t pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
	    if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
	      mapcov[mapIdx] = cov[pos];
	      mapGcContent[mapIdx] = gcContent[pos];
	      ++mapIdx;
	    }
	  }
	}

	// Find breakpoints
	typedef boost::dynamic_bitset<> TBitSet;
	TBitSet bp(mappable, false);

	// Iterate a couple CNV sizes to find breakpoint regions
	uint32_t boundarySize = 10000;
	for(uint32_t cnvSize = c.minCnvSize; cnvSize <= 1000000; cnvSize = (uint32_t) (cnvSize * 1.3)) {
	  double covsumLeft = 0;
	  double expcovLeft = 0;
	  double covsumMiddle = 0;
	  double expcovMiddle = 0;
	  double covsumRight = 0;
	  double expcovRight = 0;
	  for(int32_t pos = boundarySize; pos < (int32_t) mappable - boundarySize - cnvSize; ++pos) {
	    if (pos == boundarySize) {
	      for(int32_t i = 0; i<boundarySize; ++i) {
		covsumLeft += mapcov[i];
		expcovLeft += gcbias[mapGcContent[i]].coverage;
	      }
	      for(int32_t i = boundarySize; i < boundarySize + cnvSize; ++i) {
		covsumMiddle += mapcov[i];
		expcovMiddle += gcbias[mapGcContent[i]].coverage;
	      }
	      for(int32_t i = boundarySize + cnvSize; i < boundarySize + cnvSize + boundarySize; ++i) {
		covsumRight += mapcov[i];
		expcovRight += gcbias[mapGcContent[i]].coverage;
	      }
	    } else {
	      covsumLeft -= mapcov[pos - boundarySize - 1];
	      expcovLeft -= gcbias[mapGcContent[pos - boundarySize - 1]].coverage;
	      covsumMiddle -= mapcov[pos - 1];
	      expcovMiddle -= gcbias[mapGcContent[pos - 1]].coverage;
	      covsumRight -= mapcov[pos + cnvSize - 1];
	      expcovRight -= gcbias[mapGcContent[pos + cnvSize - 1]].coverage;
	      covsumLeft += mapcov[pos];
	      expcovLeft += gcbias[mapGcContent[pos]].coverage;
	      covsumMiddle += mapcov[pos + cnvSize];
	      expcovMiddle += gcbias[mapGcContent[pos + cnvSize]].coverage;
	      covsumRight += mapcov[pos + cnvSize + boundarySize];
	      expcovRight += gcbias[mapGcContent[pos + cnvSize + boundarySize]].coverage;
	    }
	    double cnLeft = 2 * covsumLeft / expcovLeft;
	    double cnMiddle = 2 * covsumMiddle / expcovMiddle;
	    double cnRight = 2 * covsumRight / expcovRight;
	    if ((std::abs(cnLeft - cnRight) < 0.2) && (std::abs(cnLeft - cnMiddle) > 0.8) && (std::abs(cnRight - cnMiddle) > 0.8) && ((cnMiddle < 1.5) || (cnMiddle > 2.5))) {
	      bp[pos + cnvSize / 2] = true;
	    }
	  }
	}

	// Get start and end of CNVs
	int32_t svStart = -1;
	int32_t svEnd = -1;
	for(uint32_t i = 0; i < mappable; ++i) {
	  if (bp[i]) {
	    if (svStart != -1) svEnd = i;
	    else {
	      svStart = i;
	      svEnd = i;
	    }
	  } else {
	    if (svStart != -1) {
	      if ((svStart < svEnd) && (svEnd - svStart > 50)) {
		// Start from core region
		int32_t offset = (int32_t) (0.2 * (svEnd - svStart));
		svStart += offset;
		svEnd -= offset;
		if (svStart < svEnd) {
		  double covsum = 0;
		  double expcov = 0;
		  for(int32_t pos = svStart; pos <= svEnd; ++pos) {
		    covsum += mapcov[pos];
		    expcov += gcbias[mapGcContent[pos]].coverage;
		  }
		  double cn = 2 * covsum / expcov;

		  // Use x-drop extension
		  double runningCn = cn;
		  for(int32_t extension = 0; std::abs(runningCn - cn) < 0.2; ++extension) {
		    double runningCnLeft = -1;
		    double runningCnRight = -1;
		    if (svStart > 0) {
		      double covsumLeft = covsum + mapcov[svStart - 1];
		      double expcovLeft = expcov + gcbias[mapGcContent[svStart - 1]].coverage;
		      runningCnLeft =  2 * covsumLeft / expcovLeft;
		    }
		    if (svEnd + 1 < mappable) {
		      double covsumRight = covsum + mapcov[svEnd + 1];
		      double expcovRight = expcov + gcbias[mapGcContent[svEnd + 1]].coverage;
		      runningCnRight = 2 * covsumRight / expcovRight;
		    }
		    if ((runningCnLeft != -1) && (runningCnRight != -1)) {
		      if (std::abs(runningCnLeft - cn) < (runningCnRight - cn)) {
			--svStart;
			covsum += mapcov[svStart];
			expcov += gcbias[mapGcContent[svStart]].coverage;
		      } else {
			++svEnd;
			covsum += mapcov[svEnd];
			expcov += gcbias[mapGcContent[svEnd]].coverage;
		      }
		    } else if (runningCnLeft != -1) {
		      --svStart;
		      covsum += mapcov[svStart];
		      expcov += gcbias[mapGcContent[svStart]].coverage;
		    } else if (runningCnRight != -1) {
		      ++svEnd;
		      covsum += mapcov[svEnd];
		      expcov += gcbias[mapGcContent[svEnd]].coverage;
		    } else break;
		    runningCn = 2 * covsum / expcov;
		  }

		  // Remap position
		  if ((svEnd - svStart) >= c.minCnvSize) {
		    uint32_t mapIdx = 0;
		    int32_t svStartHg = 0;
		    int32_t svEndHg = 0;
		    for(uint32_t pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
		      if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
			if (mapIdx == svStart) svStartHg = pos;
			if (mapIdx == svEnd) svEndHg = pos;
			++mapIdx;
		      }
		    }
		    std::cerr << svStartHg << '\t' << svEndHg << '\t' << cn << std::endl;
		  }
		}
	      }
	      svStart = -1;
	      svEnd = -1;
	    }
	  }
	}
	
	/*
	// Use windows
	for(uint32_t start = 0; start < hdr->target_len[refIndex]; start = start + c.window_offset) {
	  if (start + c.window_size < hdr->target_len[refIndex]) {
	    double covsum = 0;
	    double expcov = 0;
	    double obsexp = 0;
	    uint32_t winlen = 0;
	    for(uint32_t pos = start; pos < start + c.window_size; ++pos) {
	      if ((gcContent[pos] > gcbound.first) && (gcContent[pos] < gcbound.second) && (uniqContent[pos] >= c.fragmentUnique * c.meanisize)) {
		covsum += cov[pos];
		obsexp += gcbias[gcContent[pos]].obsexp;
		expcov += gcbias[gcContent[pos]].coverage;
		++winlen;
	      }
	    }
	    if (winlen >= c.fracWindow * c.window_size) {
	      obsexp /= (double) winlen;
	      double count = ((double) covsum / obsexp ) * (double) c.window_size / (double) winlen;
	      double cn = 2 * covsum / expcov;
	      dataOut << std::string(hdr->target_name[refIndex]) << "\t" << start << "\t" << (start + c.window_size) << "\t" << count << "\t" << cn << std::endl;
	    }
	  }
	}
	*/
      }
    }
	  
    // clean-up
    fai_destroy(faiRef);
    fai_destroy(faiMap);
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
    dataOut.pop();

    // Done
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Done." << std::endl;
    return 0;
  }

  
  int countReads(int argc, char **argv) {
    CountDNAConfig c;

    //cov.q10.d3.p0.0005.f0.8.e0.95.k0.5.gz
    //cov.q10.d3.p0.0005.f0.8.e0.97.k0.5.gz

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome file")
      ("sample,s", boost::program_options::value<std::string>(&c.sampleName)->default_value("NA12878"), "sample name")
      ("quality,q", boost::program_options::value<uint16_t>(&c.minQual)->default_value(10), "min. mapping quality")
      ("mappability,m", boost::program_options::value<boost::filesystem::path>(&c.mapFile), "input mappability map")
      ("minsize,z", boost::program_options::value<uint32_t>(&c.minCnvSize)->default_value(250), "min. CNV size")
      ("fragment,e", boost::program_options::value<float>(&c.fragmentUnique)->default_value(0.97), "min. fragment uniqueness [0,1]")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("cov.gz"), "coverage output file")
      ;

    boost::program_options::options_description window("Window options");
    window.add_options()
      ("window-size,i", boost::program_options::value<uint32_t>(&c.window_size)->default_value(10000), "window size")
      ("window-offset,j", boost::program_options::value<uint32_t>(&c.window_offset)->default_value(10000), "window offset")
      ("bed-intervals,b", boost::program_options::value<boost::filesystem::path>(&c.bedFile), "input BED file")
      ("fraction-window,k", boost::program_options::value<float>(&c.fracWindow)->default_value(0.5), "min. callable window fraction [0,1]")
      ;

    boost::program_options::options_description gcopt("GC options");
    gcopt.add_options()
      ("scan-window,c", boost::program_options::value<uint32_t>(&c.scanWindow)->default_value(10000), "scanning window size")
      ("fraction-unique,f", boost::program_options::value<float>(&c.uniqueToTotalCovRatio)->default_value(0.8), "uniqueness filter for scan windows [0,1]")
      ("scan-regions,r", boost::program_options::value<boost::filesystem::path>(&c.scanFile), "scanning regions in BED format")
      ("mad-cutoff,d", boost::program_options::value<uint16_t>(&c.mad)->default_value(3), "median + 3 * mad count cutoff")
      ("percentile,p", boost::program_options::value<float>(&c.exclgc)->default_value(0.0005), "excl. extreme GC fraction")
      ("no-window-selection,n", "no scan window selection")
      ;      

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("statsfile,t", boost::program_options::value<boost::filesystem::path>(&c.statsFile), "stats output file")
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
      std::cout << "Usage: coral " << argv[0] << " [OPTIONS] <aligned.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    }

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "coral ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;

    // Stats file
    if (vm.count("statsfile")) c.hasStatsFile = true;
    else c.hasStatsFile = false;

    
    // BED intervals
    if (vm.count("bed-intervals")) c.hasBedFile = true;
    else c.hasBedFile = false;

    // Scan regions
    if (vm.count("scan-regions")) c.hasScanFile = true;
    else c.hasScanFile = false;

    // Scan window selection
    if (vm.count("no-window-selection")) c.noScanWindowSelection = true;
    else c.noScanWindowSelection = false;
    
    // Open stats file
    boost::iostreams::filtering_ostream statsOut;
    if (c.hasStatsFile) {
      statsOut.push(boost::iostreams::gzip_compressor());
      statsOut.push(boost::iostreams::file_sink(c.statsFile.string().c_str(), std::ios_base::out | std::ios_base::binary));
    }

    // Library info
    LibraryInfo li;

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

      // Estimate insert size
      getLibraryParams(c, li);
      c.meanisize = ((int32_t) (li.median / 2)) * 2 + 1;
      if (c.hasStatsFile) {
	statsOut << "LP\t" << li.rs << ',' << li.median << ',' << li.mad << ',' << li.minNormalISize << ',' << li.maxNormalISize << std::endl;
      }

      // Clean-up
      bam_hdr_destroy(hdr);
      hts_idx_destroy(idx);
      sam_close(samfile);
    }

    // Scan genomic windows
    typedef std::vector<ScanWindow> TWindowCounts;
    typedef std::vector<TWindowCounts> TGenomicWindowCounts;
    TGenomicWindowCounts scanCounts(c.nchr, TWindowCounts());
    scan(c, li, scanCounts);
    
    // Select stable windows
    selectWindows(c, scanCounts);
    if (c.hasStatsFile) {
      samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
      bam_hdr_t* hdr = sam_hdr_read(samfile);
      statsOut << "SW\tchrom\tstart\tend\tselected\tcoverage\tuniqcov" <<  std::endl;
      for(uint32_t refIndex = 0; refIndex < (uint32_t) hdr->n_targets; ++refIndex) {
	for(uint32_t i = 0; i < scanCounts[refIndex].size(); ++i) {
	  statsOut << "SW\t" <<  hdr->target_name[refIndex] << '\t' << scanCounts[refIndex][i].start << '\t' << scanCounts[refIndex][i].end << '\t' << scanCounts[refIndex][i].select << '\t' << scanCounts[refIndex][i].cov << '\t' << scanCounts[refIndex][i].uniqcov << std::endl;
	}
      }
      bam_hdr_destroy(hdr);
      sam_close(samfile);
    }
    
    // GC bias estimation
    std::vector<GcBias> gcbias(c.meanisize + 1, GcBias());
    gcBias(c, scanCounts, li, gcbias);
    if (c.hasStatsFile) {
      statsOut << "GC\tgcsum\tsample\treference\tpercentileSample\tpercentileReference\tfractionSample\tfractionReference\tobsexp\tmeancoverage" << std::endl;
      for(uint32_t i = 0; i < gcbias.size(); ++i) statsOut << "GC\t" << i << "\t" << gcbias[i].sample << "\t" << gcbias[i].reference << "\t" << gcbias[i].percentileSample << "\t" << gcbias[i].percentileReference << "\t" << gcbias[i].fractionSample << "\t" << gcbias[i].fractionReference << "\t" << gcbias[i].obsexp << "\t" << gcbias[i].coverage << std::endl;
    }

    // Correctable GC range
    typedef std::pair<uint32_t, uint32_t> TGCBound;
    TGCBound gcbound = gcBound(c, gcbias);
    if (c.hasStatsFile) {
      statsOut << "BoundsGC\t" << gcbound.first << "," << gcbound.second << std::endl;
    }

    // Close stats file
    if (c.hasStatsFile) {
      statsOut.pop();
    }
    
    // Count reads
    return bamCount(c, li, gcbias, gcbound);
  }

  
}

#endif
