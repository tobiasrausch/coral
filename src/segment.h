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

#ifndef SEGMENT_H
#define SEGMENT_H

#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <fstream>
#include <complex>

#define BOOST_DISABLE_ASSERTS
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/multi_array.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <boost/unordered_map.hpp>

#include "config.h"
#include "matrix.h"
#include "gflars.h"
#include "dpseg.h"
#include "util.h"

namespace coralns
{

  
  inline int32_t
  runSegmentation(SegmentConfig const& c, std::vector<NormalizedBinCounts>& cnbc) {
    typedef SegmentConfig::TPrecision TPrecision;

    // Parse signal matrix
    uint32_t refIndex = 0;
    uint32_t row = 0;
    
    std::vector<TPrecision> lastVal(cnbc[refIndex].cols, 0.0);
    std::ifstream signalFile(c.signal.string().c_str(), std::ifstream::in);
    cnbc[refIndex].sm.resize(boost::extents[cnbc[refIndex].rows][cnbc[refIndex].cols]);
    if (signalFile.is_open()) {
      while (signalFile.good()) {
	std::string sigdata;
	getline(signalFile, sigdata);
	while ((row >= cnbc[refIndex].rows) && (refIndex + 1 < cnbc.size())) {
	  ++refIndex;
	  std::fill(lastVal.begin(), lastVal.end(), 0.0);
	  cnbc[refIndex].sm.resize(boost::extents[cnbc[refIndex].rows][cnbc[refIndex].cols]);
	  row = 0;
	}
	if (row < cnbc[refIndex].rows) {
	  typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	  boost::char_separator<char> sep(" \t,;");
	  Tokenizer tokens(sigdata, sep);
	  Tokenizer::iterator tokIter = tokens.begin();
	  if (tokIter!=tokens.end()) {
	    std::string chrName = *tokIter++;
	    if (tokIter!=tokens.end()) {
	      tokIter++;
	      if (tokIter!=tokens.end()) {
		tokIter++;
		uint32_t col = 0;
		for(;tokIter != tokens.end(); ++tokIter, ++col) {
		  if (col < cnbc[refIndex].cols) {
		    if ((*tokIter == "NaN") || (*tokIter == "NA")) cnbc[refIndex].sm[row][col] = lastVal[col]; // Carry forward old value
		    else {
		      cnbc[refIndex].sm[row][col] = boost::lexical_cast<TPrecision>(*tokIter) - 2;
		      lastVal[col] = cnbc[refIndex].sm[row][col];
		    }
		  }
		}
	      }
	    }
	  }
	  ++row;
	}
      }
      signalFile.close();
    }

    // Output matrix
    /*
    std::string outMatrix = c.outprefix + ".raw.gz";
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(outMatrix.c_str(), std::ios_base::out | std::ios_base::binary));
    for(uint32_t refIndex = 0; refIndex < cnbc.size(); ++refIndex) {
      for(uint32_t i = 0; i < cnbc[refIndex].rows; ++i) {
	dataOut << cnbc[refIndex].chr << '\t' << cnbc[refIndex].itv[i].istart << '\t' << cnbc[refIndex].itv[i].iend;
	for(uint32_t j = 0; j < cnbc[refIndex].cols; ++j) {
	  dataOut << '\t' << cnbc[refIndex].sm[i][j];
	}
	dataOut << std::endl;
      }
    }
    dataOut.pop();
    */
    
    // Iterate chromosomes
    for(uint32_t refIndex = 0; refIndex < cnbc.size(); ++refIndex) {
      // Segmentation
      Recap res;
      gflars(c, cnbc[refIndex].sm, res);
      
      // Debug
      //for(uint32_t i = 0; i < c.k; ++i ) std::cout << "RefIndex " << refIndex << ", Iteration " << i << ", " << res.lambda[i] << ", " << res.jump[i] << std::endl;
      
      // DP
      dpseg(c, cnbc[refIndex].sm, res);

      // Smooth signal
      SmoothSignal smoo;
      smoothsignal(cnbc[refIndex].sm, res.kbestjump, smoo);
      // updown stats
      updown(smoo);
      //expandpiecewiseconstant(smoo.jumps, smoo.smooth, res.smooth);
      
      // Output matrix
      std::string outFileName = c.outprefix + "." + cnbc[refIndex].chr + ".smooth.gz";
      boost::iostreams::filtering_ostream dataOutS;
      dataOutS.push(boost::iostreams::gzip_compressor());
      dataOutS.push(boost::iostreams::file_sink(outFileName.c_str(), std::ios_base::out | std::ios_base::binary));
      dataOutS << "chr\tstart\tend\tsignal\tsample\ttype" << std::endl;
      uint32_t istart = cnbc[refIndex].itv[0].istart;
      for(uint32_t i = 0; i<smoo.jumps.size(); ++i) {
	if (i) istart = cnbc[refIndex].itv[smoo.jumps[i-1]+1].istart;
	uint32_t iend = cnbc[refIndex].itv[smoo.jumps[i]].iend;
	for(uint32_t j = 0; j < cnbc[refIndex].cols; ++j) {
	  std::string type = "neutral";
	  if (smoo.smooth[i][j] > 0.5) type = "gain";
	  else if (smoo.smooth[i][j] < -0.5) type = "loss";
	  dataOutS << cnbc[refIndex].chr << '\t' << istart << '\t' << iend << '\t' << smoo.smooth[i][j] << "\tS" << j << "\t" << type << std::endl;
	}
      }
    
    
      // Debug dump
      std::cout << cnbc[refIndex].chr << std::endl;
      for(uint32_t r = 0; r<smoo.updown.shape()[0]; ++r) {
	for(uint32_t c = 0; c<smoo.updown.shape()[1]; ++c) {
	  std::cout << smoo.updown[r][c] << ',';
	}
	std::cout << std::endl;
      }
      /*
	for(uint32_t r = 0; r<smoo.smooth.shape()[0]; ++r) {
	for(uint32_t c = 0; c<smoo.smooth.shape()[1]; ++c) {
	std::cout << smoo.smooth[r][c] << ',';
	}
	std::cout << std::endl;
	}
      */
    }
  
    // End
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
    return 0;
  }


  int segment(int argc, char **argv) {
    SegmentConfig c;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("dpthreshold,d", boost::program_options::value<double>(&c.dpthreshold)->default_value(0.5), "DP threshold")
      ("epsilon,e", boost::program_options::value<double>(&c.epsilon)->default_value(1e-9), "epsilon error")
      ("kchange,k", boost::program_options::value<uint32_t>(&c.k)->default_value(100), "change points per chr")
      ("outprefix,p", boost::program_options::value<std::string>(&c.outprefix)->default_value("signal"), "output file prefix")
      ;

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.signal), "input signal matrix")
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
      std::cout << "Usage: " << argv[0] << " [OPTIONS] <signal.tsv>" << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    }
    
    // Parse matrix
    typedef std::vector<NormalizedBinCounts> TChrBinCounts;
    TChrBinCounts cnbc;
    if (!(boost::filesystem::exists(c.signal) && boost::filesystem::is_regular_file(c.signal) && boost::filesystem::file_size(c.signal))) {
      std::cerr << "Signal matrix is missing: " << c.signal.string() << std::endl;
      return 1;
    } else {
      // Get intervals & matrix dimensions for each chromosome
      std::ifstream signalFile(c.signal.string().c_str(), std::ifstream::in);
      if (signalFile.is_open()) {
	NormalizedBinCounts nbc;
	while (signalFile.good()) {
	  std::string sigdata;
	  getline(signalFile, sigdata);
	  typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	  boost::char_separator<char> sep(" \t,;");
	  Tokenizer tokens(sigdata, sep);
	  Tokenizer::iterator tokIter = tokens.begin();
	  if (tokIter!=tokens.end()) {
	    std::string chrName = *tokIter++;
	    if (tokIter!=tokens.end()) {
	      uint32_t start = boost::lexical_cast<uint32_t>(*tokIter++);
	      if (tokIter!=tokens.end()) {
		uint32_t end = boost::lexical_cast<uint32_t>(*tokIter++);
		int32_t col = 0;
		for(;tokIter != tokens.end(); ++tokIter) ++col;
		if (nbc.chr != chrName) {
		  if ((nbc.cols) && (nbc.rows)) cnbc.push_back(nbc);
		  nbc.chr = chrName;
		  nbc.cols = col;
		  //nbc.cols = 18;
		  nbc.rows = 0;
		  nbc.itv.clear();
		}
		nbc.itv.push_back(Interval(start, end));
		++nbc.rows;
	      }
	    }
	  }
	}
	if ((nbc.cols) && (nbc.rows)) cnbc.push_back(nbc);
	signalFile.close();
      }
    }
    if (cnbc.empty()) {
      std::cerr << "Signal matrix format is chr, start, end, signal, ..." << std::endl;
      return 1;
    }

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return runSegmentation(c, cnbc);
  }

}

#endif
