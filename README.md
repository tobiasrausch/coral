[![Build Status](https://travis-ci.org/tobiasrausch/coral.svg?branch=master)](https://travis-ci.org/tobiasrausch/coral)
[![GitHub license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://raw.githubusercontent.com/tobiasrausch/coral/master/LICENSE)
[![GitHub Issues](https://img.shields.io/github/issues/tobiasrausch/coral.svg)](https://github.com/tobiasrausch/coral/issues)


Installation Coral
------------------

`git clone --recursive https://github.com/tobiasrausch/coral.git`

`cd coral/`

`make all`


Reference Bundle
----------------

You then need to download the reference bundle with mappability maps and 1000 Genomes SNPs (~737Mb).

[Download Reference Bundle](https://drive.google.com/uc?export=download&id=1EfM4SdIYv4vAwz-Ri9nMxoCd4vZMMBf4)

`tar -xzf referenceBundle.tar.gz`


Running Coral
-------------

`./src/coral call -v referenceBundle/GRCh37.snps.bcf -g GRCh37.fa -m referenceBundle/Homo_sapiens.GRCh37.dna.primary_assembly.fa.r101.s501.blacklist.gz -s <SampleName> -o <outprefix> -l <control.bam> <tumor.bam>`

Segmentation plots you can generate using:

`Rscript ./R/rd.R outprefix.adaptive.cov.gz outprefix.segment.gz`


Germline/Somatic Copy-Number Calling
------------------------------------

Work-in-progress



