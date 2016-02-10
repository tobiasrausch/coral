Installation
------------

`git clone --recursive https://github.com/tobiasrausch/strandseq.git`

`cd strandseq/`

`make all`


Pre-processing Strand-Seq data
------------------------------

The pre-processing method categorizes Watson-Watson, Crick-Crick and Watson-Crick chromosomes and merges these into a ww.bam and wc.bam file. Using common SNPs (csnp.hg38.bcf) the pre-processing method also aligns Watson-Crick in phase.

`./src/preprocess -v ./data/csnp.hg38.bcf -s HG00512.ww.bam -d HG00512.wc.bam single_cell_HG00512_bam/*.bam`


