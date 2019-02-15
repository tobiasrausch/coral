library(ggplot2)
library(scales)

chrs = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20", "chr21", "chr22", "chrX")
#chrs = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X")

args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=T)
x = x[x$chr %in% chrs,]
x$chr = factor(x$chr, levels=chrs)

for(chrname in unique(x$chr)) {
 print(chrname)
 sub = x[x$chr == chrname,]
 sub$mid = (sub$end + sub$start) / 2
 p = ggplot(data=sub, aes(x=mid, y=sub[,5]))
 p = p + geom_point(pch=21, size=0.5)
 p = p + ylab("Copy-number") + xlab(chrname)
 p = p + scale_x_continuous(labels=comma)
 p = p + ylim(0,8)
 p = p + theme(axis.text.x = element_text(angle=45, hjust=1))
 ggsave(paste0("plot.", chrname, ".png"), width=24, height=6)
 print(warnings())
}
 
