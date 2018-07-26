library(ggplot2)
library(scales)
library(grid)
library(RColorBrewer)
library(reshape2)

args=commandArgs(trailingOnly=TRUE)
x= read.table(args[1], header=T)

pdf("cells.pdf", width=12, height=6)
j = 0
for(s in unique(x$sample)) {
 print(s)
 sub = x[x$sample==s, 2:5]
 sub = melt(sub, id.vars=c("pos", "chr"))
 p = ggplot(data=sub, aes(x=pos, y=value))
 p = p + ggtitle(s)
 p = p + geom_line(aes(group=variable, color=variable))
 p = p + facet_wrap(~chr, scales="free")
 print(p)
 j = j + 1
 if (j > 5) { break; }
}
dev.off()
print(warnings())
