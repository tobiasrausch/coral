library(ggplot2)
library(scales)
library(gtable)
library(grid)

chrNamesLong = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20", "chr21", "chr22", "chrX")
chrNamesShort = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X")

args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=T)
if (sum(x$chr %in% chrNamesLong) > sum(x$chr %in% chrNamesShort)) { chrs = chrNamesLong; } else { chrs = chrNamesShort; }
x = x[x$chr %in% chrs,]
x$chr = factor(x$chr, levels=chrs)
seg = read.table(args[2], header=T)
seg = seg[seg$chr %in% chrs,]
seg$chr = factor(seg$chr, levels=chrs)
plotBaf = T
if (mean(is.na(x[,7])) > 0.5) { plotBaf = F; }

p = ggplot(data=x, aes(x=start, y=x[,6]))
p = p + geom_point(pch=21, size=0.5)
p = p + geom_segment(data = seg, aes(x=start, xend=end, y=cn, yend=cn), color="lightblue")
p = p + xlab("Chromosome")
p = p + ylab("Copy-number")
p = p + scale_x_continuous(labels=comma)
p = p + facet_grid(. ~ chr, scales="free_x", space="free_x")
p = p + ylim(0,8)
p = p + theme(axis.text.x = element_text(angle=45, hjust=1))
if (plotBaf) {
 q = ggplot(data=x, aes(x=start, y=x[,7]))
 q = q + geom_point(pch=21, size=0.5)
 q = q + ylab("Obs / Exp MAF of het. SNPs") + xlab("Chromosome")
 q = q + geom_segment(data = seg, aes(x=start, xend=end, y=maf, yend=maf), color="lightblue")
 q = q + scale_x_continuous(labels=comma)
 q = q + facet_grid(. ~ chr, scales="free_x", space="free_x")
 q = q + ylim(0,1.5)
 q = q + theme(axis.text.x = element_text(angle=45, hjust=1))
 g1 = ggplotGrob(p)
 g2 = ggplotGrob(q)
 g = rbind(g1, g2, size="first")
 g$widths = unit.pmax(g1$widths, g2$widths)
 ggsave(g, file="plot.wholegenome.png", width=24, height=6)
} else {
 ggsave(p, file="plot.wholegenome.png", width=24, height=6)
}
print(warnings())

for(chrname in unique(x$chr)) {
 print(chrname)
 sub = x[x$chr == chrname,]
 local = seg[seg$chr == chrname,]
 p = ggplot(data=sub, aes(x=start, y=sub[,6]))
 p = p + geom_point(pch=21, size=0.5)
 p = p + geom_segment(data = local, aes(x=start, xend=end, y=cn, yend=cn), color="lightblue")
 p = p + ylab("Copy-number") + xlab(chrname)
 p = p + scale_x_continuous(labels=comma, breaks = scales::pretty_breaks(n=20))
 p = p + ylim(0,8)
 p = p + theme(axis.text.x = element_text(angle=45, hjust=1))
 q = ggplot(data=sub, aes(x=start, y=sub[,7]))
 #print(sd(sub[,7], na.rm=T))
 #print(mean(sub[,7], na.rm=T))
 #print(mean(is.na(sub[,7])))
 if (plotBaf) {
  q = q + geom_point(pch=21, size=0.5)
  q = q + ylab("Obs / Exp MAF of het. SNPs") + xlab(chrname)
  q = q + geom_segment(data = local, aes(x=start, xend=end, y=maf, yend=maf), color="lightblue")
  q = q + scale_x_continuous(labels=comma, breaks = scales::pretty_breaks(n=20))
  q = q + theme(axis.text.x = element_text(angle=45, hjust=1))
  q = q + ylim(0,1.5)
  g1 = ggplotGrob(p)
  g2 = ggplotGrob(q)
  g = rbind(g1, g2, size="first")
  g$widths = unit.pmax(g1$widths, g2$widths)
  ggsave(g, file=paste0("plot.", chrname, ".png"), width=24, height=6)
 } else {
  ggsave(p, file=paste0("plot.", chrname, ".png"), width=24, height=6)
 }
 print(warnings())
}

