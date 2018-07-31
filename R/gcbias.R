library(ggplot2)
library(reshape2)
library(scales)

args=commandArgs(trailingOnly=TRUE)
x= read.table(args[1], header=F)
ttl = args[2]
colnames(x) = c("gccount", "count", "type")
gcmax = nrow(x[x$type=="Reference",]) - 1
x$gc = x$gccount / gcmax

png(paste0(ttl, ".count.png"), width=800, height=400)
ggplot(data=x, aes(x=gc, y=count)) + geom_freqpoly(aes(color=type), stat="identity") + facet_wrap(~type, nrow=2, scales="free") + scale_y_continuous(labels=comma) + xlab("GC-content") + ylab("Count") + ggtitle(ttl)
dev.off()

z=cbind(x[x$type=="Sample",], refcount = x[x$type=="Reference",]$count)
z=z[z$gc>=0.25 & z$gc<=0.75,]
sc = sum(z$count) / sum(as.numeric(z$refcount))
z$div = z$count / (sc * z$refcount)
png(paste0(ttl, ".bias.png"), width=800, height=400)
ggplot(data=z, aes(x=gc, y=div)) + geom_line() + xlab("GC-content") + ylab("Observed/Expected") + ggtitle(ttl)
dev.off()
