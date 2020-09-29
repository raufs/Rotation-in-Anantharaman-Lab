library(ape)
library(phytools)
require(data.table)

args = commandArgs(trailingOnly=TRUE)

dat <- read.table(args[1], header=F, sep='\t', row.names=1)
colnames(dat) <- rownames(dat)
d <- as.dist(as.matrix(dat))
njt <- nj(d)
mid.njt <- midpoint.root(njt)
write.tree(mid.njt, file=args[2])