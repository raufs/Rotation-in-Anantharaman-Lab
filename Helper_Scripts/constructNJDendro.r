library(ape)
library(phytools)
args = commandArgs(trailingOnly=TRUE)

dat <- read.table(args[1], header=T, sep='\t', row.names=1)
d <- dist(dat)
njt <- nj(d)
mid.njt <- midpoint.root(njt)
write.tree(mid.njt, file=args[2])