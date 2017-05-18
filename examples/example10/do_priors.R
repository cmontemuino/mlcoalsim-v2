#PRIORS:
niter <- 1e2
k <- 0

ne.bott <- runif(niter,-4,0)
ne.bott <- 10^(ne.bott)
k <- k + 1
write.table(file=sprintf("prior%03d.txt",k), ne.bott,quote=F,row.names=F)

time.to.bott <- runif(niter,0,0.10)
k <- k + 1
write.table(file=sprintf("prior%03d.txt",k),time.to.bott,quote=F,row.names=F)

ne.ancestral <- runif(niter,-1,3)
ne.ancestral <- 10^(ne.ancestral)
k <- k + 1
write.table(file=sprintf("prior%03d.txt",k),ne.ancestral,quote=F,row.names=F)

time.to.outg <- runif(niter,2,12)
k <- k + 1
write.table(file=sprintf("prior%03d.txt",k),time.to.outg,quote=F,row.names=F)

sexratio <- runif(niter,-1,1)
sexratio <- 10^(sexratio)
k <- k + 1
write.table(file=sprintf("prior%03d.txt",k), sexratio,quote=F,row.names=F)

#AUTOSOMES
theta <- runif(niter,0.001,0.075)
for(i in 1:5) {
	k <- k + 1
	write.table(file=sprintf("prior%03d.txt",k),theta,quote=F,row.names=F)
}
#CHROM X
factor.theta <- runif(niter,log10(0.1),log10(2))
factor.theta <- 10^(factor.theta)
for(i in 1:5) {
	k <- k + 1
	write.table(file=sprintf("prior%03d.txt",k),theta*factor.theta,quote=F,row.names=F)
}

recombination <- runif(niter,0,0.075)
k <- k + 1
write.table(file=sprintf("prior%03d.txt",k), recombination,quote=F,row.names=F)
