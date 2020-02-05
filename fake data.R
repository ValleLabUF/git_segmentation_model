rm(list=ls(all=TRUE))
library('MCMCpack')
set.seed(2)

nobs=9000
nseg=30
tmp=runif(nseg)
prob=tmp/sum(tmp); prob
partition=rmultinom(1,size=nobs,prob=prob)
seg.index=rep(1:nseg,times=partition)

nloc=300
prob=rdirichlet(nseg,rep(0.01,nloc))

obs=rep(NA,nobs)
for (i in 1:nobs){
  tmp=rmultinom(1,size=1,prob=prob[seg.index[i],])
  obs[i]=which(tmp==1)
}
plot(obs)
tmp=cumsum(partition)
abline(v=tmp[-length(tmp)])

obs1=data.frame(loc.id=obs)
obs1$time1=1:nobs

setwd('U:\\GIT_models\\git_segmentation_model')
write.csv(obs1,'fake data.csv',row.names=F)