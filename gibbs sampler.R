rm(list=ls(all=TRUE))
set.seed(3)

setwd('U:\\GIT_models\\git_segmentation_model')
source('gibbs functions.R')
dat=read.csv('fake data.csv',as.is=T)

#priors
alpha=0.01

#useful stuff
max.time=max(dat$time1)
nloc=max(dat$loc.id)

#starting values
breakpt=floor(mean(dat$time1))

ngibbs=10000
for (i in 1:ngibbs){
  print(i)
  breakpt=samp.move(breakpt=breakpt,max.time=max.time,dat=dat,
                    alpha=alpha,nloc=nloc)   
}
length(breakpt)
abline(v=breakpt,lty=3,col='green')