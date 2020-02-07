gibbs.time.seg=function(data, ngibbs) {
  set.seed(1)
  
  uni.id=unique(data$id)
  
  #priors
  alpha=0.01
  
  #to store results
  res.brks=vector("list", ngibbs)
  res.LML=matrix(NA,1,(ngibbs+1))
  res.nbrks=matrix(NA,1,(ngibbs+1))
  store.param=matrix(NA,ngibbs,2)
    
  #re-define loc.id based only on those visited by this individual
  uni.loc=unique(data$loc.id)
  aux=data.frame(loc.id=uni.loc,loc.id1=1:length(uni.loc))
  dat1=merge(data,aux,all=T)
  dat1$loc.id=dat1$loc.id1
  dat=dat1[order(dat1$time1),c('loc.id','time1')]
    
  #useful stuff
  max.time=max(dat$time1)
  nloc=max(dat$loc.id)
  
  #starting values
  breakpt=floor(mean(dat$time1))
  
  for (i in 1:ngibbs){
    vals=samp.move(breakpt=breakpt,max.time=max.time,dat=dat,
                   alpha=alpha,nloc=nloc)  
    breakpt=vals[[1]]
    
    #store results
    res.brks[[i]]<- breakpt
    store.param[i,]=c(length(vals[[1]]), vals[[2]])  # nbrks and LML
  }
  
  tmp=store.param[,1]
  res.nbrks[1,]=c(uni.id,tmp)
  colnames(res.nbrks)<- c('id', paste0("Iter_",1:ngibbs))
  
  tmp=store.param[,2]
  res.LML[1,]=c(uni.id,tmp)
  colnames(res.LML)<- c('id', paste0("Iter_",1:ngibbs))
    
  list(breakpt=res.brks, nbrks=res.nbrks, LML=res.LML)
}



#----------------------------------------------------
space_segment=function(data, ngibbs) {
  
  tic()  #start timer
  mod<- future_map(data, function(x) gibbs.time.seg(data = x, ngibbs = ngibbs), .progress = TRUE)
  toc()  #provide elapsed time
  
  
  brkpts<- map(mod, 1)  #create list of all sets breakpoints by ID
  
  nbrks<- map_dfr(mod, 2) %>% t() %>% data.frame()  #create DF of number of breakpoints by ID
  names(nbrks)<- c('id', paste0("Iter_",1:ngibbs))
  
  LML<- map_dfr(mod, 3) %>% t() %>% data.frame()  #create DF of LML by ID
  names(LML)<- c('id', paste0("Iter_",1:ngibbs))
  
  
  list(brkpts = brkpts, nbrks = nbrks, LML = LML)
}


