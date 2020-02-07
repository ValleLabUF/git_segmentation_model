#Analyze Snail Kite Data with Bayesian Partitioning Model

library(tidyverse)
library(progress)
library(furrr)
library(tictoc)
library(viridis)

source('gibbs functions.R')
source('helper functions.R')
source('gibbs sampler.R')


###############################
#### Load and Prepare Data ####
###############################

dat<- read.csv("Snail Kite Gridded Data_TOHO.csv", header = T, sep = ",")

#remove IDs w < 3 occupied grid cells
dat.ex<- dat %>% group_by(id) %>% filter(length(unique(grid.cell)) < 3) %>% ungroup()
dat.ex.list<- df.to.list(dat = dat.ex)
dat.ex.list<- map(dat.ex.list, function(x) x %>% mutate(time1 = 1:length(x)))  #add row for obs #
dat.ex<- map_dfr(dat.ex.list, `[`) %>% mutate(tseg = 1)

dat<- dat %>% group_by(id) %>% filter(length(unique(grid.cell)) >= 3) %>% ungroup()
dat.list<- df.to.list(dat = dat)
dat.list<- map(dat.list, function(x) x %>% mutate(time1 = 1:length(x)))  #add row for obs #

#only select necessary cols and re-number grid cell IDs
dat.long<- map_dfr(dat.list, `[`) %>% dplyr::select(id, grid.cell, time1)  #create DF
names(dat.long)[2]<- "loc.id"
dat.long$loc.id<- dat.long$loc.id %>% factor()
levels(dat.long$loc.id)<- 1:length(unique(dat.long$loc.id))  #change from raw to modified cell ID
dat.long$loc.id<- dat.long$loc.id %>% as.character() %>% as.numeric()

#convert back to list
dat.list2<- df.to.list(dat.long)



#######################################
#### Run Gibbs Sampler for all IDs ####
#######################################

ngibbs = 10000

#priors
alpha=0.01

## Run Gibbs sampler
plan(multisession)  #run all MCMC chains in parallel
                    #select "multiprocess" if Unix or macOS & "multisession" if Windows
                    #refer to future::plan() for more details

dat.res<- space_segment(data = dat.list2, ngibbs = ngibbs, alpha = alpha)
###Takes 8 min to run for 10000 iterations for all IDs


## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
identity<- names(dat.list2)

traceplot(data = dat.res$nbrks, type = "nbrks", identity = identity)
traceplot(data = dat.res$LML, type = "LML", identity = identity)


## Determine maximum likelihood (ML) for selecting breakpoints
ML<- apply(dat.res$LML, 1, function(x) getML(dat = x, nburn = 500))
brkpts<- getBreakpts(dat = dat.res$brkpts, ML = ML, identity = identity)  


## Heatmaps
plot.heatmap(data = dat.list2, brkpts = brkpts, dat.res = dat.res, type = "loc")


######################################
#### Assign Spatial Time Segments ####
######################################



dat_out<- map(dat.list, assign.time.seg, brkpts = brkpts) %>% map_dfr(`[`)  #assign time seg and make as DF
dat_out<- rbind(dat_out, dat.ex)  #bring back in excluded data occupying < 3 cells
dat_out<- dat_out[order(dat_out$id, dat_out$date),]  #reorder DF by id and date

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/activcenter_subset_locations")
write.csv(dat_out, "Snail Kite Gridded Data_TOHO.csv", row.names = F)

