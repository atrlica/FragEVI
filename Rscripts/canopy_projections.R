# library(raster)
library(data.table)



## all of the simulator dump files
# load(file=paste("processed/boston/biom_street/ann.npp.street.v", vers, ".weighted.", y, ".sav", sep=""))
# load(file=paste("processed/boston/biom_street/num.trees.street.v", vers, ".weighted.", y, ".sav", sep=""))
# load(cage.dbh, file=paste("processed/boston/biom_street/dbh.street.v", vers, ".weighted.", y, ".sav", sep=""))
# load(cage.genus, file=paste("processed/boston/biom_street/genus.street.v", vers, ".weighted.", y, ".sav", sep=""))
# load(biom.track, file=paste("processed/boston/biom_street/biom.track.street.v", vers, ".weighted.", y, ".sav", sep=""))
# load(index.track, file=paste("processed/boston/biom_street/index.track.street.v", vers, ".weighted.", y, ".sav", sep=""))
# load(proc.track, file=paste("processed/boston/biom_street/proc.track.street.v", vers, ".weighted.", y, ".sav", sep=""))
# load(cage.biom.sim, file=paste("processed/boston/biom_street/biom.sim.street.v", vers, ".weighted.", y, ".sav", sep=""))
# load(attempts.track, file=paste("processed/boston/biom_street/attempts.track.street.v", vers, ".weighted.", y, ".sav", sep=""))
# load(cage.wts, file=paste("processed/boston/biom_street/wts.street.v", vers, ".weighted.", y, ".sav", sep=""))
# 

### what does the data look like?
# load("processed/boston/biom_street/dbh.street.v4.weighted.10000.sav")
# load("processed/boston/biom_street/ann.npp.street.v4.weighted.10000.sav")
# median(cage.ann.npp[[1]])
# quantile(cage.ann.npp[[1]], probs=c(0.45, 0.55))

# #### a question arises: Is the population of trees to select the one that gets closest to the pixel biomass, or the one nearest to the median NPP estimate?
# load("processed/boston/biom_street/biom.sim.street.v4.weighted.10000.sav")
# quantile(cage.biom.sim[[1]], probs=c(0.45, 0.55))
# j <- which(cage.biom.sim[[1]]>2638 & cage.biom.sim[[1]]<2650)
# load("processed/boston/biom_street/biom.track.street.v4.weighted.10000.sav")
# biom.track[1]

### mortality rate equation
mort <- function(x){(0.0008133*(x^2))-(0.0642407*x)+4.0614503}
# x <- seq(1,100)
# plot(x, mort(x))

## growth rate equation
b0 <- -2.48
b1 <- 2.4835 ## these are eastern hardwood defaults
biom.pred <- function(x){exp(b0+(b1*log(x)))}
# street <- read.csv("docs/ian/Boston_Street_Trees.csv") ## Street tree resurvey, biomass 2006/2014, species
# street <- as.data.table(street)
# street[, record.good:=0] ## ID records with good data quality
# street[is.finite(dbh.2014) & is.finite(dbh.2006) & dead.by.2014==0 & Health!="Poor", record.good:=1] #2603 records good
# street[dbh.2006<5, record.good:=0] ## filter out the handfull of truly tiny trees
# street[record.good==1, biom.2014:=biom.pred(dbh.2014)]
# street[record.good==1, biom.2006:=biom.pred(dbh.2006)]
# street[record.good==1, growth.ann:=(biom.2014-biom.2006)/8]
# street[record.good==1, growth.ann.rel:=growth.ann/biom.2006]
# street[record.good==1, dbh.ann:=(dbh.2014-dbh.2006)/8]
# street[record.good==1, dbh.ann.rel:=dbh.ann/dbh.2006]
# plot(street[record.good==1, dbh.2006], street[record.good==1, dbh.ann.rel])
# ### direct dbh growth equation from street tree records
# mod.street.dbh.nls <- nls(dbh.ann.rel ~ exp(a + b * log(dbh.2006)), data=street[record.good==1,], start=list(a=0, b=0)) 
# mm <- summary(mod.street.dbh.nls) ## RSE 0.0404
# street[record.good==1, mean(dbh.ann.rel)] ## mean 0.04
# street[record.good==1, sd(dbh.ann.rel)/sqrt(length(dbh.ann.rel))] ## SEM 0.001, fuck
## mm coefficients: had to hard code bc street data not avail on every machine
mm.a=0.2281686
mm.b=(-1.1043188)

dbhg.pred <- function(x){exp(mm.a+(mm.b*log(x)))} ## function for predicting next year's relative dbh growth from this year's dbh

### set some scenario options
# scenario = "BAU"
scenario = "oldies"

## loops for simulating growth+mortality in pixel dbh samples
## pull in all the simulator samples
vers <- 4 ## what simulator results version are we dealing with?
setwd("/projectnb/buultra/atrlica/FragEVI")
obj.list <- list.files("processed/boston/biom_street/")
obj.list <- obj.list[grep(obj.list, pattern=paste("dbh.street.v", vers, sep=""))]
obj.list <- sub('.*weighted\\.', '', obj.list)
obj.list <- (sub('\\..*', '', obj.list))

library(stringr)
### avoid processing shit that's already done
already <- list.files("processed/boston/biom_street")
already <- already[grep(already, pattern=".resim.")]
already <- already[grep(already, pattern=paste(".", scenario, ".", sep=""))]
# sub('.*\\.[[:digit:]]', '', already)
# sub('.*\\.[0-9]+', '', already)
# a <- str_extract(already, pattern='\\.[0-9]+')
# sub('.*\\.scenario.*', '', already)
# already <- (sub('\\.V1.*', '', already))

hitlist <- as.character(obj.list[!(obj.list %in%  already)])
hitlist[hitlist=="100000"] <- "1e+05"

# 
# if(scenario=="BAU"){ ## standard growth and mortality assumptions
#   mort <- function(x){(0.0008133*(x^2))-(0.0642407*x)+4.0614503}
#   b0 <- -2.48
#   b1 <- 2.4835 ## these are eastern hardwood defaults
#   biom.pred <- function(x){exp(b0+(b1*log(x)))}  
# }
# 
# if(scenario=="oldies"){ ## reduce mortality 50% in >40cm dbh (largest 25%)
#   ## mortality
#   mort <- function(x){
#     if(x<40){
#       (0.0008133*(x^2))-(0.0642407*x)+4.0614503
#       return(x)
#     } else{
#       ((0.0008133*(x^2))-(0.0642407*x)+4.0614503)*0.5
#       return(x)
#     }
#   }
#   ## growth
#   b0 <- -2.48
#   b1 <- 2.4835 ## these are eastern hardwood defaults
#   biom.pred <- function(x){exp(b0+(b1*log(x)))}  
#   
# }
# d <- 1:100
# plot(d, mort(d))
# mort(10)
# mort(50)
for(o in hitlist){
  ### initialize read-in files and save dumps
  load(paste("processed/boston/biom_street/dbh.street.v", vers, ".weighted.", o, ".sav", sep="")) # cage.dbh
  load(paste("processed/boston/biom_street/biom.sim.street.v", vers, ".weighted.", o, ".sav", sep="")) # cage.biom.sim
  load(paste("processed/boston/biom_street/ann.npp.street.v", vers, ".weighted.", o, ".sav", sep="")) # cage.ann.npp
  deaths.sav <- list()
  dbh.sav <- cage.dbh
  
  for(pix in 1:length(dbh.sav)){ ## the number of pixels in this sim chunk
    if(length(cage.biom.sim[[pix]])>=40){ ## only process if enough simulations successfully completed
      deaths.track <- numeric()
      for(a in 1:100){   ### resimulate dbh history in every pixel 100 times
        
        ## can select dbh populations based on proximity to median simulated biomass...
        # biom.lims <- quantile(cage.biom.sim[[pix]], probs=c(0.40, 0.6)) ## figure out which of the simulations to draw and modify
        # j <- sample(which(cage.biom.sim[[pix]]>=biom.lims[1] & cage.biom.sim[[pix]]<=biom.lims[2]),1) ## get a random dbh sample from the the middle 10% of simulator results close to the target biomass
        
        ## or can select dbh populations based on how close they are to median npp (not 100% overlapping)
        npp.lims <- quantile(cage.ann.npp[[pix]], probs=c(0.25, 0.75)) ## restrict which of the simulations to draw and modify
        j <- which(cage.ann.npp[[pix]]>=npp.lims[1] & cage.ann.npp[[pix]]<=npp.lims[2])
        if(length(j)>=4){
          j <- sample(j,1) ## get a random dbh sample from the the middle 10% of simulator results close to the target biomass
          }else{j <- sample(1:length(cage.biom.sim[[pix]]), 1)} ## or just sample the whole simulation collection if all simulations are nearly identical
        
        #### load up a dbh sample and resimulate each tree for 36 consecutive years
        tree.samp <- unlist(cage.dbh[[pix]][j]) ## what initial trees are present in this simulator result
        dbh.sav[[pix]][[a]] <- numeric() ## clear contents to store the new simulated results
        dbh.track <- numeric() ## keep track of the final size of each tree post growth and mortality
        deaths <- 0 ## keep track of how many times trees die in this resim
        for(y in 1:length(tree.samp)){ ## simulating each tree history separately
          d <- tree.samp[y] ## grab dbh record and go to town
          # track <- d
          for(e in 1:36){ ## test each tree for 36 years (2006-2040)
            deathwatch <- mort(d) ## recalculate death probability given dbh as it changes
            if(scenario=="oldies" & d>40){deathwatch <- deathwatch*0.5} ## in "oldies" scenario, cut mortality on large trees by 50%
            if(rbinom(1, 1, deathwatch/100)){ ## if tree doesn't make it
              # print("he's dead Jim")
              d <- 5 ## immediately replant a 5cm tree if it dies
              deaths <- deaths+1
            } else{
              # print("It's alive!!")
              d <- d*(1+dbhg.pred(d)) ## if it survives update size via growth eq.
            }
            # track <- c(track, d)
            # ptrack <- c(ptrack, deathwatch)
          } ## loop for 36 year projector
          dbh.track <- c(dbh.track, d) ## tally of the new dbhs after growth and morts applied
          # plot(track)
          # track <- numeric()
          # print(paste("deaths =", deaths))
          # print(paste("2040 tree is cm", d))
          # print(paste("2006 tree was cm", tree.samp[y]))
        } # loop for all trees in this dbh sample
        deaths.track <- c(deaths.track, deaths) ## tally of total deaths in each pixel simulation
        dbh.sav[[pix]][[a]] <- dbh.track ## record the final dbh arrangement (a vector of variable length) of each pixel simulation
      } # loop for number of resims in this pixel
      deaths.sav[[pix]] <- deaths.track ##  vector n=100 with number of deaths recorded in each pixel resimulation (i.e. how much replanting is needed in every pixel over 40 years)
      if(pix%%500 == 0){print(paste("chunk", o, "resimmed pix", pix))} ## give status updates
    } else{  ## if too few successful simulations for this pixel
      dbh.sav[[pix]][[a]] <- NA
      deaths.sav[[pix]] <- NA
      }
  } # loop for this pixel
 save(dbh.sav, file = paste("processed/boston/biom_street/dbh.resim.scenario.", scenario, ".", o, ".V1.sav", sep=""))
 save(deaths.track, file=paste("processed/boston/biom_street/death.track.scenario.", scenario, ".", o, ".V1.sav", sep=""))
}
print("finished resim run scenario 0")
# 
# load("processed/boston/biom_street/death.track.scenario-0.BAU10000.V1.sav")
# a <- unlist(death.track)
# hist(unlist(lapply(lapply(dbh.sav[[500]], biom.pred), sum)))
# hist(unlist(lapply(lapply(cage.dbh[[500]], biom.pred), sum)))
# unlist(death.track)

## how many trees do you get in each re-sim  
# num.resim.trees <- unlist(lapply(dbh.sav[[pix]], length))
# num.sim.trees <- unlist(lapply(cage.dbh[[pix]], length))
# hist(num.sim.trees)
# hist(num.resim.trees)
# 
# plot(unlist(lapply(cage.dbh[[400]], median)), unlist(lapply(dbh.sav[[400]], median)))
# abline(a=0, b=1)
