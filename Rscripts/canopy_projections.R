# library(raster)
library(data.table)
setwd("/projectnb/buultra/atrlica/FragEVI")



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

### dbh growth equation (from street tree model)
mm.a=0.2281686
mm.b=(-1.1043188)
dbhg.pred <- function(x){exp(mm.a+(mm.b*log(x)))} ## function for predicting next year's relative dbh growth from this year's dbh
## mm coefficients: had to hard code bc street data not avail on every machine

######
### set some scenario options
## available options: BAU, oldies, lowmort, highmort, lowreplant, 
scenario <- c("BAU", "highmort", "lowreplant", "oldies", "lowmort", "slowreplant")
# scen.l <- c("BAU", "highmort", "lowreplant", "oldies", "lowmort")
# scen.num <- which(scen.l==scenario) ## apply a scenario number
## record specific parameter sets
resim.vers <- 1 ## what are we labeling this round of resims?
vers <- 4 ## what simulator results version are we dealing with?
### BAU/default factors
npp.quant.range <- c(0.25, 0.75) ## what dbh samples to draw from
mort.mod <- 1 ## modification to mortality rate (scenario specific)
default.sizecutoff <- 10000 ## nothing is this big
largemort.mod <- 1
replant.factor <- 1 ## rate of replanting
delay.factor <- c(0,0) ## low/hi on years of delay from death to replanting
dbg.big <- default.sizecutoff

## set options for scenario processing
highmort.mortfactor <- 1.25
lowmort.mortfactor <- 0.5
lowreplant.factor <- 0.5 ## what fraction of morts are replanted in "lowreplant"
oldies.mortfactor <- 0.5 ## how much to reduce mortalities in large trees
oldies.sizecutoff <- 40 ## how to define "big" trees
slowreplant.delay <- c(0,2) ## 1 to 3 year delay between death and regrowth starting


# record run parameters to text file
params.list <- list(c(scenario, resim.vers, vers,
                    highmort.mortfactor, lowmort.mortfactor, 
                    npp.quant.range, lowreplant.factor, oldies.sizecutoff, 
                    oldies.mortfactor))
sink(paste0("processed/boston/biom_street/resim_params_list_V", resim.vers, ".txt"))
cat(c("Resimulation version ", resim.vers), "\n")
cat(c("Simulator results version ", vers), "\n")
cat(c("scenarios = ", scenario), "\n")
cat(c("NPP selection quantiles = ", npp.quant.range), "\n")
cat(c("highmort mortality factor = ", highmort.mortfactor), "\n")
cat(c("lowmort mortality factor = ", lowmort.mortfactor), "\n")
cat(c("oldies mortality factor = ", oldies.mortfactor), "\n")
cat(c("oldies size cutoff = ", oldies.sizecutoff), "\n")
cat(c("lowreplant replanting rate = ", lowreplant.factor), "\n")
sink()

## loops for simulating growth+mortality in pixel dbh samples
## pull in all the simulator samples of tree dbh in every pixel
obj.list <- list.files("processed/boston/biom_street/")
obj.list <- obj.list[grep(obj.list, pattern=paste("dbh.street.v", vers, sep=""))]
obj.list <- sub('.*weighted\\.', '', obj.list)
obj.list <- (sub('\\..*', '', obj.list))

## loop the simulated pixel chunks
for(s in 1:length(scenario)){
  ### modify simulator factors according to scenario construction
  if(scenario=="highmort"){
    mort.mod <- highmort.mortfactor
    dbh.big <- default.sizecutoff 
    largemort.mod <- 1
    replant.factor <- 1 ## rate of replanting
    delay.factor <- c(0,0) ## low/hi on years of delay from death to replanting
  }
  if(scenario=="lowmort"){
    mort.mod <- lowmort.mortfactor
    dbh.big <- default.sizecutoff 
    largemort.mod <- 1
    replant.factor <- 1 ## rate of replanting
    delay.factor <- c(0,0) ## low/hi on years of delay from death to replanting
  }
  if(scenario=="oldies"){
    mort.mod <- 1
    dbh.big <- oldies.sizecutoff 
    largemort.mod <- oldies.mortfactor
    replant.factor <- 1 ## rate of replanting
    delay.factor <- c(0,0) ## low/hi on years of delay from death to replanting
  }
  if(scenario=="lowreplant"){
    replant.factor=0.5
    mort.mod <- 1
    dbh.big <- default.sizecutoff 
    largemort.mod <- 1
    delay.factor <- c(0,0) ## low/hi on years of delay from death to replanting
    
  }
  if(scenario=="slowreplant"){
    delay.factor <- slowreplant.delay
    mort.mod <- 1
    dbh.big <- default.sizecutoff 
    largemort.mod <- 1
    replant.factor <- 1 ## rate of replanting
  }
  print(paste("starting resim run scenario", scenario[s]))
  ### avoid processing shit that's already done
  already <- list.files("processed/boston/biom_street")
  already <- already[grep(already, pattern=".resim.")]
  already <- already[grep(already, pattern=paste(".", scenario[s], ".", sep=""))]
  # sub('.*oldies\\.', '', already)
  already <- strsplit(already, split = "[.]")
  already <- sapply(already, "[[" ,5)
  
  hitlist <- as.character(obj.list[!(obj.list %in%  already)])
  hitlist[hitlist=="100000"] <- "1e+05"
  
  for(o in hitlist){
    ### initialize read-in files and save dumps
    load(paste("processed/boston/biom_street/dbh.street.v", vers, ".weighted.", o, ".sav", sep="")) # cage.dbh
    load(paste("processed/boston/biom_street/biom.sim.street.v", vers, ".weighted.", o, ".sav", sep="")) # cage.biom.sim
    load(paste("processed/boston/biom_street/ann.npp.street.v", vers, ".weighted.", o, ".sav", sep="")) # cage.ann.npp
    deaths.sav <- list()
    dbh.sav <- cage.dbh
    
    for(pix in 1:length(dbh.sav)){ ## the number of pixels in this sim chunk
      if(length(cage.biom.sim[[pix]])>=40){ ## only process if enough simulations successfully completed in this pixel
        deaths.track <- numeric()
        for(a in 1:100){   ### resimulate every pixel X number of times
          
          ## can select dbh populations based on proximity to median simulated biomass...
          # biom.lims <- quantile(cage.biom.sim[[pix]], probs=c(0.40, 0.6)) ## figure out which of the simulations to draw and modify
          # j <- sample(which(cage.biom.sim[[pix]]>=biom.lims[1] & cage.biom.sim[[pix]]<=biom.lims[2]),1) ## get a random dbh sample from the the middle 10% of simulator results close to the target biomass
          
          ## or can select dbh populations based on how close they are to median npp (not 100% overlapping)
          npp.lims <- quantile(cage.ann.npp[[pix]], probs=npp.quant.range) ## restrict which of the simulations to draw and modify
          j <- which(cage.ann.npp[[pix]]>=npp.lims[1] & cage.ann.npp[[pix]]<=npp.lims[2])
          if(length(j)>=4){ ## if there's very enough dbh samples that meet your criteria...
            j <- sample(j,1) ## get a random dbh sample selected range of simulator results close to the target biomass
          }else{j <- sample(1:length(cage.biom.sim[[pix]]), 1)} ## or just sample the whole simulation collection if all simulations are nearly identical
          
          #### load up a dbh sample and resimulate each tree for 36 consecutive years
          tree.samp <- cage.dbh[[pix]][[j]] ## what initial trees are present in this simulator result
          deaths <- 0 ## keep track of how many times trees die in this resim
          delay <- sample(seq(delay.factor[1],
                              delay.factor[2]), 
                          size = length(tree.samp), 
                          replace=T) ## determine initial replanting delays
          
          ## diagnostic trackers
          track <- tree.samp[1]
          del <- delay[1]
          de <- 0 ### death ID
          for(e in 1:36){ ## test each tree for 36 years (2006-2040)
            ## resample delay time
            # ring <- delay==0 ## which clocks have worn down
            # delay[ring] <- sample(seq(delay.factor[1],
            #                           delay.factor[2]), 
            #                       size = length(delay[ring]), 
            #                       replace=T)
            # delay[!ring] <- delay[!ring]-1 ## count down on clocks with time remaining
            
            ### figure mortality and determine a kill list
            deathwatch <- mort(tree.samp)*mort.mod ## standard mortality probabilities
            deathwatch[tree.samp>=dbh.big] <- deathwatch[tree.samp>=dbh.big]*largemort.mod ## adjust mortality in larger trees
            kill.list <- rbinom(length(tree.samp), 1, deathwatch/100) ## randomly kill based on mortality rate
            # kill.list <- rep(1,length(tree.samp)) ## test: kill everything

            kill.list[tree.samp==0] <- 0 ## do not rekill the previously dead
            tree.samp[kill.list==1] <- 0 ## the dead are erased
            ## update delay clock for this cycle
            delay[kill.list==1] <- sample(seq(delay.factor[1],
                                              delay.factor[2]), 
                                          size = length(delay[kill.list==1]),
                                          replace=T) # start a clock for any tree killed
            delay[kill.list==0 & delay>0] <- delay[kill.list==0 & delay>0]-1 ## count down the delay clocks that had been set before
            de <- c(de, kill.list[1])
            
            ## grow up the survivors
            tree.samp[tree.samp>0] <- tree.samp[tree.samp>0]*(1+dbhg.pred(tree.samp[tree.samp>0]))
            
            ### now to determine which of the (previously) dead are replanted this year
            replanted <- rbinom(length(tree.samp[tree.samp==0]), 1, replant.factor) ## randoly replace only some of the dead trees
            tree.samp[replanted==1 & delay==0 & tree.samp==0] <- 5 ## for trees selected and without a delay
            del <- c(del, delay[1])
            ### the logic of this scheme is that delay indicates the number of years *after the death year* that there is 0 productivity
            ### trees lose all NPP in death year (growth calc happens before replant calc)
            ### delay range of 0-2 years implies a wait-out of 1-3 years before productivity happens again

            deaths <- deaths+sum(kill.list) ## count the dead
            track <- c(track, tree.samp[1]) ## record dbh history
            
          } ## end of loop for 36 year projector
          par(mfrow=c(1,3))
          plot(track, main="dbh")
          plot(de, main="death")
          plot(del, main="delay")
          data.frame(track, de, del, seq(0,36))
          deaths
          del
          deaths.track <- c(deaths.track, deaths) ## tally of total deaths in each pixel simulation
          dbh.sav[[pix]][[a]] <- tree.samp ## updated tree sample after morts + growth
        } # loop for number of resims in this pixel
        deaths.sav[[pix]] <- deaths.track ##  vector n=100 with number of deaths recorded in each pixel resimulation (i.e. how much replanting is needed in every pixel over 40 years)
        if(pix%%500 == 0){print(paste("chunk", o, "resimmed pix", pix))} ## give status updates
      } else{  ## if too few successful simulations for this pixel
        dbh.sav[[pix]][[a]] <- NA
        deaths.sav[[pix]] <- NA
      }
    } # loop for this pixel
    save(dbh.sav, file = paste("processed/boston/biom_street/dbh.resim.scenario.", scenario[s], ".", o, ".V1.sav", sep=""))
    save(deaths.track, file=paste("processed/boston/biom_street/death.track.scenario.", scenario[s], ".", o, ".V1.sav", sep=""))
  } ## loop for pixel chunk files
  print(paste("finished resim run scenario", scenario[s]))
}


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
