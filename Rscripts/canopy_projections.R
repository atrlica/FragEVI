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

## allometric equation E hardwoods
# b0.allo <- -2.48
# b1.allo <- 2.4835 ## these are eastern hardwood defaults
# biom.pred <- function(x){exp(b0.allo+(b1.allo*log(x)))}

### upgrade: apply urban-specific allometries
street.allo <- read.csv("docs/street.biometrics.csv") ## AG wood vol (m3) as b0*DBH(cm)^b1; multiply by density to get kg-biomass

### upgrade: use urban specific crown diameter measures to get canopy coverage
street.canopy <- as.data.table(read.csv("docs/RDS-2016-0005/Data/TS6_Growth_coefficients.csv"))
street.canopy <- street.canopy[Region=="NoEast", ]
street.canopy <- street.canopy[Predicts.component=="crown dia",]
genspec <- strsplit(as.character(street.canopy$Scientific.Name), " ")
gen <- unlist(lapply(genspec, "[[", 1))
street.canopy[,genus:=gen] ## 40 genera
keepers <- c("Acer platanoides", "Aesculus hippocastanum", "Fraxinus pennsylvanica",
             "Ginkgo biloba", "Gleditsia triacanthos", "Liquidambar styraciflua", 
             "Malus spp.", "Platanus x acerifolia", "Prunus serrulata", "Pyrus calleryana",
             "Quercus rubra", "Tilia cordata", "Ulmus americana", "Zelkova serrata")
street.canopy <- street.canopy[Scientific.Name%in%keepers,]
street.canopy$a <- as.numeric(as.character(street.canopy$a))
street.canopy$b <- as.numeric(as.character(street.canopy$b))
street.canopy$c <- as.numeric(as.character(street.canopy$c))
street.canopy$d <- as.numeric(as.character(street.canopy$d))
street.canopy$e <- as.numeric(as.character(street.canopy$e))
street.canopy <- as.data.frame(street.canopy)

### dbh~biomass growth equation (from nls street tree model)
# mm.a=0.2281686
# mm.b=(-1.1043188)
# dbhg.pred <- function(x){exp(mm.a+(mm.b*log(x)))} ## function for predicting next year's relative dbh growth from this year's dbh
## mm coefficients: had to hard code bc street data not avail on every machine

### upgrade: get a list of model coefficients to run for each repetition (dbh.delta~dbh.start)
# library(lme4)
# load("processed/mod.street.dbhdelta.me.sav")
# aaa <- summary(mod.street.dbhdelta.me)
# b0.rand <- rnorm(100, coef(aaa)[1,1], coef(aaa)[1,2])
# b1.rand <- rnorm(100, coef(aaa)[2,1], coef(aaa)[2,2])
# b2.rand <- rnorm(100, coef(aaa)[3,1], coef(aaa)[3,2])
b0.hard <- c(1.234, 8.178e-02) ### have to hard code shit because R version on cluster is old
b1.hard <- c(-1.989e-02, 3.693e-03)
b2.hard <- c(1.340e-04, 3.366e-05)
b0.rand <- rnorm(100, b0.hard[1], b0.hard[2])
b1.rand <- rnorm(100, b1.hard[1], b1.hard[2])
b2.rand <- rnorm(100, b2.hard[1], b2.hard[2])

## que up the identified street tree planting polygons, filter, and figure out how much land we are dealing with

# ### now need to find out where these places are
# box@data <- merge(box@data, box.l.dat, by="OBJECTID")
# head(box@data)
# unique(box@data$num.plant)
# # biom <- raster("processed/boston/bos.biom30m.tif")
# # box.r <- rasterize(box, biom, field="num.plant", fun=sum)
# # writeRaster(box.r, "processed/boston/bos.streetplanters.tif", format="GTiff", overwrite=T)
# ## so we can randomly pick places and slowly fill them in in the "expand" scenario

## load up ancillary map data
library(raster)
library(rgdal)
biom <- raster("processed/boston/bos.biom30m.tif")
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- crop(biom, aoi)
can <- raster("processed/boston/bos.can30m.tif")
isa <- raster("processed/boston/bos.isa30m.tif")
biom.dat <- as.data.table(as.data.frame(biom))
aoi.dat <- as.data.table(as.data.frame(aoi))
can.dat <- as.data.table(as.data.frame(can))
isa.dat <- as.data.table(as.data.frame(isa))
lulc <- raster("processed/boston/bos.lulc30m.lumped.tif")
lulc.dat <- as.data.table(as.data.frame(lulc))
biom.dat <- cbind(biom.dat, aoi.dat, can.dat, isa.dat, lulc.dat)
biom.dat[,pix.ID:=seq(1:dim(biom.dat)[1])]
nonfor <- biom.dat[bos.lulc30m.lumped!=1 & bos.biom30m<20000, pix.ID] ### identify pixel ID's that are nonforest

######
### set scenario options 
scenario <- c("BAU", "highmort", "lowreplant", "oldies", "lowmort", "slowreplant", "expand")
scenario <- c("BAU")

## record specific parameter sets
resim.vers <- 2 ## what are we labeling this round of resims?
vers <- 4 ## what simulator results version are we dealing with?
### BAU/default factors
npp.quant.range <- c(0.001, 0.999) ## what dbh samples to draw from
mort.mod <- 1 ## modification to mortality rate (scenario specific)
default.sizecutoff <- 10000 ## size of tree dbh over which to screw with mortality
largemort.mod <- 1 ## multiplier for mortalities in large trees
replant.factor <- 1 ## rate of replanting
dbh.big <- default.sizecutoff
default.delay <- c(0,0) ## low/hi on years of delay from death to replanting
delay.factor <- default.delay

## set options for scenario processing
highmort.mortfactor <- 1.25
lowmort.mortfactor <- 0.5
lowreplant.factor <- 0.5 ## what fraction of morts are replanted in "lowreplant"
oldies.mortfactor <- 0.5 ## how much to reduce mortalities in large trees
oldies.sizecutoff <- 40 ## how to define "big" trees
slowreplant.delay <- c(0,2) ## 1 to 3 year delay between death and regrowth starting
# expand.timeline <- 10 ### how many years to implement the expand scenario
# expand.rate.go <- (sum(box.l.dat$num.plant, na.rm=T)/expand.timeline)/length(nonfor) ## this is the per non-forest pixel annual rate of planting needed to fill out the city
expand.rate <- 0

# record run parameters to text file
# params.list <- list(c(scenario, resim.vers, vers,
#                     highmort.mortfactor, lowmort.mortfactor, 
#                     npp.quant.range, lowreplant.factor, oldies.sizecutoff, 
#                     oldies.mortfactor))
sink(paste0("processed/boston/biom_street/resim_params_list_V", resim.vers, ".txt"))
cat(c("Resimulation version: ", resim.vers), "\n")
cat(c("Simulator results version: ", vers), "\n")
cat(c("scenarios = ", scenario), "\n")
cat(c("NPP selection quantiles = ", npp.quant.range), "\n")
cat(c("highmort mortality factor = ", highmort.mortfactor), "\n")
cat(c("lowmort mortality factor = ", lowmort.mortfactor), "\n")
cat(c("oldies mortality factor = ", oldies.mortfactor), "\n")
cat(c("oldies size cutoff = ", oldies.sizecutoff), "\n")
cat(c("lowreplant replanting rate = ", lowreplant.factor), "\n")
cat(c("slowreplant delay time (yrs) = ", slowreplant.delay), "\n")
sink()

## loops for simulating growth+mortality in pixel dbh samples
## pull in all the simulator samples of tree dbh in every pixel
obj.list <- list.files("processed/boston/biom_street/")
obj.list <- obj.list[grep(obj.list, pattern=paste("dbh.street.v", vers, sep=""))]
obj.list <- sub('.*weighted\\.', '', obj.list)
obj.list <- (sub('\\..*', '', obj.list))

## loop the simulated pixel chunks
## each pixel*scenario takes ~3.5s to run -- a single chunk will take ~2hrs per scenario, ~33hr per map for a single scenario
for(s in 1:length(scenario)){
  ### modify simulator factors according to scenario construction
  if(scenario[s]=="highmort"){
    mort.mod <- highmort.mortfactor
    dbh.big <- default.sizecutoff 
    largemort.mod <- 1
    replant.factor <- 1 ## rate of replanting
    delay.factor <- default.delay ## low/hi on years of delay from death to replanting
    expand.rate <- 0
  }
  if(scenario[s]=="lowmort"){
    mort.mod <- lowmort.mortfactor
    dbh.big <- default.sizecutoff 
    largemort.mod <- 1
    replant.factor <- 1 ## rate of replanting
    delay.factor <- default.delay ## low/hi on years of delay from death to replanting
    expand.rate <- 0
  }
  if(scenario[s]=="oldies"){
    mort.mod <- 1
    dbh.big <- oldies.sizecutoff 
    largemort.mod <- oldies.mortfactor
    replant.factor <- 1 ## rate of replanting
    delay.factor <- default.delay ## low/hi on years of delay from death to replanting
    expand.rate <- 0
  }
  if(scenario[s]=="lowreplant"){
    replant.factor=0.5
    mort.mod <- 1
    dbh.big <- default.sizecutoff 
    largemort.mod <- 1
    delay.factor <- default.delay ## low/hi on years of delay from death to replanting
    expand.rate <- 0
  }
  if(scenario[s]=="slowreplant"){
    delay.factor <- slowreplant.delay
    mort.mod <- 1
    dbh.big <- default.sizecutoff 
    largemort.mod <- 1
    replant.factor <- 1 ## rate of replanting
    delay.factor <- c(0,3)
    expand.rate <- 0
  }
  if(scenario[s]=="expand"){
    mort.mod <- 1
    dbh.big <- default.sizecutoff 
    largemort.mod <- 1
    replant.factor <- 1 ## rate of replanting
    delay.factor <- default.delay
    expand.rate <- expand.rate.go
  }
  print(paste("starting resim run scenario", scenario[s]))
  ### avoid processing shit that's already done
  already <- list.files("processed/boston/biom_street")
  already <- already[grep(already, pattern=".resim.")]
  already <- already[grep(already, pattern=paste(".", scenario[s], ".", sep=""))]
  already <- already[grep(already, pattern=paste0(".V", resim.vers))]
  already <- strsplit(already, split = "[.]")
  already <- sapply(already, "[[" ,5)
  
  hitlist <- as.character(obj.list[!(obj.list %in%  already)])
  hitlist[hitlist=="100000"] <- "1e+05"
  
  for(o in hitlist){
    ### initialize read-in files and save dumps
    load(paste("processed/boston/biom_street/dbh.street.v", vers, ".weighted.", o, ".sav", sep="")) # cage.dbh
    load(paste("processed/boston/biom_street/biom.sim.street.v", vers, ".weighted.", o, ".sav", sep="")) # cage.biom.sim
    load(paste("processed/boston/biom_street/ann.npp.street.v", vers, ".weighted.", o, ".sav", sep="")) # cage.ann.npp
    load(paste("processed/boston/biom_street/index.track.street.v", vers, ".weighted.", o, ".sav", sep="")) #index.track
    load(paste("processed/boston/biom_street/genus.street.v", vers, ".weighted.", o, ".sav", sep="")) ## comes in as "cage.wts" object
    deaths.sav <- list()
    npp.box <- list()## initialize a container for the npp and tree number results
    biom.box <- list() ## keep track of biomass through simulation
    expand.plant.num <- list() ## keeping track of the number of new plants put in
    num.box <- list() ## keep track of the number of live stems
    
    # procset <- sample(1:length(cage.dbh), size=50) ## test a random subset
    procset <- 1:length(cage.dbh) ## FULL PROCESS

    for(pix in procset){ ## the number of pixels in this sim chunk
      print(paste("working on", pix))
      if(length(cage.biom.sim[[pix]])>=40){ ## only process if enough simulations successfully completed in this pixel
        deaths.sav[[pix]] <- integer()
        npp.box[[pix]] <- list()
        biom.box[[pix]] <- list()
        expand.plant.num[[pix]] <- integer()
        num.box[[pix]] <- list()
        
        for(a in 1:100){   ### resimulate every pixel 'a' number of times
          ## can select dbh populations based on proximity to median simulated biomass...
          # biom.lims <- quantile(cage.biom.sim[[pix]], probs=c(0.40, 0.6)) ## figure out which of the simulations to draw and modify
          # j <- sample(which(cage.biom.sim[[pix]]>=biom.lims[1] & cage.biom.sim[[pix]]<=biom.lims[2]),1) ## get a random dbh sample from the the middle 10% of simulator results close to the target biomass
          
          ## ...or can select dbh populations based on how close they are to median npp (not 100% overlapping)
          npp.lims <- quantile(cage.ann.npp[[pix]], probs=npp.quant.range) ## restrict which of the simulations to draw and modify
          j <- which(cage.ann.npp[[pix]]>=npp.lims[1] & cage.ann.npp[[pix]]<=npp.lims[2])
          if(length(j)>=4){ ## if there's not enough dbh samples that meet your criteria...
            j <- sample(j,1) ## get a random dbh sample selected range of simulator results close to the target biomass
          }else{j <- sample(1:length(cage.biom.sim[[pix]]), 1)} ## or just sample the whole simulation collection if all simulations are nearly identical
          
          #### load up a dbh sample and resimulate each tree for 36 consecutive years
          tree.samp <- cage.dbh[[pix]][[j]] ## what initial trees are present in this simulator result
          
          ## tree expansion scenario
          ### this is imperfect -- but assume a small amount of new street buffer planting space is available in this spot and sim new trees
          expand.track <- rbinom(n = 10, prob=expand.rate, size=1) ## give this tree a tree plant flag for the next ten years that on average across the map gets the right number of plantings in each pixel
          
          npp.track <- numeric() ## annual tally of npp in this resim as you go up in years
          biom.track <- numeric() ## annual tally of biomass in this resim
          num.track <- integer() ## annual tally of number of trees in this resim
          newplants <- 0 ## track number of new trees planted over course of 36 years
          can.track <- numeric() ### canopy coverage in this simulation
          # delay <- sample(seq(delay.factor[1],
          #                     delay.factor[2]), 
          #                 size = length(tree.samp), 
          #                 replace=T) ## determine initial replanting delays
          delay <- rep(0, length(tree.samp))
          nogenus <- sample(1:14, size=1) ## randomly assign a canopy equation for the handful that don't match
          for(e in 1:36){ ## simulate 36 successive years of growth/mortality/replanting/plant expansion
            if(scenario[s]=="expand" & e<=10 & expand.track[e]==1 & index.track[pix]%in%nonfor){ 
              tree.samp <- c(tree.samp, 5) ## add a 5cm tree if the expand.track wins one this year (within the first 10 years)
              newplants <- newplants+1
              delay <- c(delay, 0)
            }
          deaths <- 0 ## keep track of how many times trees die in this pix resim
          
          ### figure mortality and determine a kill list
          deathwatch <- mort(tree.samp)*mort.mod ## mortality liklihood for each stem
          deathwatch[tree.samp>=dbh.big] <- deathwatch[tree.samp>=dbh.big]*largemort.mod ## adjust mortality in larger trees
          kill.list <- rbinom(length(tree.samp), 1, deathwatch/100) ## randomly kill based on mortality rate
          # kill.list <- rep(1,length(tree.samp)) ## test: kill everything

          kill.list[tree.samp==0] <- 0 ## do not rekill the previously dead
          tree.samp[kill.list==1] <- 0 ## the dead are nullified
          
          ## determine a delay clock for this stems killed this cycle and update old delay clocks
          delay[kill.list==1] <- sample(seq(delay.factor[1],
                                            delay.factor[2]), 
                                        size = length(delay[kill.list==1]),
                                        replace=T) # start a clock for any tree killed
          delay[kill.list==0 & delay>0] <- delay[kill.list==0 & delay>0]-1 ## count down the delay clocks that had been set before
          
          ## grow up the survivors using estimated growth regression
          # tree.samp <- tree.samp[tree.samp>0]*(1+dbhg.pred(tree.samp[tree.samp>0])) ## static realization of growth model
          tree.samp[tree.samp>0] <- tree.samp[tree.samp>0]+(b0.rand[a]+(b1.rand[a]*tree.samp[tree.samp>0])+(b2.rand[a]*tree.samp[tree.samp>0]^2))
          
          ### upgrade: urban specific allometries to determine biomass change
          tmp.dbh0 <- tree.samp[tree.samp>0] ## exclude dead trees from dbh/biomass change
          tmp.genus <- cage.genus[[pix]][[j]][tree.samp>0] ## genera of dbh collection

          ### fast match the genus to the specific allometric if possible, or divert to general equation at row 8
          tmp.biom0 <- street.allo[match(tmp.genus, street.allo$genus, nomatch=8), "b0"]*(tmp.dbh0^street.allo[match(tmp.genus, street.allo$genus, nomatch=8), "b1"])*street.allo[match(tmp.genus, street.allo$genus, nomatch=8), "dens"]
          tmp.dbh1 <- tmp.dbh0+(b0.rand[a]+(b1.rand[a]*tmp.dbh0)+(b2.rand[a]*(tmp.dbh0^2))) ## grow the dbh to time 1, no growth in dead trees
          tmp.biom1 <- street.allo[match(tmp.genus, street.allo$genus, nomatch=8), "b0"]*(tmp.dbh1^street.allo[match(tmp.genus, street.allo$genus, nomatch=8), "b1"])*street.allo[match(tmp.genus, street.allo$genus, nomatch=8), "dens"]
          if(length(tmp.biom0)==0){tmp.biom0 <- 0; tmp.biom1 <- 0} ## if everything is dead
          
          ### same thing for canopy coverage history
          ### get vector of equation types to use first
          eq.form <- street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "EqName"]
          ## add up total canopy area applying correct equation form to each
          tmp.can1 <- sum(((((eq.form=="quad")*(street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "a"]+
                                     (street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "b"]*tmp.dbh1)+
                                     (street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "c"]*tmp.dbh1^2)))/2)^2)*pi, na.rm=T)+
            sum(((((eq.form=="cub")*(street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "a"]+
                                       (street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "b"]*tmp.dbh1)+
                                       (street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "c"]*tmp.dbh1^2)+
                                      (street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "d"]*tmp.dbh1^3)))/2)^2)*pi, na.rm=T)+
            sum(((((eq.form=="lin")*(street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "a"]+
                                        (street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "b"]*tmp.dbh1)))/2)^2)*pi, na.rm=T)+
            sum((((eq.form=="loglogw1")*exp(street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "a"]+
                                       (street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "b"]*
                                          log(log(tmp.dbh1+1)+
                                                street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "c"]/2)))/2)^2)*pi, na.rm=T)
                                       
          
          ### update productivity and biomass for track on this year
          npp.track <- c(npp.track, sum(tmp.biom1-tmp.biom0))
          biom.track <- c(biom.track, sum(tmp.biom1))
          can.track <- c(can.track, tmp.can1)
          
          ### now determine which of the (previously) dead are replanted this year
          if(sum(tree.samp==0)>0){ ## if there are dead trees here          
            replanted <- rbinom(length(tree.samp[tree.samp==0]), 1, replant.factor) ## randomly replace only some of the dead trees
            if(sum(tree.samp==0 & delay==0)>0){
              tree.samp[tree.samp==0 & delay==0] <- replanted*5 ## for trees selected and without a delay
              }
          }
          ### the logic of this scheme is that delay indicates the number of years *after the death year* that there is 0 productivity
          ### trees lose all NPP in death year (growth calc happens before replant calc)
          ### delay range of 0-2 years implies a wait-out of 1-3 years before productivity happens again

          deaths <- deaths+sum(kill.list) ## count the dead
          num.track <- c(num.track, length(tree.samp[tree.samp>0])) ## count the living
          # track <- c(track, tree.samp[1]) ## record dbh history
          # de <- c(de, kill.list[1]) ## record death history
          #, delay[1]) ## record delay history
        } ## end of loop for 36 year projector
          # par(mfrow=c(1,3))
          # plot(track, main="dbh")
          # plot(de, main="death")
          # plot(del, main="delay")
          # data.frame(track, de, del, seq(0,36))
          
          # dbh.sav[[pix]][[a]] <- tree.samp ## updated tree sample after all 36 years of morts + growth
          npp.box[[pix]][[a]] <- npp.track ## use matrix(unlist(npp.box[[pix]]), nrow=100, byrow=T)
          biom.box[[pix]][[a]] <- biom.track
          expand.plant.num[[pix]] <- c(expand.plant.num[[pix]], newplants) ## eventually produces vector of length 100
          deaths.sav[[pix]] <- c(deaths.sav[[pix]], deaths)
          num.box[[pix]][[a]] <- num.track
          can.box[[pix]][[a]] <- can.track
        } # loop for number of resims in this pixel
        if(pix%%500 == 0){print(paste("chunk", o, "resimmed pix", pix))} ## give status updates
      } else{  ## if too few successful simulations for this pixel
        npp.box[[pix]] <- NA
        biom.box[[pix]] <- NA
        expand.plant.num[[pix]] <- NA
        deaths.sav[[pix]] <- NA
        num.box[[pix]] <- NA
        can.box[[pix]][[a]] <- NA
      }

    } # end loop for this pixel   
    
    save(deaths.sav, file=paste("processed/boston/biom_street/death.track.scenario.", scenario[s], ".", o, ".V", resim.vers, ".sav", sep=""))
    save(npp.box, file=paste("processed/boston/biom_street/npp.track.scenario.", scenario[s], ".", o, ".V", resim.vers, ".sav", sep=""))
    save(biom.box, file=paste("processed/boston/biom_street/biom.track.scenario.", scenario[s], ".", o, ".V", resim.vers, ".sav", sep=""))
    save(expand.plant.num, file=paste("processed/boston/biom_street/expand.plant.num.scenario.", scenario[s], ".", o, ".V", resim.vers, ".sav", sep=""))
    save(num.box, file=paste("processed/boston/biom_street/num.track.scenario.", scenario[s], ".", o, ".V", resim.vers, ".sav", sep=""))
  } ## loop for pixel chunk files
  print(paste("finished resim run scenario", scenario[s]))
}



# ###
# ### process resim results by scenario
# ##########
# b0 <- -2.48
# b1 <- 2.4835 ## these are eastern hardwood defaults
# biom.pred <- function(x){exp(b0+(b1*log(x)))}
# load("processed/mod.street.dbhdelta.sav") ## delta dbh predictor
# dbh.b0 <- summary(mod.street.dbhdelta)$coefficients[1]
# dbh.b1 <- summary(mod.street.dbhdelta)$coefficients[2]
# dbh.b2 <- summary(mod.street.dbhdelta)$coefficients[3]
# dbh.grow <- function(x){dbh.b0+(x*dbh.b1)+((x^2)*dbh.b2)}
# biom.grow <- function(x){biom.pred(x+dbh.grow(x))-biom.pred(x)}
# 
# ## get a general list of the resim results and the corresponding pixel index
# obj.list <- list.files("processed/boston/biom_street/")
# obj.list <- obj.list[grep(obj.list, pattern=paste("scenario",sep=""))]
# obj.list <- obj.list[grep(obj.list, pattern="dbh")]
# index.list <- list.files("processed/boston/biom_street/")
# index.list <- index.list[grep(index.list, pattern="index")]
# index.list <- index.list[grep(index.list, pattern="v4")]
# 
# ### BAU scenario
# scen <- "BAU"
# dbh.list <- obj.list[grep(obj.list, pattern=paste(scen))]
# dbh.list <- sub('.*BAU\\.', '', dbh.list)
# dbh.list <- (sub('\\..*', '', dbh.list))
# hitlist <- as.character(dbh.list)
# hitlist[hitlist=="100000"] <- "1e+05"
# preamb <- "processed/boston/biom_street/"
# dump <- data.frame()
# for(o in 1:length(hitlist)){
#   load(paste0(preamb,"dbh.resim.scenario-0.BAU.", hitlist[o], ".V1.sav")) ## as dbh.sav
#   load(paste0(preamb, "index.track.street.v4.weighted.", hitlist[o], ".sav")) ## as index.track
#   load(paste0(preamb, "dbh.street.v4.weighted.", hitlist[o], ".sav", sep="")) # as cage.dbh
#   
#   pix.track <- integer()
#   growth.track <- numeric()
#   for(pix in 1:length(dbh.sav)){ ## loop each pixel
#     if(length(cage.dbh[[pix]])>=40){ ## we only resimmed if we got enough successful retreivals
#       growth.track <- c(growth.track, median(sapply(sapply(dbh.sav[[pix]], biom.grow),sum), na.rm=T)) ## the median of the sums of each pixel sim, translated dbh to biom gain
#       pix.track <- c(pix.track, index.track[pix])
#     } else{
#       growth.track <- c(growth.track, NA)
#       pix.track <- c(pix.track, index.track[pix])
#     }
#   }
#   print(paste("finished chunk", hitlist[o], "of BAU resim"))
#   tmp <- cbind(pix.track, growth.track)
#   dump <- rbind(dump, tmp)
# }
# names(dump) <- c("pix.ID", "BAU.2040med")
# write.csv(dump, paste("processed/scenario.BAU.med2040.results.csv"))
# 
# 
# ### all other scenarios
# scen <- c("highmort", "lowmort", "lowreplant", "oldies", "slowreplant")
# for(s in 1:length(scen)){
#   dbh.list <- obj.list[grep(obj.list, pattern=paste(scen[s]))]
# #   dbh.list <- sub('.*scenario\\.', '', dbh.list)
#   dbh.list <- sub(paste0(".*", scen[s], "\\."), '', dbh.list)
#   dbh.list <- (sub('\\..*', '', dbh.list))
#   hitlist <- as.character(dbh.list)
#   hitlist[hitlist=="100000"] <- "1e+05"
#   preamb <- "processed/boston/biom_street/"
#   dump <- data.frame()
#   for(o in 1:length(hitlist)){
#     load(paste0(preamb,"dbh.resim.scenario.", scen[s], ".", hitlist[o], ".V1.sav")) ## as dbh.sav
#     load(paste0(preamb, "index.track.street.v4.weighted.", hitlist[o], ".sav")) ## as index.track
#     load(paste0(preamb, "dbh.street.v4.weighted.", hitlist[o], ".sav", sep="")) # as cage.dbh
#     
#     pix.track <- integer()
#     growth.track <- numeric()
#     for(pix in 1:length(dbh.sav)){ ## loop each pixel
#       if(length(cage.dbh[[pix]])>=40){ ## we only resimmed if we got enough successful retreivals
#         growth.track <- c(growth.track, median(sapply(sapply(dbh.sav[[pix]], biom.grow),sum), na.rm=T)) ## the median of the sums of each pixel sim, translated dbh to biom gain
#         pix.track <- c(pix.track, index.track[pix])
#       } else{
#         growth.track <- c(growth.track, NA)
#         pix.track <- c(pix.track, index.track[pix])
#       }
#     }
#     print(paste("finished chunk", hitlist[o], "of", scen[s], "resim"))
#     tmp <- cbind(pix.track, growth.track)
#     dump <- rbind(dump, tmp)
#   }
#   names(dump) <- c("pix.ID", paste0(scen[s], ".2040med"))
#   write.csv(dump, paste0("processed/scenario.", scen[s], ".med2040.results.csv"))
# }
# 
# ### put them back together with the map data
# lulc <- raster("processed/boston/bos.lulc30m.lumped.tif")
# extent(lulc)
# biom <- raster("processed/boston/bos.biom30m.tif")
# extent(biom)
# biom <- crop(biom, lulc)
# biom <- as.data.table(as.data.frame(biom))
# lulc <- as.data.table(as.data.frame(lulc))
# lulc[,pix.ID:=seq(1:dim(lulc)[1])]
# lulc[,bos.biom30m:=biom]
# aoi <- as.data.table(as.data.frame(raster("processed/boston/bos.aoi30m.tif")))
# lulc[,aoi:=aoi]
# bau <- read.csv(paste("processed/scenario.BAU.med2040.results.csv"))
# bau$X <- NULL
# lowreplant <- read.csv(paste("processed/scenario.lowreplant.med2040.results.csv"))
# lowreplant$X <- NULL ## this is coming in doubled for reasons unknown
# lowreplant <- lowreplant[!(duplicated(lowreplant$pix.ID)),]
# highmort <- read.csv(paste("processed/scenario.highmort.med2040.results.csv"))
# highmort$X <- NULL
# lowmort <- read.csv(paste("processed/scenario.lowmort.med2040.results.csv"))
# lowmort$X <- NULL
# oldies <- read.csv(paste("processed/scenario.oldies.med2040.results.csv"))
# oldies$X <- NULL
# slowreplant <- read.csv(paste("processed/scenario.slowreplant.med2040.results.csv"))
# slowreplant$X <- NULL
# 
# resim.dat <- merge(lulc, bau, by="pix.ID", all.x=T)
# # resim.dat <- resim.dat[, c(1,2,4)]
# resim.dat <- merge(resim.dat, lowreplant, by="pix.ID", all.x=T)
# # resim.dat <- resim.dat[, c(1,2,3,5)]
# resim.dat <- merge(resim.dat, highmort, by="pix.ID", all.x=T)
# # resim.dat <- resim.dat[, c(1,2,3,4,6)]
# resim.dat <- merge(resim.dat, lowmort, by="pix.ID", all.x=T)
# # resim.dat <- resim.dat[, c(1,2,3,4,5,7)]
# resim.dat <- merge(resim.dat, oldies, by="pix.ID", all.x=T)
# # resim.dat <- resim.dat[, c(1,2,3,4,5,6,8)]
# resim.dat <- merge(resim.dat, slowreplant, by="pix.ID", all.x=T)
# # resim.dat <- resim.dat[, c(1,2,3,4,6,7,9)]
# resim.dat[aoi>800 & bos.lulc30m.lumped!=1 & bos.biom30m<20000, .(median(BAU.2040med, na.rm=T),
#                                                                  median(lowreplant.2040med, na.rm=T),
#                                                                  median(highmort.2040med, na.rm=T),
#                                                                  median(lowmort.2040med, na.rm=T),
#                                                                  median(oldies.2040med, na.rm=T),
#                                                                  median(slowreplant.2040med, na.rm=T))]
# 
# resim.dat[aoi>800 & bos.lulc30m.lumped!=1 & bos.biom30m<20000, .(sum(BAU.2040med, na.rm=T)/2000,
#                                                                  sum(lowreplant.2040med, na.rm=T)/2000,
#                                                                  sum(highmort.2040med, na.rm=T)/2000,
#                                                                  sum(lowmort.2040med, na.rm=T)/2000,
#                                                                  sum(oldies.2040med, na.rm=T)/2000,
#                                                                  sum(slowreplant.2040med, na.rm=T)/2000)]
# 
# 
# resim.dat[aoi>800 & bos.lulc30m.lumped==3 & bos.biom30m<20000, .(median(BAU.2040med, na.rm=T),
#                                                                  median(lowreplant.2040med, na.rm=T),
#                                                                  median(highmort.2040med, na.rm=T),
#                                                                  median(lowmort.2040med, na.rm=T),
#                                                                  median(oldies.2040med, na.rm=T),
#                                                                  median(slowreplant.2040med, na.rm=T))]
# 
# resim.dat[aoi>800 & bos.lulc30m.lumped%in%c(3,4) & bos.biom30m<20000, .(sum(BAU.2040med, na.rm=T)/2000,
#                                                                  sum(lowreplant.2040med, na.rm=T)/2000,
#                                                                  sum(highmort.2040med, na.rm=T)/2000,
#                                                                  sum(lowmort.2040med, na.rm=T)/2000,
#                                                                  sum(oldies.2040med, na.rm=T)/2000,
#                                                                  sum(slowreplant.2040med, na.rm=T)/2000)]
# 
# resim.dat[aoi>800 & bos.lulc30m.lumped==3 & bos.biom30m<20000, .(quantile(((BAU.2040med/2000)/aoi)*1E4, na.rm=T, probs=c(0.05, 0.5, 0.95)),
#                                                                  quantile(((lowreplant.2040med/2000)/aoi)*1E4, na.rm=T, probs=c(0.05, 0.5, 0.95)),
#                                                                  quantile(((highmort.2040med/2000)/aoi)*1E4, na.rm=T, probs=c(0.05, 0.5, 0.95)),
#                                                                  quantile(((lowmort.2040med/2000)/aoi)*1E4, na.rm=T, probs=c(0.05, 0.5, 0.95)),
#                                                                  quantile(((oldies.2040med/2000)/aoi)*1E4, na.rm=T, probs=c(0.05, 0.5, 0.95)),
#                                                                  quantile(((slowreplant.2040med/2000)/aoi)*1E4, na.rm=T, probs=c(0.05, 0.5, 0.95)))]
# 
# 
# resid <- resim.dat[aoi>800 & bos.lulc30m.lumped %in%c(3,4) & bos.biom30m<20000,] ## 54k
# library(tidyr)
# library(ggplot2)
# resid.long <- gather(resid, scenario, biom.growth, BAU.2040med:slowreplant.2040med, factor_key=TRUE)
# head(resid.long)
# unique(resid.long$bos.lulc30m.lumped)
# unique(resid.long$scenario)
# resid.long <- as.data.table(resid.long)
# resid.long[,MgC.ha.yr:=((biom.growth/2000)/aoi)*1E4]
# ggplot(resid.long, aes(x=scenario, y=MgC.ha.yr))+
#   geom_boxplot(outlier.shape=NA)+
#   scale_y_continuous(limits=c(0, 2))+
#   ylab("Annual C uptake (MgC/ha/yr)")

## how many trees do you get in each re-sim  
# num.resim.trees <- unlist(lapply(dbh.sav[[pix]], length))
# num.sim.trees <- unlist(lapply(cage.dbh[[pix]], length))
# hist(num.sim.trees)
# hist(num.resim.trees)
# 
# plot(unlist(lapply(cage.dbh[[400]], median)), unlist(lapply(dbh.sav[[400]], median)))
# abline(a=0, b=1)
