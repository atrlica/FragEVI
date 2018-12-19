# library(raster)
library(data.table)
# setwd("/projectnb/buultra/atrlica/FragEVI")

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
nogenus <- 1## resolved: Acer platanoides is the most common spp -- for the handful of spp not in this list, we will use the A. platanoides canopy equation




### dbh~biomass growth equation (from nls street tree model)
# mm.a=0.2281686
# mm.b=(-1.1043188)
# dbhg.pred <- function(x){exp(mm.a+(mm.b*log(x)))} ## function for predicting next year's relative dbh growth from this year's dbh
## mm coefficients: had to hard code bc street data not avail on every machine

### upgrade: get a list of model coefficients to run for each repetition (dbh.delta~dbh.start)
library(lme4)
load("processed/mod.street.dbhdelta.me.sav")
aaa <- summary(mod.street.dbhdelta.me)
b0.rand <- rnorm(100, coef(aaa)[1,1], coef(aaa)[1,2])
b1.rand <- rnorm(100, coef(aaa)[2,1], coef(aaa)[2,2])
b2.rand <- rnorm(100, coef(aaa)[3,1], coef(aaa)[3,2])
b0.hard <- c(1.234, 8.178e-02) ### have to hard code shit because R version on cluster is old
b1.hard <- c(-1.989e-02, 3.693e-03)
b2.hard <- c(1.340e-04, 3.366e-05)
b0.rand <- rnorm(100, b0.hard[1], b0.hard[2])
b1.rand <- rnorm(100, b1.hard[1], b1.hard[2])
b2.rand <- rnorm(100, b2.hard[1], b2.hard[2])
b0.rand <- b0.hard[1]
b1.rand <- b1.hard[1]
b2.rand <- b2.hard[1]

### expand scenario ancillary data
## que up the identified street tree planting polygons, filter, and figure out how much land we are dealing with
# ### still need to figure how to locate planting areas in space viz the biomass grid
box <- as.data.table(read.csv("processed/boston/plant10m_MBG.csv")); dim(box) ## 78618
summary(box$Shape_Area); hist(box$Shape_Area)
box <- box[Shape_Area<1000,]; dim(box) ##78034 ## eliminate a fair chunk of the weird giant boxes
box[,num.trees:=((MBG_Width-1)%/%8)+1] ## this gives us 2 trees per 8m planter length with a little buffer at either end
box[,sum(num.trees)] ## 79k trees can be planted in this population of boxes

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
nonfor <- biom.dat[bos.aoi30m>800 & bos.lulc30m.lumped!=1 & bos.biom30m<20000, pix.ID] ### identify pixel ID's that are nonforest

######
### set scenario options 
# scenario <- c("BAU", "highmort", "lowreplant", "oldies", "lowmort", "slowreplant", "expand")
scenario <- c("BAU", "oldies", "expand")

## record specific parameter sets
resim.vers <- 4 ## what are we labeling this round of resims?
## 1 was prelim, 2 included full 100x pix resim per year, 3 is model mean single resim per year
vers <- 5 ## what simulator results version are we dealing with?

### default factors
npp.quant.range <- c(0.40, 0.60) ## what dbh samples to draw from
mort.mod <- 1 ## modification to mortality rate (scenario specific)
default.sizecutoff <- 10000 ## size of tree dbh over which to screw with mortality
largemort.mod <- 1 ## multiplier for mortalities in large trees
replant.factor <- 1 ## rate of replanting
dbh.big <- default.sizecutoff
default.delay <- c(0,2) ## low/hi on years of delay from death to replanting
delay.factor <- default.delay

## set options for scenario processing
highmort.mortfactor <- 1.25
lowmort.mortfactor <- 0.5
lowreplant.factor <- 0.5 ## what fraction of morts are replanted in "lowreplant"
oldies.mortfactor <- 0.5 ## how much to reduce mortalities in large trees
oldies.sizecutoff <- 40 ## how to define "big" trees
slowreplant.delay <- c(0,2) ## 1 to 3 year delay between death and regrowth starting
expand.timeline <- 10 ### how many years to implement the expand scenario
expand.rate.go <- (box[,sum(num.trees)]/expand.timeline)/length(nonfor) ## this is the per non-forest pixel annual rate of planting needed to fill out the city (assuming that only a small number of the replant boxes fall in forest pixels)

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
cat(c("default replant delay = ", delay.factor), "\n")
cat(c("highmort mortality factor = ", highmort.mortfactor), "\n")
cat(c("lowmort mortality factor = ", lowmort.mortfactor), "\n")
cat(c("oldies mortality factor = ", oldies.mortfactor), "\n")
cat(c("oldies size cutoff = ", oldies.sizecutoff), "\n")
cat(c("lowreplant replanting rate = ", lowreplant.factor), "\n")
cat(c("slowreplant delay time (yrs) = ", slowreplant.delay), "\n")
cat(c("expand per-pixel prob of tree increase/yr = ", expand.rate.go), "\n")
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
  if(scenario[s]=="BAU"){
    mort.mod <- 1 ## modification to mortality rate (scenario specific)
    replant.factor <- 1 ## rate of replanting
    dbh.big <- default.sizecutoff
    largemort.mod <- 1
    delay.factor <- default.delay ## low/hi on years of delay from death to replanting
    expand.rate <- 0
  }
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
  already <- already[grep(already, pattern=paste0(".V", resim.vers))]
  already <- already[grep(already, pattern="npp.track.")]
  already <- already[grep(already, pattern=paste(".", scenario[s], ".", sep=""))]
  already <- strsplit(already, split = "[.]")
  already <- sapply(already, "[[" ,5)
  
  # hitlist <- as.character(obj.list[!(obj.list %in%  already)])
  hitlist <- as.character(obj.list)
  hitlist[hitlist=="100000"] <- "1e+05"
  if(length(hitlist)==0){print(paste(scenario[s], "already processed"))}
  
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
    can.box <- list()
    pixID.track <- integer()
    # procset <- sample(1:length(cage.dbh), size=50) ## test a random subset
    procset <- 1:length(cage.dbh) ## FULL PROCESS
    procset <- procset[(index.track%in%nonfor)] ### only resim the non-forest
    print(paste("resimming chunk", o))
    
    for(pix in 1:length(procset)){ ## the number of pixels in this sim chunk
      # print(paste("working on pixel", index.track[procset[pix]]))
      pixID.track <- c(pixID.track, index.track[procset[pix]])
      deaths.sav[[pix]] <- integer()
      npp.box[[pix]] <- list()
      biom.box[[pix]] <- list()
      expand.plant.num[[pix]] <- integer()
      num.box[[pix]] <- list()
      can.box[[pix]] <- list()

      if(length(cage.biom.sim[[procset[pix]]])>=40){ ## only process if enough simulations successfully completed in this pixel

        for(a in 1:1){   ### resimulate every pixel 'a' number of times

          # can select dbh populations based on proximity to median simulated biomass...
          biom.lims <- quantile(cage.biom.sim[[procset[pix]]], probs=c(0.40, 0.6)) ## figure out which of the simulations to draw and modify
          j <- sample(which(cage.biom.sim[[procset[pix]]]>=biom.lims[1] & cage.biom.sim[[procset[pix]]]<=biom.lims[2]),1) ## get a random dbh sample from the the middle 10% of simulator results close to the target biomass

          ## ...or can select dbh populations based on how close they are to median npp (not 100% overlapping)
          # npp.lims <- quantile(cage.ann.npp[[pix]], probs=npp.quant.range) ## restrict which of the simulations to draw and modify
          # j <- which(cage.ann.npp[[pix]]>=npp.lims[1] & cage.ann.npp[[pix]]<=npp.lims[2])
          # if(length(j)>=4){ ## if there's not enough dbh samples that meet your criteria...
          #   j <- sample(j,1) ## get a random dbh sample selected range of simulator results close to the target biomass
          # }else{j <- sample(1:length(cage.biom.sim[[pix]]), 1)} ## or just sample the whole simulation collection if all simulations are nearly identical

          #### load up a dbh sample and resimulate each tree for 36 consecutive years
          tree.samp <- cage.dbh[[procset[pix]]][[j]] ## what initial trees are present in this simulator result

          ## tree expansion scenario
          ### this is imperfect -- but assume a small amount of new street buffer planting space is available in this spot and sim new trees
          expand.track <- rbinom(n = 10, prob=expand.rate, size=1) ## give this tree a tree plant flag for the next ten years that on average across the map gets the right number of plantings in each pixel

          npp.track <- numeric() ## annual tally of npp in this resim as you go up in years
          biom.track <- numeric() ## annual tally of biomass in this resim
          num.track <- integer() ## annual tally of number of trees in this resim
          newplants <- 0 ## track number of new trees planted over course of 36 years
          deaths <- 0 ## keep track of how many times trees die in this pix resim
          can.track <- numeric() ### canopy coverage in this simulation
          ## determine initial replanting delays
          delay <- sample(seq(delay.factor[1],
                              delay.factor[2]),
                          size = length(tree.samp),
                          replace=T) 

          ### figure out biomass and NPP at start of simulation without any effects of growth/mort
          tmp.dbh0 <- tree.samp 
          tmp.genus <- cage.genus[[procset[pix]]][[j]]
          tmp.biom0 <- street.allo[match(tmp.genus, street.allo$genus, nomatch=8), "b0"]*(tmp.dbh0^street.allo[match(tmp.genus, street.allo$genus, nomatch=8), "b1"])*street.allo[match(tmp.genus, street.allo$genus, nomatch=8), "dens"]
          tmp.dbh1 <- tmp.dbh0+(b0.rand[a]+(b1.rand[a]*tmp.dbh0)+(b2.rand[a]*(tmp.dbh0^2))) ## grow the dbh to time 1, no growth in dead trees
          tmp.biom1 <- street.allo[match(tmp.genus, street.allo$genus, nomatch=8), "b0"]*(tmp.dbh1^street.allo[match(tmp.genus, street.allo$genus, nomatch=8), "b1"])*street.allo[match(tmp.genus, street.allo$genus, nomatch=8), "dens"]
          
          ### canopy: get vector of equation types to use first
          eq.form <- street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "EqName"]
          ## add up total canopy area applying correct equation form to each
          tmp.can0 <- sum(((((eq.form=="quad")*(street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "a"]+
                                                  (street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "b"]*tmp.dbh0)+
                                                  (street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "c"]*tmp.dbh0^2)))/2)^2)*pi, na.rm=T)+
            sum(((((eq.form=="cub")*(street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "a"]+
                                       (street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "b"]*tmp.dbh0)+
                                       (street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "c"]*tmp.dbh0^2)+
                                       (street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "d"]*tmp.dbh0^3)))/2)^2)*pi, na.rm=T)+
            sum(((((eq.form=="lin")*(street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "a"]+
                                       (street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "b"]*tmp.dbh0)))/2)^2)*pi, na.rm=T)+
            sum((((eq.form=="loglogw1")*exp(street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "a"]+
                                              (street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "b"]*
                                                 log(log(tmp.dbh0+1)+
                                                       street.canopy[match(tmp.genus, street.canopy$genus, nomatch=nogenus), "c"]/2)))/2)^2)*pi, na.rm=T)
          
          ## record start conditions
          npp.track <- c(npp.track, sum(tmp.biom1-tmp.biom0))
          num.track <- c(num.track, length(tmp.dbh0))
          biom.track <- c(biom.track, tmp.biom0)
          can.track <- c(can.track, tmp.can0)
          
          ## Now simulate 36 successive years of growth/mortality/replanting/plant expansion 2007-2043
          for(e in 1:33){
            if(scenario[s]=="expand" & e<=expand.timeline){
              tree.samp <- c(tree.samp, 5*expand.track[e]) ## add a 5cm tree if the expand.track wins one this year (within the first 10 years)
              newplants <- newplants+1
              delay <- c(delay, sample(seq(delay.factor[1],
                                           delay.factor[2]),
                                       size = 1,
                                       replace=T)) ## give this new tree its own delay clock
            }

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
          delay[kill.list==0 & tree.samp==0 & delay>0] <- delay[kill.list==0 & tree.samp==0 & delay>0]-1 ## count down the delay clocks that had been set before

          ## grow up the survivors using estimated growth regression
          # tree.samp <- tree.samp[tree.samp>0]*(1+dbhg.pred(tree.samp[tree.samp>0])) ## static realization of growth model
          tree.samp[tree.samp>0] <- tree.samp[tree.samp>0]+(b0.rand[a]+(b1.rand[a]*tree.samp[tree.samp>0])+(b2.rand[a]*tree.samp[tree.samp>0]^2))

          ### upgrade: urban specific allometries to determine biomass change
          tmp.dbh0 <- tree.samp[tree.samp>0] ## exclude dead trees from dbh/biomass change
          tmp.genus <- tmp.genus[tree.samp>0] ## genera of dbh collection

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
            replanted <- rbinom(length(tree.samp[tree.samp==0 & delay==0]), 1, replant.factor) ## decide which of the dead trees without a delay timer get a chance to be replaced
            tree.samp[tree.samp==0 & delay==0] <- replanted*5 ## replace selected dead trees that win the replacement lottery
          }
          ### the logic of this scheme is that delay indicates the number of years *after the death year* that there is 0 productivity
          ### trees lose all NPP in death year (growth calc happens before replant calc)
          ### delay range of 0-2 years implies a wait-out of 1-3 years before productivity happens again

          deaths <- deaths+sum(kill.list) ## count the dead
          num.track <- c(num.track, length(tree.samp[tree.samp>0])) ## count the living

        } ## end of loop for 36 year projector
          
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
        npp.box[[pix]][[a]] <- NA
        biom.box[[pix]][[a]] <- NA
        expand.plant.num[[pix]] <- NA
        deaths.sav[[pix]] <- NA
        num.box[[pix]][[a]] <- NA
        can.box[[pix]][[a]] <- NA
        print(paste("pix", index.track[procset[pix]], "failed, too few sims"))
      }

    } # end loop for this pixel   
    save(pixID.track, file=paste("processed/boston/biom_street/index.track.scenario.", scenario[s], ".", o, ".V", resim.vers, ".resim.sav", sep=""))
    save(deaths.sav, file=paste("processed/boston/biom_street/death.track.scenario.", scenario[s], ".", o, ".V", resim.vers, ".resim.sav", sep=""))
    save(npp.box, file=paste("processed/boston/biom_street/npp.track.scenario.", scenario[s], ".", o, ".V", resim.vers, ".resim.sav", sep=""))
    save(biom.box, file=paste("processed/boston/biom_street/biom.track.scenario.", scenario[s], ".", o, ".V", resim.vers, ".resim.sav", sep=""))
    save(expand.plant.num, file=paste("processed/boston/biom_street/expand.plant.num.scenario.", scenario[s], ".", o, ".V", resim.vers, ".resim.sav", sep=""))
    save(num.box, file=paste("processed/boston/biom_street/num.track.scenario.", scenario[s], ".", o, ".V", resim.vers, ".resim.sav", sep=""))
    save(can.box, file=paste("processed/boston/biom_street/can.track.scenario.", scenario[s], ".", o, ".V", resim.vers, ".resim.sav", sep=""))
    print("finished chunk")
  } ## loop for pixel chunk files
  print(paste("finished resim run scenario", scenario[s]))
}


###
### process resim results by scenario
#####
## get a general list of the resim results and the corresponding pixel index
sum.na <- function(x){sum(x, na.rm=T)}
scenario <- c("BAU", "oldies", "expand")
resim.vers <- 4
preamb <- "processed/boston/biom_street/"
## containers for the histories as they unfurl themselves in the resims

for(s in 1:length(scenario)){
  obj.list <- list.files("processed/boston/biom_street/")
  obj.list <- obj.list[grep(obj.list, pattern=paste("scenario", scenario[s], sep="."))]
  obj.list <- obj.list[grep(obj.list, pattern=paste0("V", resim.vers))]
  index.list <- obj.list[grep(obj.list, pattern="index")]
  chunks <- strsplit(index.list, split = "[.]")
  chunks <- sapply(chunks, "[[" ,5)
  npp.contain <- data.frame()
  biom.contain <- data.frame()
  can.contain <- data.frame()
  for(o in 1:length(chunks)){
    print(paste("processing chunk", chunks[o]))
    npp.mat <- matrix()
    can.mat <- matrix()
    # num.mat <- matrix()
    # death.mat <- matrix()
    biom.mat <- matrix()
    expand.mat <- matrix()
    load(paste0(preamb, "npp.track.scenario.", scenario[s], ".", chunks[o], ".V", resim.vers, ".resim.sav", sep="")) ## npp.box
    load(paste0(preamb, "biom.track.scenario.", scenario[s], ".", chunks[o], ".V", resim.vers, ".resim.sav", sep="")) ## biom.box
    load(paste0(preamb, "num.track.scenario.", scenario[s], ".", chunks[o], ".V", resim.vers, ".resim.sav", sep="")) ## num.box
    load(paste0(preamb, "can.track.scenario.", scenario[s], ".", chunks[o], ".V", resim.vers, ".resim.sav", sep="")) ## can.box
    load(paste0(preamb, "death.track.scenario.", scenario[s], ".", chunks[o], ".V", resim.vers, ".resim.sav", sep="")) ## deaths.sav
    load(paste0(preamb, "expand.plant.num.scenario.", scenario[s], ".", chunks[o], ".V", resim.vers, ".resim.sav", sep="")) ## expand.plant.num
    load(paste0(preamb, "index.track.scenario.", scenario[s], ".", chunks[o], ".V", resim.vers, ".resim.sav", sep="")) ## pixID.track
    
    kill <- sapply(npp.box, is.na) ## cancel the non-resimmed entreid
    npp.box[kill] <- NULL
    biom.box[kill] <- NULL
    num.box[kill] <- NULL
    can.box[kill] <- NULL
    expand.plant.num[kill] <- NULL
    pixID.track <- pixID.track[!kill]
    
    npp.mat <- matrix(unlist(npp.box), ncol=37, byrow = TRUE) ## this is year-by-column results for all pix in this chunk
    biom.mat <- matrix(unlist(biom.box), ncol=37, byrow = TRUE)
    # num.mat ## this is missing a tooth
    can.mat <- matrix(unlist(can.box), ncol=37, byrow = TRUE)
    
    npp.mat <- cbind(pixID.track, npp.mat)
    biom.mat <- cbind(pixID.track, biom.mat)
    can.mat <- cbind(pixID.track, can.mat)
    # bu <- apply(npp.mat[,2:38], MARGIN = 2, FUN=sum.na)
    # plot(1:37, bu/2000)
    npp.contain <- rbind(npp.contain, npp.mat)
    biom.contain <- rbind(biom.contain, biom.mat)
    can.contain <- rbind(can.contain, can.mat)
  }
  print(paste("writing collated", scenario[s], "to disk"))
  write.csv(npp.contain, paste0("processed/results/", scenario[s], ".V", resim.vers, ".npp.trendmap.csv"))
  write.csv(biom.contain, paste0("processed/results/", scenario[s], ".V", resim.vers, ".biom.trendmap.csv"))
  write.csv(can.contain, paste0("processed/results/", scenario[s], ".V", resim.vers, ".can.trendmap.csv"))
}
#####


## Exploratory of resim results
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
nonfor <- biom.dat[bos.aoi30m>800 & bos.lulc30m.lumped!=1 & bos.biom30m<20000, pix.ID] ### identify pixel ID's that are nonforest
# hyb <- as.data.table(read.csv("processed/results/hybrid.results.V6.csv"))
hyb.t <- as.data.table(as.data.frame(raster("processed/results/hybrid.V6.median.tif")))
hyb.t[,pix.ID:=seq(1, dim(hyb.t)[1])]
biom.dat <- merge(biom.dat, hyb.t, by="pix.ID")
dim(biom.dat[bos.aoi30m>800,]) ## 136667 pix in the AOI

## BAU scenario
biom.dat[bos.aoi30m>800, sum(hybrid.V6.median, na.rm=T)/2000/1000] ### 10.8 ktC in ~2007
bau.npp <- read.csv("processed/results/BAU.V3.npp.trendmap.csv")
bau.can <- read.csv("processed/results/BAU.V3.can.trendmap.csv")
bau.biom <- read.csv("processed/results/BAU.V3.biom.trendmap.csv")

## BAU NPP start/finish
biom.dat <- merge(biom.dat, bau.npp[,c("pixID.track", "V2", "V38")], by.x="pix.ID", by.y="pixID.track", all.x=T)
dim(biom.dat[bos.aoi30m>800,]) ## still 136667
names(biom.dat)[8:9] <- c("BAU.start.npp", "BAU.finish.npp")
dim(biom.dat[bos.aoi30m>800 & !(is.na(BAU.start.npp))]) ## 92076
plot(biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), hybrid.V6.median], 
     biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), BAU.start.npp])
abline(a=0, b=1, col="red")
plot(biom.dat[bos.aoi30m>800 & !is.na(BAU.finish.npp), hybrid.V6.median], 
     biom.dat[bos.aoi30m>800 & !is.na(BAU.finish.npp), BAU.finish.npp])
abline(a=0, b=1, col="blue")
biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), sum(hybrid.V6.median, na.rm=T)/2000/1000] ## 8.9 ktC IN THE RESIMMED AREA
biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), sum(BAU.start.npp, na.rm=T)/2000/1000] ## 8.3 ktC
biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), sum(BAU.finish.npp, na.rm=T)/2000/1000] ## 9.0 ktC  ### 8% change over time

## BAU canopy start/finish
biom.dat <- merge(biom.dat, bau.can[,c("pixID.track", "V2", "V38")], by.x="pix.ID", by.y="pixID.track", all.x=T)
names(biom.dat)[10:11] <- c("BAU.start.can", "BAU.finish.can")
plot(biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), bos.can30m*bos.aoi30m], 
     biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), BAU.start.can])
abline(a=0, b=1, col="red")
## not even close, way overpredicted by the model
plot(biom.dat[bos.aoi30m>800 & !is.na(BAU.finish.npp), bos.can30m*bos.aoi30m], 
     biom.dat[bos.aoi30m>800 & !is.na(BAU.finish.npp), BAU.finish.can])
abline(a=0, b=1, col="blue")
biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), sum(bos.can30m*bos.aoi30m, na.rm=T)/1E4] ## 2.8 kha in the resimmed area
biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), sum(bos.aoi30m, na.rm=T)/1E4] ## 8.3 kha total resimmed area = 33.7% canopy
biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), sum(BAU.start.can, na.rm=T)/1E4] ## 4.9 kha start in the resimmed area
biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), sum(BAU.finish.can, na.rm=T)/1E4] ## 4.9 kha finish in the resimmed area, slight decline


biom.dat[bos.aoi30m>800, sum(bos.can30m*bos.aoi30m, na.rm=T)]/biom.dat[bos.aoi30m>800, sum(bos.aoi30m)] ## 31.8% canopy in the AOI
biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), length(bos.aoi30m)] ## 92076 pix resimmed
biom.dat[bos.aoi30m>800, length(bos.aoi30m)] ## 136667 pix with data

## BAU biomass start/finish
biom.dat <- merge(biom.dat, bau.biom[,c("pixID.track", "V2", "V38")], by.x="pix.ID", by.y="pixID.track", all.x=T)
names(biom.dat)[12:13] <- c("BAU.start.biom", "BAU.finish.biom")
plot(biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), bos.biom30m], 
     biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), BAU.start.biom])
abline(a=0, b=1, col="red") ## LIKE A FUCKING GLOVE
plot(biom.dat[bos.aoi30m>800 & !is.na(BAU.finish.npp), bos.biom30m], 
     biom.dat[bos.aoi30m>800 & !is.na(BAU.finish.npp), BAU.finish.biom])
abline(a=0, b=1, col="blue")
biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), sum(bos.biom30m, na.rm=T)/2000/1000] ## 209 ktC in the resimmed area compare 357 ktC in AOI
biom.dat[bos.aoi30m>800, sum(BAU.finish.biom, na.rm=T)/2000/1000] ## 266 ktC by the end of simulation, 27% increase

## time series of BAU 
bau.npp.t <- apply(bau.npp[,3:39], MARGIN = 2, FUN = sum.na)
plot(1:37, bau.npp.t/2000)
last(bau.npp.t)/2000; first(bau.npp.t)/2000

bau.can.t <- apply(bau.can[,3:39], MARGIN = 2, FUN = sum.na)
plot(1:37, bau.can.t/2000)
last(bau.can.t)/1E4; first(bau.can.t)/1E4

bau.biom.t <- apply(bau.biom[,3:39], MARGIN = 2, FUN = sum.na)
plot(1:37, bau.biom.t/2000/1000)
last(bau.biom.t)/2000/1000; first(bau.biom.t)/2000/1000


### Oldies scenario
biom.dat[bos.aoi30m>800, sum(hybrid.V6.median, na.rm=T)/2000/1000] ### 10.8 ktC in ~2007
old.npp <- read.csv("processed/results/oldies.V3.npp.trendmap.csv")
old.can <- read.csv("processed/results/oldies.V3.can.trendmap.csv")
old.biom <- read.csv("processed/results/oldies.V3.biom.trendmap.csv")

## oldies NPP start/finish
biom.dat <- merge(biom.dat, old.npp[,c("pixID.track", "V2", "V38")], by.x="pix.ID", by.y="pixID.track", all.x=T)
dim(biom.dat[bos.aoi30m>800,]) ## still 136667
names(biom.dat)[14:15] <- c("old.start.npp", "old.finish.npp")
dim(biom.dat[bos.aoi30m>800 & !(is.na(old.start.npp))]) ## 92076
plot(biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), hybrid.V6.median], 
     biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), old.start.npp])
abline(a=0, b=1, col="red")
plot(biom.dat[bos.aoi30m>800 & !is.na(old.finish.npp), hybrid.V6.median], 
     biom.dat[bos.aoi30m>800 & !is.na(old.finish.npp), old.finish.npp])
abline(a=0, b=1, col="blue")
biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), sum(hybrid.V6.median, na.rm=T)/2000/1000] ## 8.9 ktC IN THE RESIMMED AREA
biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), sum(old.start.npp, na.rm=T)/2000/1000] ## 8.3 ktC
biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), sum(old.finish.npp, na.rm=T)/2000/1000] ## 10.2 ktC  ### 23% change over time

## oldies canopy start/finish
biom.dat <- merge(biom.dat, old.can[,c("pixID.track", "V2", "V38")], by.x="pix.ID", by.y="pixID.track", all.x=T)
names(biom.dat)[16:17] <- c("old.start.can", "old.finish.can")
plot(biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), bos.can30m*bos.aoi30m], 
     biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), old.start.can])
abline(a=0, b=1, col="red")
## not even close, way overpredicted by the model
plot(biom.dat[bos.aoi30m>800 & !is.na(old.finish.npp), bos.can30m*bos.aoi30m], 
     biom.dat[bos.aoi30m>800 & !is.na(old.finish.npp), old.finish.can])
abline(a=0, b=1, col="blue")
biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), sum(bos.can30m*bos.aoi30m, na.rm=T)/1E4] ## 2.8 kha in the resimmed area
biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), sum(bos.aoi30m, na.rm=T)/1E4] ## 8.3 kha total resimmed area = 33.7% canopy
biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), sum(old.start.can, na.rm=T)/1E4] ## 4.9 kha start in the resimmed area
biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), sum(old.finish.can, na.rm=T)/1E4] ## 6.1 kha finish in the resimmed area, slight decline

## oldies biomass start/finish
biom.dat <- merge(biom.dat, old.biom[,c("pixID.track", "V2", "V38")], by.x="pix.ID", by.y="pixID.track", all.x=T)
names(biom.dat)[18:19] <- c("old.start.biom", "old.finish.biom")
plot(biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), bos.biom30m], 
     biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), old.start.biom])
abline(a=0, b=1, col="red") ## LIKE A FUCKING GLOVE
plot(biom.dat[bos.aoi30m>800 & !is.na(old.finish.npp), bos.biom30m], 
     biom.dat[bos.aoi30m>800 & !is.na(old.finish.npp), old.finish.biom])
abline(a=0, b=1, col="blue")
biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), sum(bos.biom30m, na.rm=T)/2000/1000] ## 209 ktC in the resimmed area compare 357 ktC in AOI
biom.dat[bos.aoi30m>800, sum(old.finish.biom, na.rm=T)/2000/1000] ## 361 ktC by the end of simulation, 76% increase

## time series of oldies 
old.npp.t <- apply(old.npp[,3:39], MARGIN = 2, FUN = sum.na)
plot(1:37, old.npp.t/2000)
last(old.npp.t)/2000; first(old.npp.t)/2000

old.can.t <- apply(old.can[,3:39], MARGIN = 2, FUN = sum.na)
plot(1:37, old.can.t/2000)
last(old.can.t)/1E4; first(old.can.t)/1E4

old.biom.t <- apply(old.biom[,3:39], MARGIN = 2, FUN = sum.na)
plot(1:37, old.biom.t/2000/1000)
last(old.biom.t)/2000/1000; first(old.biom.t)/2000/1000


### Expand scenario
biom.dat[bos.aoi30m>800, sum(hybrid.V6.median, na.rm=T)/2000/1000] ### 10.8 ktC in ~2007
exp.npp <- read.csv("processed/results/expand.V3.npp.trendmap.csv")
exp.can <- read.csv("processed/results/expand.V3.can.trendmap.csv")
exp.biom <- read.csv("processed/results/expand.V3.biom.trendmap.csv")

## oldies NPP start/finish
biom.dat <- merge(biom.dat, exp.npp[,c("pixID.track", "V2", "V38")], by.x="pix.ID", by.y="pixID.track", all.x=T)
dim(biom.dat[bos.aoi30m>800,]) ## still 136667
names(biom.dat)[20:21] <- c("exp.start.npp", "exp.finish.npp")
dim(biom.dat[bos.aoi30m>800 & !(is.na(exp.start.npp))]) ## 92076
plot(biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), hybrid.V6.median], 
     biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), exp.start.npp])
abline(a=0, b=1, col="red")
plot(biom.dat[bos.aoi30m>800 & !is.na(exp.finish.npp), hybrid.V6.median], 
     biom.dat[bos.aoi30m>800 & !is.na(exp.finish.npp), exp.finish.npp])
abline(a=0, b=1, col="blue")
biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), sum(hybrid.V6.median, na.rm=T)/2000/1000] ## 8.9 ktC IN THE RESIMMED AREA
biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), sum(exp.start.npp, na.rm=T)/2000/1000] ## 8.3 ktC
biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), sum(exp.finish.npp, na.rm=T)/2000/1000] ## 9.5 ktC  ### 14% change over time

## expand canopy start/finish
biom.dat <- merge(biom.dat, exp.can[,c("pixID.track", "V2", "V38")], by.x="pix.ID", by.y="pixID.track", all.x=T)
names(biom.dat)[22:23] <- c("exp.start.can", "exp.finish.can")
plot(biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), bos.can30m*bos.aoi30m], 
     biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), exp.start.can])
abline(a=0, b=1, col="red")
## not even close, way overpredicted by the model
plot(biom.dat[bos.aoi30m>800 & !is.na(exp.finish.npp), bos.can30m*bos.aoi30m], 
     biom.dat[bos.aoi30m>800 & !is.na(exp.finish.npp), exp.finish.can])
abline(a=0, b=1, col="blue")
biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), sum(bos.can30m*bos.aoi30m, na.rm=T)/1E4] ## 2.8 kha in the resimmed area
biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), sum(bos.aoi30m, na.rm=T)/1E4] ## 8.3 kha total resimmed area = 33.7% canopy
biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), sum(exp.start.can, na.rm=T)/1E4] ## 4.9 kha start in the resimmed area
biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), sum(exp.finish.can, na.rm=T)/1E4] ## 5.1 kha finish in the resimmed area, slight decline

## Expand biomass start/finish
biom.dat <- merge(biom.dat, exp.biom[,c("pixID.track", "V2", "V38")], by.x="pix.ID", by.y="pixID.track", all.x=T)
names(biom.dat)[24:25] <- c("exp.start.biom", "exp.finish.biom")
plot(biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), bos.biom30m], 
     biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), exp.start.biom])
abline(a=0, b=1, col="red") ## LIKE A FUCKING GLOVE
plot(biom.dat[bos.aoi30m>800 & !is.na(exp.finish.npp), bos.biom30m], 
     biom.dat[bos.aoi30m>800 & !is.na(exp.finish.npp), exp.finish.biom])
abline(a=0, b=1, col="blue")
biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), sum(bos.biom30m, na.rm=T)/2000/1000] ## 209 ktC in the resimmed area compare 357 ktC in AOI
biom.dat[bos.aoi30m>800, sum(exp.finish.biom, na.rm=T)/2000/1000] ## 275 ktC by the end of simulation, 31% increase

## time series of expand 
exp.npp.t <- apply(exp.npp[,3:39], MARGIN = 2, FUN = sum.na)
plot(1:37, exp.npp.t/2000)
last(exp.npp.t)/2000; first(exp.npp.t)/2000

exp.can.t <- apply(exp.can[,3:39], MARGIN = 2, FUN = sum.na)
plot(1:37, exp.can.t/2000)
last(exp.can.t)/1E4; first(exp.can.t)/1E4

exp.biom.t <- apply(exp.biom[,3:39], MARGIN = 2, FUN = sum.na)
plot(1:37, exp.biom.t/2000/1000)
last(exp.biom.t)/2000/1000; first(exp.biom.t)/2000/1000












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
