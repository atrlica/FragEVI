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

### apply urban-specific allometries
street.allo <- read.csv("docs/street.biometrics.csv") ## AG wood vol (m3) as b0*DBH(cm)^b1; multiply by density to get kg-biomass

### use urban specific crown diameter measures to get canopy coverage
# street.canopy <- as.data.table(read.csv("docs/RDS-2016-0005/Data/TS6_Growth_coefficients.csv"))
# street.canopy <- street.canopy[Region=="NoEast", ]
# street.canopy <- street.canopy[Predicts.component=="crown dia",]
# genspec <- strsplit(as.character(street.canopy$Scientific.Name), " ")
# gen <- unlist(lapply(genspec, "[[", 1))
# street.canopy[,genus:=gen] ## 40 genera
# keepers <- c("Acer platanoides", "Aesculus hippocastanum", "Fraxinus pennsylvanica",
#              "Ginkgo biloba", "Gleditsia triacanthos", "Liquidambar styraciflua", 
#              "Malus spp.", "Platanus x acerifolia", "Prunus serrulata", "Pyrus calleryana",
#              "Quercus rubra", "Tilia cordata", "Ulmus americana", "Zelkova serrata")
# street.canopy <- street.canopy[Scientific.Name%in%keepers,]
# street.canopy$a <- as.numeric(as.character(street.canopy$a))
# street.canopy$b <- as.numeric(as.character(street.canopy$b))
# street.canopy$c <- as.numeric(as.character(street.canopy$c))
# street.canopy$d <- as.numeric(as.character(street.canopy$d))
# street.canopy$e <- as.numeric(as.character(street.canopy$e))
# street.canopy <- as.data.frame(street.canopy)
nogenus <- 1 ## resolved: Acer platanoides is the most common spp -- for the handful of spp not in this list, we will use the A. platanoides canopy equation
# write.csv(street.canopy, "docs/street.canopy.csv")
street.canopy <- read.csv("docs/street.canopy.csv")

## load up street tree survey data for use in pieces of this thing
street <- as.data.table(read.csv("processed/boston/street.trees.dbh.csv"))

###  get a list of model coefficients to run for each repetition (dbh.delta~dbh.start)
# library(lme4)  ### won't work on cluster, have to hard code shit because R version on cluster is old
# load("processed/mod.street.dbhdelta.me.sav")
# aaa <- summary(mod.street.dbhdelta.me)
# b0.rand <- rnorm(100, coef(aaa)[1,1], coef(aaa)[1,2])
# b1.rand <- rnorm(100, coef(aaa)[2,1], coef(aaa)[2,2])
# b2.rand <- rnorm(100, coef(aaa)[3,1], coef(aaa)[3,2])
# b0.hard <- c(1.234, 7.905e-02) ### these are the sampling distributions of the dbh growth predictor coefficients
# b1.hard <- c(-1.982e-02, 3.698e-03)
# b2.hard <- c(1.329e-04, 3.369e-05)
# b0.rand <- rnorm(100, b0.hard[1], b0.hard[2])
# b1.rand <- rnorm(100, b1.hard[1], b1.hard[2])
# b2.rand <- rnorm(100, b2.hard[1], b2.hard[2])
# write.csv(b0.rand, "processed/boston/biom_street/b0.rand.csv") ## export a stable collection of coefficients
# write.csv(b1.rand, "processed/boston/biom_street/b1.rand.csv")
# write.csv(b2.rand, "processed/boston/biom_street/b2.rand.csv")
# b0.rand <- b0.hard[1]
# b1.rand <- b1.hard[1]
# b2.rand <- b2.hard[1]
b0.rand <- read.csv("processed/boston/biom_street/b0.rand.csv")
b0.rand <- as.vector(b0.rand[,2])
b1.rand <- read.csv("processed/boston/biom_street/b1.rand.csv")
b1.rand <- as.vector(b1.rand[,2])
b2.rand <- read.csv("processed/boston/biom_street/b2.rand.csv")
b2.rand <- as.vector(b2.rand[,2])

###
### expand scenario data file processing
#####
# library(rgeos)
# library(rgdal)
# library(raster)
# rbuff <- readOGR("processed/boston/newplanting/road2345_buffsD.shp") ## this is the prepared road buffers polygon
# isa <- raster("processed/boston/bos.isa.RR2.tif")
# aoi <- raster("processed/boston/bos.aoi.tif")
# # rbuff.t <- spTransform(rbuff, CRSobj = crs(isa))
# # writeOGR(rbuff.t, layer="road2345_buffsUTM", dsn = "processed/boston/newplanting/road2345_buffsUTM.shp", driver = "ESRI Shapefile") ## this gets rasterized to the isa.RR2 grid in arc
# rbuff.r <- raster("processed/boston/newplanting/road2345_buffsR.tif") ## this is the 1m raster of the road buffers, on the isa grid
# rbuff.r <- crop(rbuff.r, aoi)
# rbuff.r <- extend(rbuff.r, aoi)
# isa <- crop(isa, aoi)
# #### this function ID's road buffer pixels that are also pervious
# id.perv <- function(buff, perv, aoi, filename) {
#   out <- raster(buff)
#   bs <- blockSize(out)
#   out <- writeStart(out, filename, overwrite=TRUE, format="GTiff")
#   for (i in 1:bs$n) {
#     b <- getValues(buff, row=bs$row[i], nrows=bs$nrows[i])
#     p <- getValues(perv, row=bs$row[i], nrows=bs$nrows[i])
#     a <- getValues(aoi, row=bs$row[i], nrows=bs$nrows[i])
#     g <- rep(0, length(b))
#     g[b==0 & p==0 & a==1] <- 1 ## mark all inside-aoi non-pervious buffer pix as good
#     g[is.na(a)] <- NA
#     g[b!=0] <- NA ## NA out non buffers and buffers over impervious
#     out <- writeValues(out, g, bs$row[i])
#     print(paste("finished block", i, "of", bs$n))
#   }
#   out <- writeStop(out)
#   return(out)
# }
# s <- id.perv(rbuff.r, isa, aoi, filename="processed/boston/newplanting/plantable.tif")
# 
# holes <- read.csv("processed/boston/newplanting/road2345_donut_line_rec.csv")
# holes$plantable.length <- holes$Shape_Length/2
# holes$num.trees <- floor(holes$plantable.length/8)+1
# hist(holes$num.trees); boxplot(holes$num.trees, outline=F) ## bulk of them are 8 and below
# sum(holes$num.trees) ## 170147 tree spaces plantable
# sum(holes$plantable.length)/1000 ## 1197 km of plantable road buffer


### this was the old way of doing the processing to determine street tree plantable area
######
# box <- as.data.table(read.csv("processed/boston/plant10m_MBG.csv")); dim(box) ## 78618
# summary(box$Shape_Area); hist(box$Shape_Area)
# box <- box[Shape_Area<1000,]; dim(box) ##78034 ## eliminate a fair chunk of the weird giant boxes
# box[,num.trees:=((MBG_Width-1)%/%8)+1] ## this gives us 2 trees per 8m planter length with a little buffer at either end
# box[,sum(num.trees)] ## 79k trees can be planted in this population of boxes
######
#####

## load up ancillary map data
library(raster)
library(rgdal)
library(data.table)
biom <- raster("processed/boston/bos.biom30m.tif")
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- crop(biom, aoi)
can <- raster("processed/boston/bos.can.redux30m.tif")
isa <- raster("processed/boston/bos.isa30m.tif")
biom.dat <- as.data.table(as.data.frame(biom))
aoi.dat <- as.data.table(as.data.frame(aoi))
can.dat <- as.data.table(as.data.frame(can))
isa.dat <- as.data.table(as.data.frame(isa))
lulc <- raster("processed/boston/bos.lulc30m.lumped.tif")
lulc.dat <- as.data.table(as.data.frame(lulc))
biom.dat <- cbind(biom.dat, aoi.dat, can.dat, isa.dat, lulc.dat)
biom.dat[,pix.ID:=seq(1:dim(biom.dat)[1])]
nonfor <- biom.dat[bos.aoi30m>800 & !(bos.lulc30m.lumped %in% c(1,4,5,6)) & bos.biom30m<20000, pix.ID] ### identify pixel ID's that are dev, hdres, ldres

### actual simulator machine here
######
### set scenario options 
# scenario <- c("BAU", "highmort", "lowreplant", "oldies", "lowmort", "slowreplant", "expand")
scenario <- c("BAU", "oldies", "expand")
# scenario <- "expand"
## record specific parameter sets
resim.vers <- 8 ## what are we labeling this round of resims?
vers <- 6 ## what simulator results version are we dealing with?

### default factors
npp.quant.range <- c(0.025, 0.975) ## what dbh samples to draw from
mort.mod <- 1 ## default modification to mortality rate (scenario specific)
default.sizecutoff <- 10000 ## size of tree dbh over which to screw with mortality
largemort.mod <- 1 ## multiplier for mortalities in large trees
replant.factor <- 1 ## rate of replanting
dbh.big <- default.sizecutoff
default.delay <- c(0,2) ## low/hi on years of delay from death to replanting
delay.factor <- default.delay
yr.run <- 34 ## how many years past year 0 to run the scenario
realize <- 20 ## how many realizations to put every pixel resim through (model+simulator error)

## set options for scenario processing
highmort.mortfactor <- 1.25
lowmort.mortfactor <- 0.5
lowreplant.factor <- 0.5 ## what fraction of morts are replanted in "lowreplant"
oldies.mortfactor <- 0.5 ## how much to reduce mortalities in large trees
oldies.sizecutoff <- 40 ## how to define "big" trees
slowreplant.delay <- c(0,2) ## 1 to 3 year delay between death and regrowth starting
expand.timeline <- 10 ### how many years to implement the expand scenario
# expand.rate.go <- (sum(holes$num.trees)/expand.timeline)/length(nonfor) ## new buffer calc, this is the per-year rate of new stem appearance in NONFOREST pixels needed to get us to full expand coverage after number of years = expand.timeline
expand.rate.go <- (170147/expand.timeline)/length(nonfor) ## hard code the fucker
# expand.rate.go*length(nonfor)*10 ## i.e. planting expand rate/yr in each pixel over ten years gives you the target tree number
# yr1 <- rbinom(n=length(nonfor), size=1, prob = expand.rate.go)
# (sum(yr1)*10)/170147 ## ok this should get us about what we need, close to the total required

# record run parameters to text file
# params.list <- list(c(scenario, resim.vers, vers,
#                     highmort.mortfactor, lowmort.mortfactor, 
#                     , lowreplant.factor, oldies.sizecutoff, 
#                     oldies.mortfactor))
sink(paste0("processed/boston/biom_street/resim_params_list_V", resim.vers, ".txt"))
cat(c("Resimulation version: ", resim.vers), "\n")
cat(c("Simulator results version: ", vers), "\n")
cat(c("scenarios = ", scenario), "\n")
cat(c("pixels resimmed = ", length(nonfor)), "\n")
cat(c("resims per pixel = ", realize), "\n")
cat(c("years resimulated = ", yr.run), "\n")
cat(c("NPP selection quantiles = ", npp.quant.range), "\n")
cat(c("default replant delay yrs = ", delay.factor), "\n")
cat(c("highmort mortality factor = ", highmort.mortfactor), "\n")
cat(c("lowmort mortality factor = ", lowmort.mortfactor), "\n")
cat(c("oldies mortality factor = ", oldies.mortfactor), "\n")
cat(c("oldies size cutoff = ", oldies.sizecutoff), "\n")
cat(c("lowreplant replanting rate = ", lowreplant.factor), "\n")
cat(c("slowreplant delay time yrs = ", slowreplant.delay), "\n")
cat(c("expand per-pixel prob of tree increase/yr = ", expand.rate.go), "\n")
cat(c("expand timeline yrs =", expand.timeline), "\n")
sink()

## loops for simulating growth+mortality in pixel dbh samples
## pull in all the pixel simulator samples of tree dbh in every pixel
obj.list <- list.files("processed/boston/biom_street/")
obj.list <- obj.list[grep(obj.list, pattern=paste("dbh.street.v", vers, sep=""))]
obj.list <- sub('.*weighted\\.', '', obj.list)
obj.list <- (sub('\\..*', '', obj.list))

## check existing resim files, find next chunk to write
check <- list.files("processed/boston/biom_street")
check <- check[grep(check, pattern=paste("index.track.scenario"))] ### version label here
check <- check[grep(check, pattern=paste0("V", resim.vers))]
library(stringr)
already <- str_match(check, paste0("BAU.(.*?).V", resim.vers))
if(length(already)!=0){
  already <- already[,2]
  already <- already[!is.na(already)]} else{already <- integer()}

notyet <- obj.list[!(obj.list%in%already)]

## if any chunks are not fully processed, next check to see if they're being worked on currently
## the one weakness of this (besides requiring manual launch of each chunk)
## is that if a script times out or otherwise aborts before successfully completing after the step below
## it can't release the unfinished job from being "in process" and will keep the other scripts from running that chunk
if(length(notyet)!=0){ 
  ## make an empty record file if no file exists
  if(!file.exists(paste("processed/boston/biom_street/resim.atwork", resim.vers, "csv", sep="."))){ ## if it isn't there start a file of who is working now
    l <- data.frame(at.work=integer())
    write.csv(l, file=paste("processed/boston/biom_street/resim.atwork", resim.vers, "csv", sep="."))
  }
  ## read the at work file and see which jobs are running
  atwork <- read.csv(paste("processed/boston/biom_street/resim.atwork", resim.vers, "csv", sep="."))
  atwork <- atwork$at.work
  ## set target for what isn't completed and isn't being worked on
  y <- as.numeric(min(notyet[!(notyet%in%atwork)])) 
  
  if(length(y)==0 | !is.finite(y)){stop("all resim pixels currently finished or in process")}
  print(paste("going to work on chunk", y)) 
  
  ## update the at.work file to warn other job instances
  atwork <- c(atwork, y)
  l <- data.frame(at.work=atwork)
  write.csv(l, file=paste("processed/boston/biom_street/resim.atwork", resim.vers, "csv", sep="."))
  
}else{stop("all pixels already processed")}

## upgrade -- parallelize
## just assign static value for chunk to hitlist and feed in
hitlist <- y

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
  print(paste("starting resim run scenario", scenario[s], "with", realize, "realizations per pixel"))
  
  # ### avoid processing shit that's already done
  # already <- list.files("processed/boston/biom_street")
  # already <- already[grep(already, pattern=".resim.")]
  # already <- already[grep(already, pattern=paste0(".V", resim.vers))]
  # already <- already[grep(already, pattern="npp.track.")]
  # already <- already[grep(already, pattern=paste(".", scenario[s], ".", sep=""))]
  # already <- strsplit(already, split = "[.]")
  # already <- sapply(already, "[[" ,5)
  
# #   hitlist <- as.character(obj.list)
#   hitlist <- as.character(obj.list[!(obj.list %in%  already)])
#   hitlist[hitlist=="100000"] <- "1e+05"
#   if(length(hitlist)==0){print(paste(scenario[s], "already processed"))}
  
  if(length(hitlist)>0){
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
      
      for(pix in 1:length(procset)){ ## the number of pixels in this sim chunk
        # print(paste("working on pixel", index.track[procset[pix]]))
        pixID.track <- c(pixID.track, index.track[procset[pix]])
        deaths.sav[[pix]] <- integer()
        npp.box[[pix]] <- list()
        biom.box[[pix]] <- list()
        expand.plant.num[[pix]] <- integer()
        num.box[[pix]] <- list()
        can.box[[pix]] <- list()

          for(a in 1:realize){   ### resimulate every pixel realize number of times; about 13 hours per scenario at 10x resim realizations
            if(length(cage.biom.sim[[procset[pix]]])>=40){ ## only process if enough simulations successfully completed in this pixel
            # print(paste("resimming pix", index.track[pix], "realization", a))
            ### first select which pixel simulation you're drawing from in this pixel*realization
            # can select dbh populations based on proximity to median simulated biomass...
            biom.lims <- quantile(cage.biom.sim[[procset[pix]]], probs=npp.quant.range) ## figure out which of the simulations to draw and modify
            j <- sample(which(cage.biom.sim[[procset[pix]]]>=biom.lims[1] & cage.biom.sim[[procset[pix]]]<=biom.lims[2]),1) ## get a random dbh sample from the the middle 10% of simulator results close to the target biomass
            ## just in case...
            if(length(j)==0){j <- sample(1:length(cage.biom.sim[[procset[pix]]]), 1)} ## if you come up empty just sample any fucking thing
            
              ## this sampling method won't work for later pixel sims because we did not calculate npp in the simulator
  #           # ...or can select dbh populations based on how close they are to median npp (not 100% overlapping)
  #           npp.lims <- quantile(cage.ann.npp[[pix]], probs=npp.quant.range) ## restrict which of the simulations to draw and modify
  #           j <- which(cage.ann.npp[[pix]]>=npp.lims[1] & cage.ann.npp[[pix]]<=npp.lims[2])
  #           if(length(j)>=4){ ## if there's not enough dbh samples that meet your criteria...
  #             j <- sample(j,1) ## get a random dbh sample selected range of simulator results close to the target biomass
  #           }else{j <- sample(1:length(cage.biom.sim[[pix]]), 1)} ## or just sample the whole simulation collection if all simulations are nearly identical
  
            #### load up this dbh sample and resimulate each tree for yr.run consecutive years
            tree.samp <- cage.dbh[[procset[pix]]][[j]] ## what initial trees are present in this simulator result
  
            ## initialize trackers for this resim
            expand.track <- rbinom(n = expand.timeline, prob=expand.rate, size=1) ## give this tree an annual tree plant flag for the next ten years that on average across the map gets the right number of plantings in each pixel
            # expand.track=c(rep(1, 2),rep(0, 10))
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
            tmp.genus <- as.character(cage.genus[[procset[pix]]][[j]])
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
            
            ## record start conditions (year 0)
            npp.track <- c(npp.track, sum(tmp.biom1-tmp.biom0))
            num.track <- c(num.track, length(tmp.dbh0))
            biom.track <- c(biom.track, sum(tmp.biom0))
            can.track <- c(can.track, tmp.can0)
            
            ## Now simulate yr.run successive years of growth/mortality/replanting/plant expansion 2007-20???
            for(e in 1:yr.run){ ## yr.run=34 gets you from 2006 to 2040
              # print(paste("begin year", e))
              if(scenario[s]=="expand" & e<=expand.timeline & expand.track[e]==1){
                tree.samp <- c(tree.samp, 5) ## add a 5cm tree if the expand.track wins one this year (within the first 10 years)
                newplants <- newplants+1
                delay <- c(delay, sample(seq(delay.factor[1],
                                             delay.factor[2]),
                                         size = 1,
                                         replace=T)) ## give this new tree its own delay clock
                tmp.genus <- c(tmp.genus, as.character(sample(street[record.good==1, genus], 1)))
#                 print(paste("added expand tree of genus", tail(tmp.genus, n=1)))
              }
  
            ### figure mortality and determine a kill list
            deathwatch <- mort(tree.samp)*mort.mod ## mortality liklihood for each stem
            deathwatch[tree.samp>=dbh.big] <- deathwatch[tree.samp>=dbh.big]*largemort.mod ## adjust mortality in larger trees
            kill.list <- rbinom(length(tree.samp), 1, deathwatch/100) ## randomly kill based on mortality rate
            # kill.list <- rep(1,length(tree.samp)) ## test: kill everything
            kill.list[tree.samp==0] <- 0 ## do not rekill the previously dead
            tree.samp[kill.list==1] <- 0 ## the dead are nullified
#             if(sum(kill.list>0)){
#               print(paste(sum(kill.list), "tree(s) of genus", tmp.genus[kill.list==1], "have died"))
#             }
  
            ## determine a delay clock for the stems killed this cycle
            delay[kill.list==1] <- sample(seq(delay.factor[1],
                                              delay.factor[2]),
                                          size = length(delay[kill.list==1]),
                                          replace=T) # start a clock for any tree killed this year
  
            ## figure biomass and forward project npp in survivors using estimated growth regression
            tmp.dbh0.live <- tree.samp[tree.samp>0] ## exclude dead trees from dbh/biomass change
            tmp.genus.live <- tmp.genus[tree.samp>0] ## genera of live dbh collection
            
            tmp.dbh0.live <- tmp.dbh0.live+(b0.rand[a]+(b1.rand[a]*tmp.dbh0.live)+(b2.rand[a]*tmp.dbh0.live^2))
            tmp.biom0.live <- street.allo[match(tmp.genus.live, street.allo$genus, nomatch=8), "b0"]*(tmp.dbh0.live^street.allo[match(tmp.genus.live, street.allo$genus, nomatch=8), "b1"])*street.allo[match(tmp.genus.live, street.allo$genus, nomatch=8), "dens"]
            tree.samp[tree.samp>0] <- tmp.dbh0.live ## update the master list for trees that are alive and growing
            
            ## forward project next year's npp with this year's biomass 
            tmp.dbh1 <- tmp.dbh0.live+(b0.rand[a]+(b1.rand[a]*tmp.dbh0.live)+(b2.rand[a]*(tmp.dbh0.live^2))) ## grow the dbh to time 1, no growth in dead trees
            tmp.biom1 <- street.allo[match(tmp.genus.live, street.allo$genus, nomatch=8), "b0"]*(tmp.dbh1^street.allo[match(tmp.genus.live, street.allo$genus, nomatch=8), "b1"])*street.allo[match(tmp.genus.live, street.allo$genus, nomatch=8), "dens"]
            if(length(tmp.biom0.live)==0){tmp.biom0.live <- 0; tmp.biom1 <- 0} ## if everything is dead
  
            ### Determine canopy cover in THIS YEAR based on live stem dbh and matching street tree allometrics by genus
            ### get vector of equation types to use first
            eq.form <- street.canopy[match(tmp.genus.live, street.canopy$genus, nomatch=nogenus), "EqName"]
            ## add up total canopy area applying correct equation form to each
            tmp.can0 <- sum(((((eq.form=="quad")*(street.canopy[match(tmp.genus.live, street.canopy$genus, nomatch=nogenus), "a"]+
                                       (street.canopy[match(tmp.genus.live, street.canopy$genus, nomatch=nogenus), "b"]*tmp.dbh0.live)+
                                       (street.canopy[match(tmp.genus.live, street.canopy$genus, nomatch=nogenus), "c"]*tmp.dbh0.live^2)))/2)^2)*pi, na.rm=T)+
              sum(((((eq.form=="cub")*(street.canopy[match(tmp.genus.live, street.canopy$genus, nomatch=nogenus), "a"]+
                                         (street.canopy[match(tmp.genus.live, street.canopy$genus, nomatch=nogenus), "b"]*tmp.dbh0.live)+
                                         (street.canopy[match(tmp.genus.live, street.canopy$genus, nomatch=nogenus), "c"]*tmp.dbh0.live^2)+
                                        (street.canopy[match(tmp.genus.live, street.canopy$genus, nomatch=nogenus), "d"]*tmp.dbh0.live^3)))/2)^2)*pi, na.rm=T)+
              sum(((((eq.form=="lin")*(street.canopy[match(tmp.genus.live, street.canopy$genus, nomatch=nogenus), "a"]+
                                          (street.canopy[match(tmp.genus.live, street.canopy$genus, nomatch=nogenus), "b"]*tmp.dbh0.live)))/2)^2)*pi, na.rm=T)+
              sum((((eq.form=="loglogw1")*exp(street.canopy[match(tmp.genus.live, street.canopy$genus, nomatch=nogenus), "a"]+
                                         (street.canopy[match(tmp.genus.live, street.canopy$genus, nomatch=nogenus), "b"]*
                                            log(log(tmp.dbh0.live+1)+
                                                  street.canopy[match(tmp.genus.live, street.canopy$genus, nomatch=nogenus), "c"]/2)))/2)^2)*pi, na.rm=T)
  
            ### update track of productivity/biomass/canopy for year e
            npp.track <- c(npp.track, sum(tmp.biom1)-sum(tmp.biom0.live))
            biom.track <- c(biom.track, sum(tmp.biom0.live))
            can.track <- c(can.track, tmp.can0)
            
            ### now determine which of the previously dead are replanted this year
            if(sum(tree.samp==0 & delay==0)>0){ ## if there are dead trees here that qualify for replanting
              replanted <- rbinom(length(tree.samp), 1, replant.factor) ## figure which get the green light for replanting
              replanted[delay!=0] <- 0 ## only replant where delay clock is up
              replanted[tree.samp!=0] <- 0 ## only replant trees that are dead
              tmp.genus[replanted==1] <- as.character(sample(street[record.good==1, genus], sum(replanted==1))) ## swap in a new genus chosen at random
              tree.samp <- tree.samp+(replanted*5) ## replace selected dead trees that win the replacement lottery
#               print(paste("replanted", sum(replanted), "dead trees of genus", tmp.genus[replanted==1]))
              delay[replanted==1] <- sample(seq(delay.factor[1],delay.factor[2]),size = sum(replanted),replace=T) ## set a new delay clock for the replanted trees
            }
            
            ## update death tracking and wind down/reset the replant delay for the non-killed trees
            deaths <- deaths+sum(kill.list) ## count the dead
            num.track <- c(num.track, length(tree.samp[tree.samp>0])) ## count the living
            delay[kill.list==0 & delay==0] <- sample(seq(delay.factor[1],delay.factor[2]),size = sum(kill.list==0 & delay==0),replace=T) ## reset clocks that have wound down
            delay[kill.list==0 & delay>0] <- delay[kill.list==0 & delay>0]-1 ## count down the replant delay clock for things not killed this year
#             print(paste("there are", length(tree.samp[tree.samp>0]), "live trees by the end of year", e))
  
          } ## end of loop for successive years in the resim
            
            ### save out the results in lists
            # dbh.sav[[pix]][[a]] <- tree.samp ## updated tree sample after all 36 years of morts + growth
            npp.box[[pix]][[a]] <- npp.track
            biom.box[[pix]][[a]] <- biom.track
            expand.plant.num[[pix]] <- c(expand.plant.num[[pix]], newplants) ## just a vector of length yr.run how many newplantings you got in this resim
            deaths.sav[[pix]] <- c(deaths.sav[[pix]], deaths)
            num.box[[pix]][[a]] <- num.track
            can.box[[pix]][[a]] <- can.track
            
          } ## if enough pixel simulations exist
          else{  ## if too few successful simulations for this pixel
          npp.box[[pix]][[a]] <- NA
          biom.box[[pix]][[a]] <- NA
          expand.plant.num[[pix]] <- NA
          deaths.sav[[pix]] <- NA
          num.box[[pix]][[a]] <- NA
          can.box[[pix]][[a]] <- NA
#           print(paste("pix", index.track[procset[pix]], "failed, too few sims"))
            } ## end else statement
          } ## end pixel resim iteration loop
        
        if(pix%%500 == 0){print(paste("chunk", o, "resimmed pix", pix))} ## give status updates
        
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
  } ## if check on whether file has been processed
} ## scenario loop

## update the atwork file to release this job
check <- list.files("processed/boston/biom_street")
check <- check[grep(check, pattern=paste("index.track.scenario"))] ### version label here
check <- check[grep(check, pattern=paste0("V", resim.vers))]
library(stringr)
already <- str_match(check, paste0("BAU.(.*?).V", resim.vers))
already <- already[,2]
already <- already[!is.na(already)]
notyet <- obj.list[!(obj.list%in%already)]

## update atwork to remove any in-progress files that now have finished files on disk
atwork <- read.csv(paste("processed/boston/biom_street/resim.atwork", resim.vers, "csv", sep="."))
atwork <- atwork$at.work
atwork <- atwork[atwork%in%notyet] ## clear job records that have completed files on disk
l <- data.frame(at.work=atwork)
write.csv(l, file=paste("processed/boston/biom_street/resim.atwork", resim.vers, "csv", sep="."))
#####


###
### process resim results by scenario
#####
## get a general list of the resim results and the corresponding pixel index
sum.na <- function(x){sum(x, na.rm=T)}
scenario <- c("BAU", "oldies", "expand")
vers <- 6
resim.vers <- 8
realize <- 20
yr.run <- 34
preamb <- "processed/boston/biom_street/"
# preamb <- "/projectnb/buultra/atrlica/FragEVI/processed/boston/biom_street/"
library(abind) ## lets us build array boxes containing the maps/years/realizations organized well

### loop each scenario's results containers and process
for(s in 1:length(scenario)){
  print(paste("processing resim results for scenario", scenario[s]))
  obj.list <- list.files("processed/boston/biom_street/")
  obj.list <- obj.list[grep(obj.list, pattern=paste("scenario", scenario[s], sep="."))]
  obj.list <- obj.list[grep(obj.list, pattern=paste0("V", resim.vers))]
  index.list <- obj.list[grep(obj.list, pattern="index")]
  chunks <- strsplit(index.list, split = "[.]")
  chunks <- sapply(chunks, "[[" ,5)
  
  ### set up results containers
  biom.contain <- vector("list", realize)
  can.contain <- vector("list", realize)
  num.contain <- vector("list", realize)
  npp.contain <- vector("list", realize) ## each list element will be a matrix of pixIDx35 years, each top list element is a realization
  deaths.contain <- vector() 
  expand.contain <- vector() 
  pixID.record <- integer()
  npp.maps <- list()
  for(o in 1:length(chunks)){
    load(paste0(preamb, "npp.track.scenario.", scenario[s], ".", chunks[o], ".V", resim.vers, ".resim.sav", sep="")) ## npp.box
    load(paste0(preamb, "biom.track.scenario.", scenario[s], ".", chunks[o], ".V", resim.vers, ".resim.sav", sep="")) ## biom.box
    load(paste0(preamb, "num.track.scenario.", scenario[s], ".", chunks[o], ".V", resim.vers, ".resim.sav", sep="")) ## num.box
    load(paste0(preamb, "can.track.scenario.", scenario[s], ".", chunks[o], ".V", resim.vers, ".resim.sav", sep="")) ## can.box
    load(paste0(preamb, "death.track.scenario.", scenario[s], ".", chunks[o], ".V", resim.vers, ".resim.sav", sep="")) ## deaths.sav
    load(paste0(preamb, "expand.plant.num.scenario.", scenario[s], ".", chunks[o], ".V", resim.vers, ".resim.sav", sep="")) ## expand.plant.num
    load(paste0(preamb, "index.track.scenario.", scenario[s], ".", chunks[o], ".V", resim.vers, ".resim.sav", sep="")) ## pixID.track

#     dum <- sapply(npp.box, "[[", 1)
#     ded <- sapply(dum, is.null) ## use realization 1 to detect which resims didn't go through (appears stable across realizations)
    ## some of these bastards get an NA entry instead of a pinche NULL entry...
    for(r in 1:realize){ ## for every realization of the resimmed map
#       pix1r1 <- npp.box[[1]][[1]]
#       pix1r2 <- npp.box[[1]][[2]]
#       pix1r20 <- npp.box[[1]][[20]]
#       ss1 <- sapply(npp.box, "[[", 1) ## this is a list of realization 1 for all pix
#       ss2 <- sapply(npp.box, "[[", 2) ## this is a list of realization 2 for all pix
#       ss20 <- sapply(npp.box, "[[", 20) ## this is a list of realization 2 for all pix
#       plot(pix1r1, ss1[[1]])
#       plot(pix1r2, ss2[[1]])
#       plot(pix1r20, ss20[[1]]) ## these are all identical vectors, each is realization r of this chunk
      rel.npp <- sapply(npp.box, "[[", r) ## pick out realization r from the list of pixels*years
      rel.biom <- sapply(biom.box, "[[", r)
      rel.can <- sapply(can.box, "[[", r)
      rel.num <- sapply(num.box, "[[", r)
      if(r==1){ ## if this is the first realization processed in this chunk...
        ### a robust means of picking out and recording which pixel resimms worked (examine first realization as an archetype)
        ded1 <- sapply(rel.npp, is.null)
        ded1 <- which(ded1)
        ded2 <- sapply(rel.npp, "[[", 1)
        ded2 <- which(is.na(ded2))
        ded <- c(ded1, ded2) ## ID everything that is either NA or NULL
        pixID.record <- c(pixID.record, pixID.track[-ded]) ## update the vector of pixel IDs included in the resim record
      }
      ### filter bad pixels in this realization list
      ### unlist to single vector
      ### bind to whole-map vector for this realization as element in in ***.contain
      rel.npp <- rel.npp[-ded]
      npp.raw <- unlist(rel.npp) ## just pure vector of pixel resim histories, no entries for NULL
      npp.raw <- npp.raw[!is.na(npp.raw)] ## pick out any NAs that appear (do not have length 35)
      npp.contain[[r]] <- c(npp.contain[[r]], npp.raw)
      rel.biom <- rel.biom[-ded]
      biom.raw <- unlist(rel.biom)
      biom.raw <- biom.raw[!is.na(npp.raw)]
      biom.contain[[r]] <- c(biom.contain[[r]], biom.raw)
      rel.can <- rel.can[-ded]
      can.raw <- unlist(rel.can)
      can.raw <- can.raw[!is.na(npp.raw)]
      can.contain[[r]] <- c(can.contain[[r]], can.raw)
      rel.num <- rel.num[-ded]
      num.raw <- unlist(rel.num)
      num.raw <- num.raw[!is.na(npp.raw)]
      num.contain[[r]] <- c(num.contain[[r]], num.raw)  
    } ## loop for different realizations r of resimulator
  ## deaths/expand are sums per pixel across all resim years, recorded for each realization
  deaths.contain <- c(deaths.contain, deaths.sav[-ded])
  expand.contain <- c(expand.contain, expand.plant.num[-ded]) ## just get a raw vector and reconstitute to matrix later
  print(paste("finished parsing all resim realizations of chunk", chunks[o]))
  } ## loop for pixel chunks[o]

  ### sanity check: what do these "raw" list containers contain
  # length(npp.contain)
  # length(npp.contain[[1]])/35 ## this is all 35 year maps of a single realization  
  # unlist(lapply(npp.contain, length)) ## consistent lengths
  # unlist(lapply(biom.contain, length)) ## consistent but about 100 extra pix compared to npp ressims
  
  ### repackage raw vectors as pixelXyear matrices and write to disk
  mapmaker <- function(x){ ## take the raw vector and reconstitute into a map pixXyear matrix with pixID markers
    ana <- matrix(x, ncol=yr.run+1, byrow = T)
    ana <- cbind(pixID.record, ana)
    return(ana)
    } 
  npp.maps <- lapply(npp.contain, mapmaker) ## now we have a list of the map matricies with a list element for each realization
  biom.maps <- lapply(biom.contain, mapmaker)
  can.maps <- lapply(can.contain, mapmaker)
  num.maps <- lapply(num.contain, mapmaker)
  deaths.rel <- matrix(unlist(deaths.contain), nrow=length(pixID.record), byrow=T) ## rows=pixels, col=realization(sum all years)
  expand.rel <- matrix(unlist(expand.contain), nrow=length(pixID.record), byrow=T)
  
  print(paste("writing collated results of", scenario[s], "to disk"))
  save(npp.maps, file=paste0(preamb, "results/", scenario[s], ".V", resim.vers, ".npp.maps.sav"))
  save(biom.maps, file=paste0(preamb, "results/", scenario[s], ".V", resim.vers, ".biom.maps.sav"))
  save(can.maps, file=paste0(preamb, "results/", scenario[s], ".V", resim.vers, ".can.maps.sav"))
  save(num.maps, file=paste0(preamb, "results/", scenario[s], ".V", resim.vers, ".num.maps.sav"))
  write.csv(deaths.rel, file=paste0(preamb, "results/", scenario[s], ".V", resim.vers, ".deaths.maps.csv"))
  write.csv(expand.rel, file=paste0(preamb, "results/", scenario[s], ".V", resim.vers, ".expand.maps.csv"))
} ## loop for scenario[s]
#####


# ## Exploratory of resim results
# #####
# ## load up ancillary map data
# library(raster)
# library(rgdal)
# biom <- raster("processed/boston/bos.biom30m.tif")
# aoi <- raster("processed/boston/bos.aoi30m.tif")
# biom <- crop(biom, aoi)
# can <- raster("processed/boston/bos.can.redux30m.tif")
# isa <- raster("processed/boston/bos.isa30m.tif")
# biom.dat <- as.data.table(as.data.frame(biom))
# aoi.dat <- as.data.table(as.data.frame(aoi))
# can.dat <- as.data.table(as.data.frame(can))
# isa.dat <- as.data.table(as.data.frame(isa))
# lulc <- raster("processed/boston/bos.lulc30m.lumped.tif")
# lulc.dat <- as.data.table(as.data.frame(lulc))
# biom.dat <- cbind(biom.dat, aoi.dat, can.dat, isa.dat, lulc.dat)
# biom.dat[,pix.ID:=seq(1:dim(biom.dat)[1])]
# nonfor <- biom.dat[bos.aoi30m>800 & !(bos.lulc30m.lumped %in% c(1,4,5,6)) & bos.biom30m<20000, pix.ID] ### identify pixel ID's that are nonforest
# # hyb <- as.data.table(read.csv("processed/results/hybrid.results.V6.csv"))
# hyb.t <- as.data.table(as.data.frame(raster("processed/results/hybrid.V7.median.tif")))
# hyb.t[,pix.ID:=seq(1, dim(hyb.t)[1])]
# biom.dat <- merge(biom.dat, hyb.t, by="pix.ID")
# dim(biom.dat[bos.aoi30m>800,]) ## 136667 pix in the AOI
# 
# ## BAU scenario
### process list containers
npp.box <- abind(npp.maps, along=3) ## stack up each list element (matrix)
apply(npp.box, MARGIN = 3, sum)/2000/1000 ## total 35 years of productivity by realization
apply(npp.box, MARGIN=2, mean)/2000 ### mean per-pixel npp by year
apply(npp.box, MARGIN=1, mean) ## mean NPP across years, by pixel
dd <- apply(npp.box, MARGIN=c(2,3), sum) ## total full-map NPP in each yearXrealization
for(i in 1:20){
  plot(dd[,i]/2000/1000, main=paste("map NPP, realization", i)) ## why is realization 1 so weird?
}

# biom.dat[bos.aoi30m>800, sum(hybrid.V7.median, na.rm=T)/2000/1000] ### 9.5 ktC in ~2007
# bau.npp <- read.csv("processed/results/BAU.V7.npp.trendmap.csv")
# bau.can <- read.csv("processed/results/BAU.V7.can.trendmap.csv")
# bau.biom <- read.csv("processed/results/BAU.V7.biom.trendmap.csv")
# bau.num <- read.csv("processed/results/BAU.V7.num.trendmap.csv")
# 
# ## BAU NPP start/finish
# biom.dat <- merge(biom.dat, bau.npp[,c("pixID.track", "V2", "V36")], by.x="pix.ID", by.y="pixID.track", all.x=T)
# dim(biom.dat[bos.aoi30m>800,]) ## still 136667
# names(biom.dat)[8:9] <- c("BAU.start.npp", "BAU.finish.npp")
# dim(biom.dat[bos.aoi30m>800 & !(is.na(BAU.start.npp))]) ## 77479 (excluded all non dev/hdres/ldres from resim)
# ### does the predicted npp at the start of the resim resemble our estimated hybrid npp 
# plot(biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), hybrid.V7.median], 
#      biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), BAU.start.npp])
# abline(a=0, b=1, col="red") ## fairly close
# plot(biom.dat[bos.aoi30m>800 & !is.na(BAU.finish.npp), hybrid.V7.median], 
#      biom.dat[bos.aoi30m>800 & !is.na(BAU.finish.npp), BAU.finish.npp]) ## a possible slight decline?
# abline(a=0, b=1, col="blue")
# biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), sum(hybrid.V7.median, na.rm=T)/2000/1000] ## 6.0 ktC in hybrid npp IN THE RESIMMED PIX
# biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), sum(BAU.start.npp, na.rm=T)/2000/1000] ## 6.0 ktC in start of resim
# biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), sum(BAU.finish.npp, na.rm=T)/2000/1000] ## 5.9 ktC, slight decline over time
# 
# ## BAU canopy start/finish
# biom.dat <- merge(biom.dat, bau.can[,c("pixID.track", "V2", "V36")], by.x="pix.ID", by.y="pixID.track", all.x=T)
# names(biom.dat)[10:11] <- c("BAU.start.can", "BAU.finish.can")
# plot(biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), bos.can.redux30m*bos.aoi30m], 
#      biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), BAU.start.can])
# abline(a=0, b=1, col="red")
# ## not even close, way overpredicted by the model
# plot(biom.dat[bos.aoi30m>800 & !is.na(BAU.finish.npp), bos.can.redux30m*bos.aoi30m], 
#      biom.dat[bos.aoi30m>800 & !is.na(BAU.finish.npp), BAU.finish.can])
# abline(a=0, b=1, col="blue") ## the simulator distinctly overpredicts canopy cover
# biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), sum(bos.can.redux30m*bos.aoi30m, na.rm=T)/1E4] ## 1.7 kha of canopy in the resimmed pix
# biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), sum(bos.aoi30m, na.rm=T)/1E4] ## 7.0 kha total resimmed area = ~24.3% canopy cover in resimmed pix
# biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), sum(BAU.start.can, na.rm=T)/1E4] ## 3.7 kha canopy area start in the resimmed area
# biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), sum(BAU.finish.can, na.rm=T)/1E4] ## 3.7 kha canopy finish in the resimmed area
# 
# biom.dat[bos.aoi30m>800, sum(bos.can.redux30m*bos.aoi30m, na.rm=T)]/biom.dat[bos.aoi30m>800, sum(bos.aoi30m)] ## 25.6% canopy in the AOI
# biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), length(bos.aoi30m)] ## 77479 pix resimmed
# biom.dat[bos.aoi30m>800, length(bos.aoi30m)] ## 136667 pix with data
# 
# ## BAU biomass start/finish
# biom.dat <- merge(biom.dat, bau.biom[,c("pixID.track", "V2", "V36")], by.x="pix.ID", by.y="pixID.track", all.x=T)
# names(biom.dat)[12:13] <- c("BAU.start.biom", "BAU.finish.biom")
# plot(biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), bos.biom30m],
#      biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), BAU.start.biom])
# abline(a=0, b=1, col="red") ## LIKE A FUCKING GLOVE
# plot(biom.dat[bos.aoi30m>800 & !is.na(BAU.finish.npp), bos.biom30m],
#      biom.dat[bos.aoi30m>800 & !is.na(BAU.finish.npp), BAU.finish.biom])
# abline(a=0, b=1, col="blue")
# biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), sum(bos.biom30m, na.rm=T)/2000/1000] ## 162 ktC measured in the resimmed area, compare 357 ktC in AOI
# biom.dat[bos.aoi30m>800, sum(BAU.finish.biom, na.rm=T)/2000/1000] ## 174 ktC by the end of simulation, 7% increase
# 
# ## tree number start/finish
# biom.dat <- merge(biom.dat, bau.num[,c("pixID.track", "V2", "V36")], by.x="pix.ID", by.y="pixID.track", all.x=T)
# names(biom.dat)[14:15] <- c("BAU.start.num", "BAU.finish.num")
# # plot(biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), bos.biom30m],
# #      biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), BAU.start.biom])
# # abline(a=0, b=1, col="red") ## LIKE A FUCKING GLOVE
# ## can only compare start/finish, haven't exported median tree numbers for the map
# plot(biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), BAU.start.num],
#      biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), BAU.finish.num])
# abline(a=0, b=1, col="blue") ## a little marginal loss in there
# biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), sum(BAU.start.num, na.rm=T)] ## 581879 simmed trees to start
# biom.dat[bos.aoi30m>800, sum(BAU.finish.num, na.rm=T)] ## 552854 ktC by the end of simulation, ~5% decrease 
# 
# ## time series of BAU
# sum.na <- function(x){sum(x, na.rm=T)}
# bau.npp.t <- apply(bau.npp[,3:37], MARGIN = 2, FUN = sum.na)
# plot(1:35, bau.npp.t/2000)
# last(bau.npp.t)/2000; first(bau.npp.t)/2000
# 
# bau.can.t <- apply(bau.can[,3:37], MARGIN = 2, FUN = sum.na)
# plot(1:35, bau.can.t/1E4)
# last(bau.can.t)/1E4; first(bau.can.t)/1E4
# 
# bau.biom.t <- apply(bau.biom[,3:37], MARGIN = 2, FUN = sum.na)
# plot(1:35, bau.biom.t/2000/1000)
# last(bau.biom.t)/2000/1000; first(bau.biom.t)/2000/1000
# 
# bau.num.t <- apply(bau.num[,3:37], MARGIN = 2, FUN = sum.na)
# plot(1:35, bau.num.t) ## after stabilizing, settles into an equilibrium
# last(bau.num.t)/1000; first(bau.num.t)/1000
# ### interesting: There's about a 5% loss in he first 3 years as part of the "burn in" phase for equilibrating mortality and 
# ### delayed replanting: net loss at first as the mortalities kick in (about ~2%/yr) and then you have a large enough
# ### pool of dead-awaiting-replanting that the numbers then stable off as the new replants each year begin to replace the new deaths. 
# ### The model reads as an "excess" of mortalities in the first three years because it doesn't see any pool of dead-awaiting-replanting
# ### to draw from until it has recorded a large enough pool of mortalities for itself.
# 
# 
# ### Oldies scenario
# biom.dat[bos.aoi30m>800, sum(hybrid.V7.median, na.rm=T)/2000/1000] ### 10.8 ktC in ~2007
# old.npp <- read.csv("processed/results/oldies.V7.npp.trendmap.csv")
# old.can <- read.csv("processed/results/oldies.V7.can.trendmap.csv")
# old.biom <- read.csv("processed/results/oldies.V7.biom.trendmap.csv")
# old.num <- read.csv("processed/results/oldies.V7.num.trendmap.csv")
# 
# ## oldies NPP start/finish
# biom.dat <- merge(biom.dat, old.npp[,c("pixID.track", "V2", "V36")], by.x="pix.ID", by.y="pixID.track", all.x=T)
# dim(biom.dat[bos.aoi30m>800,]) ## still 136667
# names(biom.dat)[16:17] <- c("old.start.npp", "old.finish.npp")
# dim(biom.dat[bos.aoi30m>800 & !(is.na(old.start.npp))]) ## 77479
# plot(biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), hybrid.V7.median],
#      biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), old.start.npp])
# abline(a=0, b=1, col="red")
# plot(biom.dat[bos.aoi30m>800 & !is.na(old.finish.npp), hybrid.V7.median],
#      biom.dat[bos.aoi30m>800 & !is.na(old.finish.npp), old.finish.npp])
# abline(a=0, b=1, col="blue") ## maybe an increase
# biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), sum(hybrid.V7.median, na.rm=T)/2000/1000] ## 6.0 ktC IN THE RESIMMED AREA
# biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), sum(old.start.npp, na.rm=T)/2000/1000] ## 6.0 ktC
# biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), sum(old.finish.npp, na.rm=T)/2000/1000] ## 6.6 ktC  ### a ~10% increase
# 
# ## oldies canopy start/finish
# biom.dat <- merge(biom.dat, old.can[,c("pixID.track", "V2", "V36")], by.x="pix.ID", by.y="pixID.track", all.x=T)
# names(biom.dat)[18:19] <- c("old.start.can", "old.finish.can")
# plot(biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), bos.can.redux30m*bos.aoi30m],
#      biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), old.start.can])
# abline(a=0, b=1, col="red") ## overpredicted
# plot(biom.dat[bos.aoi30m>800 & !is.na(old.finish.npp), bos.can.redux30m*bos.aoi30m],
#      biom.dat[bos.aoi30m>800 & !is.na(old.finish.npp), old.finish.can])
# abline(a=0, b=1, col="blue")
# biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), sum(bos.can.redux30m*bos.aoi30m, na.rm=T)/1E4] ## 1.7 kha in the resimmed area
# biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), sum(old.start.can, na.rm=T)/1E4] ## 3.7 kha start in the resimmed area
# biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), sum(old.finish.can, na.rm=T)/1E4] ## 4.6 kha finish in the resimmed area, 25% increase
# 
# ## oldies biomass start/finish
# biom.dat <- merge(biom.dat, old.biom[,c("pixID.track", "V2", "V36")], by.x="pix.ID", by.y="pixID.track", all.x=T)
# names(biom.dat)[20:21] <- c("old.start.biom", "old.finish.biom")
# plot(biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), bos.biom30m],
#      biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), old.start.biom])
# abline(a=0, b=1, col="red") ## LIKE A FUCKING GLOVE
# plot(biom.dat[bos.aoi30m>800 & !is.na(old.finish.npp), bos.biom30m],
#      biom.dat[bos.aoi30m>800 & !is.na(old.finish.npp), old.finish.biom])
# abline(a=0, b=1, col="blue")
# biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), sum(bos.biom30m, na.rm=T)/2000/1000] ## 162 ktC in the resimmed area compare 357 ktC in AOI
# biom.dat[bos.aoi30m>800, sum(old.finish.biom, na.rm=T)/2000/1000] ## 235 ktC by the end of simulation, 45% increase!! (slow down rate you are sending trees to the chipper...)
# 
# ## tree number start/finish
# biom.dat <- merge(biom.dat, old.num[,c("pixID.track", "V2", "V36")], by.x="pix.ID", by.y="pixID.track", all.x=T)
# names(biom.dat)[22:23] <- c("old.start.num", "old.finish.num")
# # plot(biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), bos.biom30m],
# #      biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), BAU.start.biom])
# # abline(a=0, b=1, col="red") ## LIKE A FUCKING GLOVE
# ## can only compare start/finish, haven't exported median tree numbers for the map
# plot(biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), old.start.num],
#      biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), old.finish.num])
# abline(a=0, b=1, col="blue") ## a little marginal loss in there
# biom.dat[bos.aoi30m>800 & !is.na(old.start.npp), sum(old.start.num, na.rm=T)] ## 581879 simmed trees to start
# biom.dat[bos.aoi30m>800, sum(old.finish.num, na.rm=T)] ## 558854 ktC by the end of simulation, ~4% decrease, so about the same number, slightly more alive than BAU
# 
# ## time series of oldies
# old.npp.t <- apply(old.npp[,3:37], MARGIN = 2, FUN = sum.na)
# plot(1:35, old.npp.t/2000)
# last(old.npp.t)/2000; first(old.npp.t)/2000
# 
# old.can.t <- apply(old.can[,3:37], MARGIN = 2, FUN = sum.na)
# plot(1:35, old.can.t/1E4)
# last(old.can.t)/1E4; first(old.can.t)/1E4
# 
# old.biom.t <- apply(old.biom[,3:37], MARGIN = 2, FUN = sum.na)
# plot(1:35, old.biom.t/2000/1000)
# last(old.biom.t)/2000/1000; first(old.biom.t)/2000/1000
# 
# old.num.t <- apply(old.num[,3:37], MARGIN = 2, FUN = sum.na)
# plot(1:35, old.num.t) # after stabilizing, slow draw upwards as more trees age into protection zone
# last(old.num.t)/1000; first(old.num.t)/1000
# 
# 
# ### Expand scenario
# biom.dat[bos.aoi30m>800, sum(hybrid.V7.median, na.rm=T)/2000/1000] ### 9.5 ktC in ~2007
# exp.npp <- read.csv("processed/results/expand.V7.npp.trendmap.csv")
# exp.can <- read.csv("processed/results/expand.V7.can.trendmap.csv")
# exp.biom <- read.csv("processed/results/expand.V7.biom.trendmap.csv")
# exp.num <- read.csv("processed/results/expand.V7.num.trendmap.csv")
# 
# ## expand NPP start/finish
# biom.dat <- merge(biom.dat, exp.npp[,c("pixID.track", "V2", "V36")], by.x="pix.ID", by.y="pixID.track", all.x=T)
# dim(biom.dat[bos.aoi30m>800,]) ## still 136667
# names(biom.dat)[24:25] <- c("exp.start.npp", "exp.finish.npp")
# dim(biom.dat[bos.aoi30m>800 & !(is.na(exp.start.npp))]) ## 77479
# plot(biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), hybrid.V7.median],
#      biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), exp.start.npp])
# abline(a=0, b=1, col="red")
# plot(biom.dat[bos.aoi30m>800 & !is.na(exp.finish.npp), hybrid.V7.median],
#      biom.dat[bos.aoi30m>800 & !is.na(exp.finish.npp), exp.finish.npp])
# abline(a=0, b=1, col="blue")
# biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), sum(hybrid.V7.median, na.rm=T)/2000/1000] ## 6.0 ktC IN THE RESIMMED AREA
# biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), sum(exp.start.npp, na.rm=T)/2000/1000] ## 6.0 ktC
# biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), sum(exp.finish.npp, na.rm=T)/2000/1000] ## 6.9 ktC  ### 15% increase over time (more than large tree perservation)
# 
# ## expand canopy start/finish
# biom.dat <- merge(biom.dat, exp.can[,c("pixID.track", "V2", "V36")], by.x="pix.ID", by.y="pixID.track", all.x=T)
# names(biom.dat)[26:27] <- c("exp.start.can", "exp.finish.can")
# plot(biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), bos.can.redux30m*bos.aoi30m],
#      biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), exp.start.can])
# abline(a=0, b=1, col="red")
# ## not even close, way overpredicted by the model
# plot(biom.dat[bos.aoi30m>800 & !is.na(exp.finish.npp), bos.can.redux30m*bos.aoi30m],
#      biom.dat[bos.aoi30m>800 & !is.na(exp.finish.npp), exp.finish.can])
# abline(a=0, b=1, col="blue")
# biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), sum(bos.can.redux30m*bos.aoi30m, na.rm=T)/1E4] ## 1.7 kha in the resimmed area
# biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), sum(bos.aoi30m, na.rm=T)/1E4] ## 7.0 kha total resimmed area = 24.3% canopy
# biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), sum(exp.start.can, na.rm=T)/1E4] ## 3.7 kha start in the resimmed area
# biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), sum(exp.finish.can, na.rm=T)/1E4] ## 4.2 kha finish in the resimmed area, 13% increase, less than large tree preservation
# 
# ## Expand biomass start/finish
# biom.dat <- merge(biom.dat, exp.biom[,c("pixID.track", "V2", "V36")], by.x="pix.ID", by.y="pixID.track", all.x=T)
# names(biom.dat)[28:29] <- c("exp.start.biom", "exp.finish.biom")
# plot(biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), bos.biom30m],
#      biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), exp.start.biom])
# abline(a=0, b=1, col="red") ## LIKE A FUCKING GLOVE
# plot(biom.dat[bos.aoi30m>800 & !is.na(exp.finish.npp), bos.biom30m],
#      biom.dat[bos.aoi30m>800 & !is.na(exp.finish.npp), exp.finish.biom])
# abline(a=0, b=1, col="blue")
# biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), sum(bos.biom30m, na.rm=T)/2000/1000] ## 162 ktC in the resimmed area compare 357 ktC in AOI
# biom.dat[bos.aoi30m>800, sum(exp.finish.biom, na.rm=T)/2000/1000] ## 191 ktC by the end of simulation, 18% increase (less than large tree preservation)
# 
# ## tree number start/finish
# biom.dat <- merge(biom.dat, exp.num[,c("pixID.track", "V2", "V36")], by.x="pix.ID", by.y="pixID.track", all.x=T)
# names(biom.dat)[30:31] <- c("exp.start.num", "exp.finish.num")
# # plot(biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), bos.biom30m],
# #      biom.dat[bos.aoi30m>800 & !is.na(BAU.start.npp), BAU.start.biom])
# # abline(a=0, b=1, col="red") ## LIKE A FUCKING GLOVE
# ## can only compare start/finish, haven't exported median tree numbers for the map
# plot(biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), exp.start.num],
#      biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), exp.finish.num])
# abline(a=0, b=1, col="blue") ## marginal gain in the pixels, as expected
# biom.dat[bos.aoi30m>800 & !is.na(exp.start.npp), sum(exp.start.num, na.rm=T)] ## 580768 simmed trees to start
# biom.dat[bos.aoi30m>800, sum(exp.finish.num, na.rm=T)] ## 673189 ktC by the end of simulation, 16% gain -- alright this is the marginal ability of the city to add trees to its mix via street buffers
# 673189-580768 ## packed an additional 92k trees in, well short of goal
# 
# 
# ## time series of expand
# exp.npp.t <- apply(exp.npp[,3:37], MARGIN = 2, FUN = sum.na)
# plot(1:35, exp.npp.t/2000)
# last(exp.npp.t)/2000; first(exp.npp.t)/2000
# plot(1:35, exp.npp.t/2000, col="blue", ylim=c(5500, 7000))
# points(1:35, old.npp.t/2000, col="black")
# points(1:35, bau.npp.t/2000, col="red")
# 
# exp.can.t <- apply(exp.can[,3:37], MARGIN = 2, FUN = sum.na)
# plot(1:35, exp.can.t/2000)
# last(exp.can.t)/1E4; first(exp.can.t)/1E4
# plot(1:35, old.can.t/1E4)
# points(1:35, exp.can.t/1E4, col="blue")
# points(1:35, bau.can.t/1E4, col="red")
# 
# exp.biom.t <- apply(exp.biom[,3:37], MARGIN = 2, FUN = sum.na)
# plot(1:35, exp.biom.t/2000/1000)
# last(exp.biom.t)/2000/1000; first(exp.biom.t)/2000/1000
# plot(1:35, old.biom.t/2000/1000)
# points(1:35, exp.biom.t/2000/1000, col="blue")
# points(1:35, bau.biom.t/2000/1000, col="red")
# 
# exp.num.t <- apply(exp.num[,3:37], MARGIN = 2, FUN = sum.na)
# plot(1:35, exp.num.t) ## as expected, jumps in the begining up to a cruising altitude
# last(old.num.t)/1000; first(old.num.t)/1000
#####



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
#####
