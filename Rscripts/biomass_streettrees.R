################
#########
### final script for estimating tree distribution based on street tree records and running on the cluster
library(data.table)
library(raster)
library(rgdal)
library(rgeos)

setwd("/projectnb/buultra/atrlica/FragEVI/")

######
###### Approach 5: Reconfigure street tree analysis for more effective search
## Look at cell biomass sample smartly from dbh data, bounded to stop excessive density or BA to reach the biomass total
## street tree simulator
################
#########
### biomass equation biom~dbh
b0 <- -2.48
b1 <- 2.4835 ## these are eastern hardwood defaults
biom.pred <- function(x){exp(b0+(b1*log(x)))}

## get the street tree record prepped
street <- read.csv("docs/ian/Boston_Street_Trees.csv") ## Street tree resurvey, biomass 2006/2014, species
street <- as.data.table(street)
street[is.na(Species), Species:="Unknown unknown"]
street[Species=="Unknown", Species:="Unknown unknown"]
street[Species=="unknown", Species:="Unknown unknown"]
street[Species=="Purpleleaf plum", Species:="Prunus cerasifera"]
street[Species=="Bradford pear", Species:="Pyrus calleryana"]
street[Species=="Liriondendron tulipifera", Species:="Liriodendron tulipifera"]
street[Species=="Thornless hawthorn", Species:="Crataegus crus-galli"]
genspec <- strsplit(as.character(street$Species), " ")
gen <- unlist(lapply(genspec, "[[", 1))
street[,genus:=gen] ## 40 genera
street[, record.good:=0] ## ID records with good data quality
street[is.finite(dbh.2014) & is.finite(dbh.2006) & dead.by.2014==0 & Health!="Poor", record.good:=1] #2603 records good
street[record.good==1, biom.2014:=biom.pred(dbh.2014)]
street[record.good==1, biom.2006:=biom.pred(dbh.2006)]
street[record.good==1, npp.ann:=(biom.2014-biom.2006)/8]
street[record.good==1, npp.ann.rel:=npp.ann/biom.2006]
street[record.good==1,range(dbh.2006, na.rm=T)] 
street[record.good==1, range(npp.ann, na.rm=T)] ## 202 records show negative growth -- questionable 2006 dbh record, cull
street[record.good==1 & npp.ann<0, record.good:=0] ## 2401

### biomass increment model based on dbh (exponential function)
# mod.biom.rel <- summary(lm(log(npp.ann.rel)~log(dbh.2006), data=street[record.good==1 & npp.ann!=0])) ## exclude 4 no-growth records
# plot(street[record.good==1, log(biom.2006)], street[record.good==1, log(npp.ann.rel)]) # r2 = 0.54, slope significant
## there's <10 records below 5cm with a lot of leverage, cut them
mod.biom.rel <- summary(lm(log(npp.ann.rel)~log(dbh.2006), data=street[record.good==1 & npp.ann!=0 & dbh.2006>=5,])) ## exclude 4 no-growth records
# plot(street[record.good==1 & dbh.2006>=5, log(biom.2006)], street[record.good==1 & dbh.2006>=5, log(npp.ann.rel)]) # r2 = 0.55, looks nicer
street[dbh.2006<5, record.good:=0] ## filter out the handfull of truly tiny trees


## read in biomass and canopy data
biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m kg-biomass to 30m pixel
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- crop(biom, aoi) ## biomass was slightly buffered, need to clip to match canopy fraction raster
biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[,aoi:=as.vector(getValues(aoi))]
can <- raster("processed/boston/bos.can30m.tif")
can.dat <- as.data.table(as.data.frame(can))
biom.dat <- cbind(biom.dat, can.dat)
biom.dat[,index:=1:dim(biom.dat)[1]] ## master pixel index for stack of biom/aoi/can, 354068 pix
# p <- ecdf(biom.dat[!is.na(bos.biom30m) & bos.biom30m>10,bos.biom30m]) ## cdf of cell biomass

## prep street tree data and biomass data for processing (small biomass first)
ba.pred <- function(x){(x/2)^2*pi*0.0001} ## get BA per stump from dbh
clean <- street[record.good==1,] # get a good street tree set ready
clean[, ba:=ba.pred(dbh.2006)]
setkey(clean, biom.2006) # 2390 records in final selection, dbh range 5-112 cm
clean <- clean[order(clean$dbh.2006),]
clean[,rank:=seq(from=1, to=dim(clean)[1])] ## all the trees have a fixed size rank now (used in adjusting sampling weights)

### set up data
runme <- biom.dat[!is.na(bos.biom30m) & bos.biom30m>10 & !is.na(aoi) & !is.na(bos.can30m) & bos.can30m>0,] ## 108k

## test chunks
# runme <- runme[bos.biom30m>30000,]
# runme <- runme[sample(1:dim(runme)[1], size=1000),]

### NOTE To parallelize this process (for the lower biomass pixels), the script checks for already-written chunks of results, and then tries to produce the next chunk
### calling the script multiple times will result in multiple successive chunks of pixels being run at the same time

runme.x <- runme[1:10000,] ## initialize first chunk

## parallel process: check to see if any containers have been written to disk, if not queue up the next chunk
check <- list.files("processed/boston/biom_street")
check <- check[grep(check, pattern="ann.npp.street.v3")]
already <- sub(".*weighted\\.", "", check)
if(length(check)!=0){ ## ie if you detect that results have already been written to disk, take the next chunk
  y <- as.numeric(max(already)) ## figure out what's been written to disk already
  runme.x <- runme[(y+1):(y+10000),] ## reset the next chunk
  if((y+10000)>(dim(runme)[1])){ ## if you're at the end of the file, only grab up to the last row
    runme.x <- runme[(y+1):dim(runme)[1],]
  }
}

## set up containers
cage.num.trees <- list()
cage.ann.npp <- list()
cage.dbh <- list()
cage.genus <- list()
index.track <- rep(9999, dim(runme.x)[1])
biom.track <- rep(9999, dim(runme.x)[1])
cage.wts <- list()
cage.biom.sim <- list()
attempts.track <- rep(9999, dim(runme.x)[1])
proc.track <- rep(9999, dim(runme.x)[1])

## create an empty save file to warn the script next time that this chunk is being worked on
if(length(check)!=0){
  stor <- (y+10000) ## will name the end files up to its max, but it might not be that long
  save(cage.ann.npp, file=paste("processed/boston/biom_street/ann.npp.street.v3.weighted", stor, sep=".")) 
} else{
  stor <- 10000
  save(cage.ann.npp, file=paste("processed/boston/biom_street/ann.npp.street.v3.weighted", stor, sep="."))
}

## for dynamic weighting
incr=300 ## incrememnt to change weighting
k=200 ## constant, arbitrary top of prob weights, higher-->steeper change in sampling weight
D=k/incr ## drop increment, sensitive to number of tries to make while adjusting the sampling weights

## loop each row
for(t in 1:dim(runme.x)[1]){
  ann.npp <- numeric()
  num.trees <- numeric()
  biom.sim.track <- numeric()
  wts.track <- integer()
  cage.genus[[t]] <- list()
  cage.dbh[[t]] <- list()
  
  x <- 0 ## count successful simulations
  q <- 0 ## count attempts since last successful sample
  g <- 0 ## weight adjustment start
  z <- 0 ## total number of sample iterations
  jeez <- 0 ## if it's grinding through a slow one
  
  ### this approach uses a moving/shrinking window on the street tree records as biomass starts to get bigger
  ## has the disadvantage of supressing the number of trees at higher biomass (start to oversample very large trees), would need some tuning to get the right sampling window
  # ## 0 figure where and how much to sample
  # cell.q <- p(runme[t, bos.biom30m]) ## this is the percentile of the cell biomass
  # spread <- min((1-cell.q), 0.5) ## interval to sample in dbh records, cap it off if the spread gets bigger than 50%
  # window.l <- max(0, (cell.q-spread))
  # window.h <- cell.q+spread
  # if((spread)==0.5){window.h <- 1} ## basic scheme is that for biomass>median, restrict to (rarer) higher window, but if biom<median, sample whole window
  # dbh.lims <- quantile(clean$dbh.2006, probs=c(window.l, window.h)) ## window of dbh to select based on cell biomass
  # grasp <- clean[biom.2006<runme[t, bos.biom30m] &dbh.2006<=dbh.lims[2] & dbh.2006>=dbh.lims[1], .(dbh.2006, biom.2006, ba)] ## this is a street tree record restricted based on cell biomass
  
  dens.lim <- ceiling(0.106*runme.x[t,aoi*bos.can30m]) ## maximum number to select, probably only restrictive in low-canopy cells
  
  ### this approach takes the full street trees record and samples it with differing weights if the algorithm fails
  grasp <- clean[biom.2006<(1.1*runme.x[t, bos.biom30m]), .(dbh.2006, biom.2006, ba, rank, genus)] ## excluding single trees too big to fit cell biomass
  
  ### attempt to sample street trees to fill the cells
  while(x<100 & q<1000){ ## select 100 workable samples, or quit after q attempts with no success
    # ## grab a sample of grasp (no weighting)
    # samp <- grasp[sample(dim(grasp)[1], size=dens.lim, replace=T),] ## sample grasp with replacement only up to density limit
    
    ## OR: figure out the dynamic weighting to use based on previous failures
    wts=(k-(g*D))+(((g*D)/(dim(grasp)[1]-1))*((grasp$rank)-1))
    # sample grasp with replacement only up to density limit, with specified weights based on iterative reweighting
    samp <- grasp[sample(dim(grasp)[1], size=dens.lim, replace=F, prob = wts),] 
    w=samp[1, biom.2006] ## keep cummulative tally of biomass
    d=1 # keep track of the number of trees
    while(w<(0.9*runme.x[t, bos.biom30m]) & d<dens.lim){ ## keep adding trees until you just get over the target biomass or run out of records
      w=w+samp[d+1, biom.2006]
      d=d+1
    }
    ### if this is too many trees too tightly packed (over 40m2/ha BA) or not enough biomass in the sample, readjust weights
    if(((samp[1:d, sum(ba)]/runme.x[t, aoi*bos.can30m])*1E4)>40 | w<(0.90*runme.x[t, bos.biom30m])){
      if(g<incr){
        g=g+1 ## readjust the sample weights
        if(g==300){print(paste("pixel", runme[t, index], "weights at maximum"))}
        q=0 ## reset attempt timeout clock
      }
    }
    if(((samp[1:d, sum(ba)]/runme.x[t, aoi*bos.can30m])*1E4)<40 & w<(1.10*runme.x[t, bos.biom30m])){ ## if the BA density is low enough & didn't overshoot biomass too much
      x <- x+1 ## record successful sample
      ann.npp <- c(ann.npp, sum(samp[1:d, biom.2006]*exp((mod.biom.rel$coefficients[2]*log(samp[1:d, dbh.2006]))+mod.biom.rel$coefficients[1])))
      num.trees <- c(num.trees, d)
      biom.sim.track <- c(biom.sim.track, w) ## simulated biomass in this sample
      cage.dbh[[t]][[x]] <- samp[1:d, dbh.2006] ## which dbhs did you select
      cage.genus[[t]][[x]] <- samp[1:d, genus] # running list of which genera you select
      wts.track <- c(wts.track, g) ## weights used to get this sample
      if(q>300){print(paste("finally got sample", x, "after", q, "tries")); jeez <- 1} ## to give status reports for difficult/long-process samples
      q=0 ## reset counter if you've found some solution
      if(jeez==1 & x%%10==0){print(paste("patience, just got sample", x))}
      if(x>30 & g>250){q <- (-500)} ## give more room if this is a hard one and you've managed by chance to get a pretty good sample going
    }
    q=q+1 ## record number of attempts since last succesful simulation
    if(q>600 & q%%100==0){print(paste("attempt clock out in", (1000-q)/100))}
    z=z+1 ## record total number of sample iterations
  }
  cage.ann.npp[[t]] <- ann.npp
  cage.num.trees[[t]] <- num.trees
  cage.wts[[t]] <- wts.track
  biom.track <- c(biom.track, runme[t, bos.biom30m])
  index.track <- c(index.track, runme[t, index])
  cage.biom.sim[[t]] <- biom.sim.track
  attempts.track <- c(attempts.track, z)
  if(q<1000){
    proc.track <- c(proc.track, 1)
    print(paste("finished pixel", runme[t, index], "=", t))
  } else{
    proc.track <- c(proc.track, 0) ## record as incomplete processing
    print(paste("pixel", runme[t, index], "failed"))
  }
}

## when complete dump everything back into the save file
save(cage.ann.npp, file=paste("processed/boston/biom_street/ann.npp.street.v3.weighted", stor, sep="."))
save(cage.num.trees, file=paste("processed/boston/biom_street/num.trees.street.v3.weighted", stor, sep="."))
save(cage.dbh, file=paste("processed/boston/biom_street/dbh.street.v3.weighted", stor, sep="."))
save(cage.genus, file=paste("processed/boston/biom_street/genus.street.v3.weighted", stor, sep="."))
save(biom.track, file=paste("processed/boston/biom_street/biom.track.street.v3.weighted", stor, sep="."))
save(index.track, file=paste("processed/boston/biom_street/index.track.street.v3.weighted", stor, sep="."))
save(proc.track, file=paste("processed/boston/biom_street/index.track.street.v3.weighted", stor, sep="."))
save(cage.biom.sim, file=paste("processed/boston/biom_street/cage.biom.sim.street.v3.weighted", stor, sep="."))
save(attempts.track, file=paste("processed/boston/biom_street/attempts.track.street.v3.weighted", stor, sep="."))





#### OLD SCRIPT: Treats <20k and >20k pixels with different sampling distributions
# ### biomass equation biom~dbh
# b0 <- -2.48
# b1 <- 2.4835 ## these are eastern hardwood defaults
# biom.pred <- function(x){exp(b0+(b1*log(x)))}
# 
# ## get the street tree record prepped
# street <- read.csv("docs/ian/Boston_Street_Trees.csv") ## Street tree resurvey, biomass 2006/2014, species
# street$ann.npp <- (street$biomass.2014-street$biomass.2006)/8
# street <- as.data.table(street)
# street[is.na(Species), Species:="Unknown unknown"]
# street[Species=="Unknown", Species:="Unknown unknown"]
# street[Species=="unknown", Species:="Unknown unknown"]
# street[Species=="Purpleleaf plum", Species:="Prunus cerasifera"]
# street[Species=="Bradford pear", Species:="Pyrus calleryana"]
# street[Species=="Liriondendron tulipifera", Species:="Liriodendron tulipifera"]
# street[Species=="Thornless hawthorn", Species:="Crataegus crus-galli"]
# genspec <- strsplit(as.character(street$Species), " ")
# gen <- unlist(lapply(genspec, "[[", 1))
# street[,genus:=gen] ## 40 genera
# # a <- street[,length(ann.npp)/dim(street)[1], by=genus]
# # a <- a[order(a$V1, decreasing = T),]
# # # sum(a$V1[1:12]) ### top 12 is 94% of trees present
# # focus <- a$genus[1:12]
# street[, delta.dbh:=dbh.2014-dbh.2006]
# street[, record.good:=0]
# street[is.finite(dbh.2014) & is.finite(dbh.2006) & dead.by.2014==0 & delta.dbh>=0, record.good:=1] # exclude fishy looking records
# street[record.good==1, at.biom.2014:=biom.pred(dbh.2014)]
# street[record.good==1, at.biom.2006:=biom.pred(dbh.2006)]
# street[record.good==1, at.delta.biom:=at.biom.2014-at.biom.2006]
# street[record.good==1, at.npp.ann:=at.delta.biom/8]
# street[record.good==1, at.npp.ann.rel:=at.npp.ann/at.biom.2006]
# street[record.good==1,range(dbh.2006, na.rm=T)]
# street[dbh.2006<5, record.good:=0] ## filter for tiny trees, full record is now 2455
# 
# ### biomass increment model based on dbh (exponential function)
# mod.biom.rel <- summary(lm(log(at.npp.ann.rel)~log(dbh.2006), data=street[record.good==1 & ann.npp!=0])) #r2=0.54! so relative %change biomass per year is the most predictive so far
# # plot(street[record.good==1, log(at.biom.2006)], street[record.good==1, log(at.npp.ann.rel)])
# 
# ## read in biomass and canopy data
# biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m kg-biomass to 30m pixel
# aoi <- raster("processed/boston/bos.aoi30m.tif")
# biom <- crop(biom, aoi) ## biomass was slightly buffered, need to clip to match canopy fraction raster
# biom.dat <- as.data.table(as.data.frame(biom))
# biom.dat[,aoi:=as.vector(getValues(aoi))]
# biom.dat[aoi<800, bos.biom30m:=NA] ### cancel values in cells that are not at least ~90% complete coverage in AOI
# can <- raster("processed/boston/bos.can30m.tif")
# can.dat <- as.data.table(as.data.frame(can))
# biom.dat <- cbind(biom.dat, can.dat)
# biom.dat[,index:=1:dim(biom.dat)[1]]
# 
# ## prep street tree data and biomass data for processing (small biomass first)
# # clean <- street[record.good==1 & dbh.2006>=5,] # get a good street tree set ready
# clean <- street[record.good==1,] # get a good street tree set ready
# setkey(clean, at.biom.2006)
# 
# ### first part: process the low biomass with distributions drawn from the unmodified street set
# ### NOTE To parallelize this process (for the lower biomass pixels), the script checks for already-written chunks of results, and then tries to produce the next chunk
# ### calling the script multiple times will result in multiple successive chunks of pixels being run at the same time
# 
# ### code commented to exclusively run the "big tree" code modification below
# runme <- biom.dat[!is.na(bos.biom30m) & !is.na(bos.can30m) & bos.biom30m<20000 & bos.biom30m>5 & bos.can30m>0,] #99k, filter for very small biomass and 0 canopy
# runme.x <- runme[1:10000,] ## initialize first chunk
# 
# ## parallel process: check to see if any containers have been written to disk, if not queue up the next chunk
# check <- list.files("processed/boston/biom_street")
# check <- check[grep(check, pattern="ann.npp.street.small")]
# already <- substr(check, start = 22, stop=26)
# if(length(check)!=0){ ## ie if you detect that results have already been written to disk, take the next chunk
#   y <- as.numeric(max(already)) ## figure out what's been written to disk already
#   runme.x <- runme[(y+1):(y+10000),] ## reset the next chunk
#   if((y+10000)>(dim(runme)[1])){ ## if you're at the end of the file, only grab up to the last row
#     runme.x <- runme[(y+1):dim(runme)[1],]
#   }
# }
# 
# ## set up containers
# cage.num.trees <- list()
# cage.ann.npp <- list()
# cage.dbh <- list()
# cage.genus <- list()
# index.track <- rep(9999, dim(runme.x)[1])
# biom.track <- rep(9999, dim(runme.x)[1])
# 
# ## create an empty save file to warn the script next time that this chunk is being worked on
# if(length(check)!=0){
#   stor <- (y+10000) ## will name the end files up to its max, but it might not be that long
#   save(cage.ann.npp, file=paste("processed/boston/biom_street/ann.npp.street.small", stor, sep=".")) 
# } else{
#   stor <- 10000
#   save(cage.ann.npp, file=paste("processed/boston/biom_street/ann.npp.street.small", stor, sep="."))
# }
# 
# ## loop each row (pixel) of the chunk
# for(t in 1:dim(runme.x)[1]){
#   ann.npp <- numeric()
#   num.trees <- numeric()
#   cage.genus[[t]] <- list()
#   cage.dbh[[t]] <- list()
#   x <- 0
#   q <- 0
#   while(x<100 & q<2000){ ## select x workable samples, or quit after q attempts
#     ## grab a random number of randomly selected trees out of the street tree database
#     grasp <- clean[at.biom.2006<(1.10*runme.x[t, bos.biom30m]), .(dbh.2006, at.biom.2006, genus)] # don't select any trees that are bigger biomass than the pixel total
#     n <- min(c(80, dim(grasp)[1])) ## if pixel biomass is tiny, don't try to sample too many rows
#     grasp <- grasp[sample(dim(grasp)[1], size=n),]
#     w=grasp[1, at.biom.2006] ## cummulative tally of biomass
#     d=1 # track of the number of trees
#     while(w<(0.9*runme.x[t, bos.biom30m])){ ## keep adding trees until you just get over the target biomass
#       w=w+grasp[d+1, at.biom.2006]
#       d=d+1
#     }
#     ### check if you've got too much biomass or if you've packed them in too tight
#     if(grasp[1:d, sum((((dbh.2006/2)^2)*pi)/1E4)]<(runme.x[t, bos.can30m]*3.6) & w<(1.10*runme.x[t, bos.biom30m])){ ## if the BA density is low enough & didn't overshoot biomass too much
#       ann.npp <- c(ann.npp, sum(grasp[1:d, at.biom.2006]*exp((mod.biom.rel$coefficients[2]*log(grasp[1:d, dbh.2006]))+mod.biom.rel$coefficients[1])))
#       num.trees <- c(num.trees, d)
#       x <- x+1
#       cage.dbh[[t]][[x]] <- grasp[1:d, dbh.2006] ## which dbhs did you select
#       cage.genus[[t]][[x]] <- grasp[1:d, genus] # running list of which genera you select
#       #       print(paste("recorded", x))
#     }
#     q=q+1 ## record this as an attempt
# #     print(paste("attempted", q))
#   }
#   if(q<2000){
#     cage.ann.npp[[t]] <- ann.npp
#     cage.num.trees[[t]] <- num.trees
#     biom.track[t] <- runme.x[t, bos.biom30m]
#     index.track[t] <- runme.x[t, index]
#     print(paste("finished pixel", t))
#   } else{
#     cage.ann.npp[[t]] <- 9999 ## could not find a solution
#     cage.num.trees[[t]] <- 9999
#     biom.track[t] <- runme.x[t, bos.biom30m]
#     index.track[t] <- runme.x[t, index]
#     print(paste("pixel", t, "error"))
#   }
# }
# 
# ## when complete dump everything back into the save file
# save(cage.ann.npp, file=paste("processed/boston/biom_street/ann.npp.street.small", stor, sep="."))
# save(cage.num.trees, file=paste("processed/boston/biom_street/num.trees.street.small", stor, sep="."))
# save(cage.dbh, file=paste("processed/boston/biom_street/dbh.street.small", stor, sep="."))
# save(cage.genus, file=paste("processed/boston/biom_street/genus.street.small", stor, sep="."))
# save(biom.track, file=paste("processed/boston/biom_street/biom.track.street.small", stor, sep="."))
# save(index.track, file=paste("processed/boston/biom_street/index.track.street.small", stor, sep="."))
# 



# ####
# ### now specificially tackle the high-biomass areas
# ### modify the street tree record to amplify high-dbh trees
# clean <- street[record.good==1 & dbh.2006>=5,]
# bar <- clean[,quantile(dbh.2006, probs= 0.80)] ## what part of the high end to amplify
# bar2 <- clean[,quantile(dbh.2006, probs= 0.90)]
# clean <- rbind(clean, clean[dbh.2006>bar,], clean[dbh.2006>bar,], clean[dbh.2006>bar,], clean[dbh.2006>bar2,], clean[dbh.2006>bar2,]) ## add the top 20% in another 3 times, the top 10% x2
# setkey(clean, at.biom.2006)
# 
# ### first part: process the low biomass with distributions drawn from the unmodified street set
# runme <- biom.dat[!is.na(bos.biom30m) & !is.na(bos.can30m) & bos.biom30m>=20000 & bos.biom30m<30000,] #6k, ditch the small fraction that are too large to ever solve
# runme.x <- runme
# 
# ## set up containers
# cage.num.trees <- list()
# cage.ann.npp <- list()
# cage.dbh <- list()
# cage.genus <- list()
# index.track <- rep(9999, dim(runme.x)[1])
# biom.track <- rep(9999, dim(runme.x)[1])
# 
# ## loop each row
# for(t in 1:dim(runme.x)[1]){
#   ann.npp <- numeric()
#   num.trees <- numeric()
#   cage.genus[[t]] <- list()
#   cage.dbh[[t]] <- list()
#   x <- 0
#   q <- 0
#   while(x<100 & q<3000){ ## select 100 workable samples, or quit after 4000 attempts
#     ## grab a random number of randomly selected trees
#     grasp <- clean[, .(dbh.2006, at.biom.2006, genus)]
#     n <- min(c(80, dim(grasp)[1])) ## if you get a really tiny biomass it might try to sample too many rows
#     grasp <- grasp[sample(dim(grasp)[1], size=n),]
#     w=grasp[1, at.biom.2006] ## keep cummulative tally of biomass
#     d=1 # keep track of the number of trees
#     while(w<(0.9*runme[t, bos.biom30m])){ ## keep adding trees until you just get over the target biomass
#       w=w+grasp[d+1, at.biom.2006]
#       d=d+1
#     }
#     ### if you've gone too far or if you've packed them in too tight, ditch this sample
#     if(grasp[1:d, sum((((dbh.2006/2)^2)*pi)/1E4)]<(runme.x[t, bos.can30m]*3.6) & w<(1.10*runme.x[t, bos.biom30m])){ ## if the BA density is low enough & didn't overshoot biomass too much
#       ann.npp <- c(ann.npp, sum(grasp[1:d, at.biom.2006]*exp((mod.biom.rel$coefficients[2]*log(grasp[1:d, dbh.2006]))+mod.biom.rel$coefficients[1])))
#       num.trees <- c(num.trees, d)
#       x <- x+1
#       cage.dbh[[t]][[x]] <- grasp[1:d, dbh.2006] ## which dbhs did you select
#       cage.genus[[t]][[x]] <- grasp[1:d, genus] # running list of which genera you select
#     }
#     q=q+1 ## if you can't find a combination that works, try q times to get develop a sample and if you can't fuck this pixel
#   }
#   if(q<3000){
#     cage.ann.npp[[t]] <- ann.npp
#     cage.num.trees[[t]] <- num.trees
#     biom.track[t] <- runme.x[t, bos.biom30m]
#     index.track[t] <- runme.x[t, index]
#     print(paste("finished pixel", t))
#   } else{
#     cage.ann.npp[[t]] <- 9999 ## could not find a solution
#     cage.num.trees[[t]] <- 9999
#     biom.track[t] <- runme.x[t, bos.biom30m]
#     index.track[t] <- runme.x[t, index]
#     print(paste("pixel", t, "error"))
#   }
# }
# 
# save(cage.ann.npp, file=paste("processed/boston/biom_street/ann.npp.street.big"))
# save(cage.num.trees, file=paste("processed/boston/biom_street/num.trees.street.big"))
# save(cage.dbh, file=paste("processed/boston/biom_street/dbh.street.big"))
# save(biom.track, file=paste("processed/boston/biom_street/biom.track.street.big"))
# save(index.track, file=paste("processed/boston/biom_street/index.track.street.big"))
# save(cage.genus, file=paste("processed/boston/biom_street/genus.street.big"))



