################
#########
### final script for estimating tree distribution based on street tree records and running on the cluster
library(data.table)
library(raster)
library(rgdal)
library(rgeos)

setwd("/projectnb/buultra/atrlica/FragEVI/")

### biomass equation biom~dbh
b0 <- -2.48
b1 <- 2.4835 ## these are eastern hardwood defaults
biom.pred <- function(x){exp(b0+(b1*log(x)))}

## get the street tree record prepped
street <- read.csv("docs/ian/Boston_Street_Trees.csv") ## Street tree resurvey, biomass 2006/2014, species
street$ann.npp <- (street$biomass.2014-street$biomass.2006)/8
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
# a <- street[,length(ann.npp)/dim(street)[1], by=genus]
# a <- a[order(a$V1, decreasing = T),]
# # sum(a$V1[1:12]) ### top 12 is 94% of trees present
# focus <- a$genus[1:12]
street[, delta.dbh:=dbh.2014-dbh.2006]
street[, record.good:=0]
street[is.finite(dbh.2014) & is.finite(dbh.2006) & dead.by.2014==0 & delta.dbh>=0, record.good:=1] 
street[record.good==1, at.biom.2014:=biom.pred(dbh.2014)]
street[record.good==1, at.biom.2006:=biom.pred(dbh.2006)]
street[record.good==1, at.delta.biom:=at.biom.2014-at.biom.2006]
street[record.good==1, at.npp.ann:=at.delta.biom/8]
street[record.good==1, at.npp.ann.rel:=at.npp.ann/at.biom.2006]
street[record.good==1,range(dbh.2006, na.rm=T)]
street[dbh.2006<5, record.good:=0] ## filter for tiny trees, full record is now 2455

### biomass increment model based on dbh (exponential function)
mod.biom.rel <- summary(lm(log(at.npp.ann.rel)~log(dbh.2006), data=street[record.good==1 & ann.npp!=0])) #r2=0.54! so relative %change biomass per year is the most predictive so far
# plot(street[record.good==1, log(at.biom.2006)], street[record.good==1, log(at.npp.ann.rel)])

## read in biomass and canopy data
biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m kg-biomass to 30m pixel
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- crop(biom, aoi) ## biomass was slightly buffered, need to clip to match canopy fraction raster
biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[,aoi:=as.vector(getValues(aoi))]
biom.dat[aoi<800, bos.biom30m:=NA] ### cancel values in cells that are not at least ~90% complete coverage in AOI
can <- raster("processed/boston/bos.can30m.tif")
can.dat <- as.data.table(as.data.frame(can))
biom.dat <- cbind(biom.dat, can.dat)
biom.dat[,index:=1:dim(biom.dat)[1]]

## prep street tree data and biomass data for processing (small biomass first)
# clean <- street[record.good==1 & dbh.2006>=5,] # get a good street tree set ready
clean <- street[record.good==1,] # get a good street tree set ready
setkey(clean, at.biom.2006)

### first part: process the low biomass with distributions drawn from the unmodified street set
runme <- biom.dat[!is.na(bos.biom30m) & !is.na(bos.can30m) & bos.biom30m<20000 & bos.biom30m>5 & bos.can30m>0,] #99k, filter for very small biomass and 0 canopy

## set up containers
cage.num.trees <- list()
cage.ann.npp <- list()
cage.dbh <- list()
index.track <- integer()
biom.track <- numeric()

## loop each row
for(t in 1:dim(runme)[1]){
  ann.npp <- numeric()
  num.trees <- numeric()
  x <- 0
  q <- 0
  while(x<100 & q<2000){ ## select 100 workable samples, or quit after 4000 attempts
    ## grab a random number of randomly selected trees
    grasp <- clean[at.biom.2006<runme[t, bos.biom30m], .(dbh.2006, at.biom.2006)]
    n <- min(c(80, dim(grasp)[1])) ## if you get a really tiny biomass it might try to sample too many rows
    grasp <- grasp[sample(dim(grasp)[1], size=n),]
    w=grasp[1, at.biom.2006] ## keep cummulative tally of biomass
    d=1 # keep track of the number of trees
    while(w<(0.9*runme[t, bos.biom30m])){ ## keep adding trees until you just get over the target biomass
      w=w+grasp[d+1, at.biom.2006]
      d=d+1
    }
    ### if you've gone too far or if you've packed them in too tight, ditch this sample
    if(grasp[1:d, sum((((dbh.2006/2)^2)*pi)/1E4)]<(runme[t, bos.can30m]*3.6) & w<(1.10*runme[t, bos.biom30m])){ ## if the BA density is low enough & didn't overshoot biomass too much
      ann.npp <- c(ann.npp, sum(grasp[1:d, at.biom.2006]*exp((mod.biom.rel$coefficients[2]*log(grasp[1:d, dbh.2006]))+mod.biom.rel$coefficients[1])))
      num.trees <- c(num.trees, d)
      x <- x+1
#       print(paste("recorded", x))
    }
    q=q+1 ## if you can't find a combination that works, try q times to get develop a sample and if you can't fuck this pixel
#     print(paste("attempted", q))
  }
  if(q<2000){
    cage.ann.npp[[t]] <- ann.npp
    cage.num.trees[[t]] <- num.trees
    cage.dbh[[t]] <- grasp[1:d, dbh.2006]
    biom.track <- c(biom.track, runme[t, bos.biom30m])
    index.track <- c(index.track, runme[t, index])
    print(paste("finished pixel", t))
  } else{
    cage.ann.npp[[t]] <- 9999 ## could not find a solution
    cage.num.trees[[t]] <- 9999
    cage.dbh[[t]] <- 9999
    biom.track <- c(biom.track, runme[t, bos.biom30m])
    index.track <- c(index.track, runme[t, index])
    print(paste("pixel", t, "error"))
  }
}

## appears working
save(cage.ann.npp, file=paste("processed/boston/ann.npp.street.small"))
save(cage.num.trees, file=paste("processed/boston/num.trees.street.small"))
save(cage.dbh, file=paste("processed/boston/dbh.street.small"))
save(biom.track, file=paste("processed/boston/biom.track.street.small"))
save(cage.ann.npp, file=paste("processed/boston/index.track.street.small"))





####
### now specificially tackle the high-biomass areas
### modify the street tree record to amplify high-dbh trees
clean <- street[record.good==1 & dbh.2006>=5,]
bar <- clean[,quantile(dbh.2006, probs= 0.80)] ## what part of the high end to amplify
bar2 <- clean[,quantile(dbh.2006, probs= 0.90)]
clean <- rbind(clean, clean[dbh.2006>bar,], clean[dbh.2006>bar,], clean[dbh.2006>bar,], clean[dbh.2006>bar,], clean[dbh.2006>bar,]) ## add the top 20% in another 3 times
setkey(clean, at.biom.2006)

### first part: process the low biomass with distributions drawn from the unmodified street set
runme <- biom.dat[!is.na(bos.biom30m) & !is.na(bos.can30m) & bos.biom30m>=20000 & bos.biom30m<30000,] #6k

## set up containers
cage.num.trees <- list()
cage.ann.npp <- list()
cage.dbh <- list()
index.track <- integer()
biom.track <- numeric()

## loop each row
for(t in 1:dim(runme)[1]){
  ann.npp <- numeric()
  num.trees <- numeric()
  x <- 0
  q <- 0
  while(x<100 & q<1000){ ## select 100 workable samples, or quit after 4000 attempts
    ## grab a random number of randomly selected trees
    grasp <- clean[, .(dbh.2006, at.biom.2006)]
    n <- min(c(80, dim(grasp)[1])) ## if you get a really tiny biomass it might try to sample too many rows
    grasp <- grasp[sample(dim(grasp)[1], size=n),]
    w=grasp[1, at.biom.2006] ## keep cummulative tally of biomass
    d=1 # keep track of the number of trees
    while(w<(0.9*runme[t, bos.biom30m])){ ## keep adding trees until you just get over the target biomass
      w=w+grasp[d+1, at.biom.2006]
      d=d+1
    }
    ### if you've gone too far or if you've packed them in too tight, ditch this sample
    if(grasp[1:d, sum((((dbh.2006/2)^2)*pi)/1E4)]<(runme[t, bos.can30m]*3.6) & w<(1.10*runme[t, bos.biom30m])){ ## if the BA density is low enough & didn't overshoot biomass too much
      ann.npp <- c(ann.npp, sum(grasp[1:d, at.biom.2006]*exp((mod.biom.rel$coefficients[2]*log(grasp[1:d, dbh.2006]))+mod.biom.rel$coefficients[1])))
      num.trees <- c(num.trees, d)
      x <- x+1
    }
    q=q+1 ## if you can't find a combination that works, try q times to get develop a sample and if you can't fuck this pixel
  }
  if(q<1000){
    cage.ann.npp[[t]] <- ann.npp
    cage.num.trees[[t]] <- num.trees
    cage.dbh[[t]] <- grasp[1:d, dbh.2006]
    biom.track <- c(biom.track, runme[t, bos.biom30m])
    index.track <- c(index.track, runme[t, index])
    print(paste("finished pixel", t))
  } else{
    cage.ann.npp[[t]] <- 9999 ## could not find a solution
    cage.num.trees[[t]] <- 9999
    cage.dbh[[t]] <- 9999
    biom.track <- c(biom.track, runme[t, bos.biom30m])
    index.track <- c(index.track, runme[t, index])
    print(paste("pixel", t, "error"))
  }
}

# hist(unlist(cage.ann.npp[[3]]))
# hist(unlist(cage.dbh[[3]]))
# hist(unlist(cage.num.trees[[3]]))
# biom.track
## appears working
save(cage.ann.npp, file=paste("processed/boston/ann.npp.street.big"))
save(cage.num.trees, file=paste("processed/boston/num.trees.street.big"))
save(cage.dbh, file=paste("processed/boston/dbh.street.big"))
save(biom.track, file=paste("processed/boston/biom.track.street.big"))
save(cage.ann.npp, file=paste("processed/boston/index.track.street.big"))