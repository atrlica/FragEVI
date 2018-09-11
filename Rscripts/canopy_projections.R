library(raster)
library(data.table)



vers <- 4 ## which model run are we picking at
# 
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


### street tree reconstruction of results
# #### import results objects pulled from parallel processing on the cluster (chunks of 10k pixels)
# obj.dump <- list.files("processed/boston/biom_street/")
# npp.dump <- obj.dump[grep(obj.dump, pattern = paste("ann.npp.street.v", vers, ".weighted*", sep=""))]
# npp.dump.chunks <- sub('.*weighted.', '', npp.dump)
# npp.dump.chunks <- sub("\\.sav.*", "", npp.dump.chunks) ## later version get named .sav

load("processed/boston/biom_street/dbh.street.v4.weighted.10000.sav")
load("processed/boston/biom_street/ann.npp.street.v4.weighted.10000.sav")
median(cage.ann.npp[[1]])
quantile(cage.ann.npp[[1]], probs=c(0.45, 0.55))

#### a question arises: Is the population of trees to select the one that gets closest to the pixel biomass, or the one nearest to the median NPP estimate?
load("processed/boston/biom_street/biom.sim.street.v4.weighted.10000.sav")
quantile(cage.biom.sim[[1]], probs=c(0.45, 0.55))
j <- which(cage.biom.sim[[1]]>2638 & cage.biom.sim[[1]]<2650)
load("processed/boston/biom_street/biom.track.street.v4.weighted.10000.sav")
biom.track[1]

### mortality rate equation
mort <- function(x){(0.0008133*(x^2))-(0.0642407*x)+4.0614503}
x <- seq(1,100)
plot(x, mort(x))
# deathwatch <- mort(unlist(cage.dbh[[1]][j[1]]))
# rbinom(10, 100, 0.1)
# reap.me <- 0
# g <- 0
# while(g==0){
#   g <- rbinom(1, 1, deathwatch[1]/100)
#   reap.me=reap.me+1
# }
# reap.me
# 

## growth rate equation
b0 <- -2.48
b1 <- 2.4835 ## these are eastern hardwood defaults
biom.pred <- function(x){exp(b0+(b1*log(x)))}
street <- read.csv("docs/ian/Boston_Street_Trees.csv") ## Street tree resurvey, biomass 2006/2014, species
street <- as.data.table(street)
street[, record.good:=0] ## ID records with good data quality
street[is.finite(dbh.2014) & is.finite(dbh.2006) & dead.by.2014==0 & Health!="Poor", record.good:=1] #2603 records good
street[dbh.2006<5, record.good:=0] ## filter out the handfull of truly tiny trees
street[record.good==1, biom.2014:=biom.pred(dbh.2014)]
street[record.good==1, biom.2006:=biom.pred(dbh.2006)]
street[record.good==1, growth.ann:=(biom.2014-biom.2006)/8]
street[record.good==1, growth.ann.rel:=growth.ann/biom.2006]
street[record.good==1, dbh.ann:=(dbh.2014-dbh.2006)/8]
street[record.good==1, dbh.ann.rel:=dbh.ann/dbh.2006]
plot(street[record.good==1, dbh.2006], street[record.good==1, dbh.ann.rel])
mod.street.dbh.nls <- nls(dbh.ann.rel ~ exp(a + b * log(dbh.2006)), data=street[record.good==1,], start=list(a=0, b=0)) 
mm <- summary(mod.street.dbh.nls) ## RSE is pretty high 0.33
dbhg.pred <- function(x){exp(mm$coefficients[1]+(mm$coefficients[2]*log(x)))} ## function for predicting next year's relative dbh growth from this year's dbh



## loops for simulating growth+mortality in pixel dbh samples
## pull in all the simulator samples
vers <- 4 ## what version are we dealing with?
obj.list <- list.files("processed/boston/biom_street/")
obj.list <- obj.list[grep(obj.list, pattern=paste("dbh.street.v", vers, sep=""))]
obj.list <- sub('.*weighted\\.', '', obj.list)
obj.list <- as.integer(sub('\\..*', '', obj.list))

for(o in obj.list){
  load(paste("processed/boston/biom_street/dbh.street.v", vers, ".weighted.", o, ".sav", sep=""))
  dbh.sav <- cage.dbh
  for(pix in 1:length(cage.dbh)){
    if(!is.na(cage.biom.sim[[pix]])){ ## if a successful sample was simulated for this pixel
      for(a in 1:100){   ### resimulate dbh history in every pixel 100 times
        # biom.lims <- quantile(cage.biom.sim[[pix]], probs=c(0.40, 0.6)) ## figure out which of the simulations to draw and modify
        # j <- sample(which(cage.biom.sim[[pix]]>=biom.lims[1] & cage.biom.sim[[pix]]<=biom.lims[2]),1) ## get a random dbh sample from the the middle 10% of simulator results close to the target biomass
        npp.lims <- quantile(cage.ann.npp[[pix]], probs=c(0.40, 0.6)) ## figure out which of the simulations to draw and modify
        j <- sample(which(cage.ann.npp[[pix]]>=npp.lims[1] & cage.ann.npp[[pix]]<=npp.lims[2]),1) ## get a random dbh sample from the the middle 10% of simulator results close to the target biomass
        tree.samp <- unlist(cage.dbh[[pix]][j]) ## what initial trees are present in the pixel simulator result
        dbh.sav[[pix]][[a]] <- numeric() ## clear contents to store the new simulated results
        dbh.track <- numeric()
        for(y in 1:length(tree.samp)){ ## simulate each tree history separately
          d <- tree.samp[y]
          # track <- d[y]
          for(e in 1:36){ ## test each tree for 36 years (2006-2040)
            deathwatch <- mort(d) ## recalculate death probability
            if(rbinom(1, 1, deathwatch/100)){ ## if tree doesn't make it
              # print("he's dead Jim")
              d <- 5
            } else{
              # print("It's alive!!")
              d <- d*(1+dbhg.pred(d))
            }
            # track <- c(track, d)
            # ptrack <- c(ptrack, deathwatch)
          }
          dbh.track <- c(dbh.track, d)
          # plot(track)
          # track <- numeric()
          # track <- numeric()
          # plot(ptrack)
          # ptrack <- numeric()
        }
        dbh.sav[[pix]][[a]] <- dbh.track
      }
      print(paste("chunk", o, "resimmed pix", pix))
    } else{dbh.sav[[pix]][[a]] <- NA}
  }
 save(dbh.sav, paste("processed/boston/biom_street/dbh.resim.scenario-0.BAU.", o, ".sav", sep="")) 
}


## how many trees do you get in each re-sim  
num.resim.trees <- unlist(lapply(dbh.sav[[pix]], length))
num.sim.trees <- unlist(lapply(cage.dbh[[pix]], length))
hist(num.sim.trees)
hist(num.resim.trees)

plot(unlist(lapply(cage.dbh[[400]], median)), unlist(lapply(dbh.sav[[400]], median)))
abline(a=0, b=1)
