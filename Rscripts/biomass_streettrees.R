#########
### final script for estimating tree distribution based on street tree records and running on the cluster
#### NOTE for running on cluster, be sure to comment out Part 1 (analysis) and Part 3 (reconstruction) and let just the simulator run.

library(data.table)
library(raster)
library(rgdal)
library(rgeos)
library(lme4)

# setwd("/projectnb/buultra/atrlica/FragEVI/")


#### PART 0: TESTING THE SENSE OF THINGS
#####
# ### Some ground truthing using the street trees at corner of Granby and Bay State
# ## dbh measurements made 9/12/2018 (contrast biomass based on 2006/2007 data)
# biom1m <- raster("data/dataverse_files/bostonbiomass_1m.tif")
# bs1 <- readOGR("processed/boston/BayState1.shp")
# bs1.biom <- sum(unlist(extract(biom1m, bs1)))
# bs1.biom # 162 kg
# biom.inv(162) # 21cm dbh predicted, actual = 36
# 
# bs2 <- readOGR("processed/boston/BayState2.shp")
# bs2.biom <- sum(unlist(extract(biom1m, bs2)))
# biom.inv(bs2.biom) ##32 cm predicted, actual 38.4
# 
# bs3 <- readOGR("processed/boston/BayState3.shp")
# bs3.biom <- sum(unlist(extract(biom1m, bs3)))
# biom.inv(bs3.biom) #48 cm predicted, actual 46.5
# 
# bs4 <- readOGR("processed/boston/BayState4.shp")
# bs4.biom <- sum(unlist(extract(biom1m, bs4)))
# biom.inv(bs4.biom) #46 cm predited, actual 41.9
# 
# bs5 <- readOGR("processed/boston/BayState5.shp")
# bs5.biom <- sum(unlist(extract(biom1m, bs5)))
# biom.inv(bs5.biom) ## 21 cm predicted, actual 30
# 
# bs6 <- readOGR("processed/boston/BayState6.shp")
# bs6.biom <- sum(unlist(extract(biom1m, bs6)))
# biom.inv(bs6.biom) ## 20cm predicted, actual 33.8
# 
# bs7 <- readOGR("processed/boston/BayState7.shp")
# bs7.biom <- sum(unlist(extract(biom1m, bs7)))
# biom.inv(bs7.biom) ##22 cm predicted, actual 29.5
# 
# bs8 <- readOGR("processed/boston/BayState8.shp")
# bs8.biom <- sum(unlist(extract(biom1m, bs8)))
# biom.inv(bs8.biom) ## 14 cm predicted, 20.4cm actual
# 
# park1 <- readOGR("processed/boston/Park1.shp")
# park1.biom <- sum(unlist(extract(biom1m,park1)))
# biom.inv(park1.biom) ## 40 cm predcited, 50.3 cm actual
# 
# g1 <- readOGR("processed/boston/Granby1.shp")
# g1.biom <- sum(unlist(extract(biom1m, g1)))
# biom.inv(g1.biom) ## 40cm predicted, 52.5 cm actual

# ### get tree survey coordinates and translate to UTM
# pts <- readOGR("processed/boston/street_trees_06_coords.shp")
# pts <- pts[pts@data$V1!=0,] ## remove bad coords
# pts <- spTransform(pts, crs(biom1m))
# # plot(biom1m)
# # plot(pts, add=T) ### tight tight tight they overlay
# writeOGR(pts, "processed/boston/street_trees_06_coords_UTM.shp", 
#          layer="pts",
#          driver="ESRI Shapefile")
### visually, I'm not sure how knowing the street trees in 30m pixel will get us to a numerical check. 
## Not all trees in any pixel are street trees, not clear that survey is exhaustive in any pixel
#####


######### Part 1: Develop prediction for growth~dbh in street trees
## read data and clean up
#####
street <- read.csv("docs/ian/Boston_Street_Trees.csv") ## Street tree resurvey, biomass 2006/2014, species
street <- as.data.table(street)
street[, record.good:=0] ## ID records with good data quality
street[is.finite(dbh.2014) & is.finite(dbh.2006) & dead.by.2014==0 & Health!="Poor", record.good:=1] #2603 records good
street[record.good==1, biom.2014:=biom.pred(dbh.2014)]
street[record.good==1, biom.2006:=biom.pred(dbh.2006)]
street[record.good==1, npp.ann:=(biom.2014-biom.2006)/8]
street[record.good==1, npp.ann.rel:=npp.ann/biom.2006]
street[dbh.2006<5, record.good:=0] ## filter out the handfull of truly tiny trees
street[,delta.diam:=dbh.2014-dbh.2006]
street[,diam.rate:=delta.diam/8] ## annualized

### fix taxa and assign to categories of use
street[is.na(Species), Species:="Unknown unknown"]
street[Species=="Unknown", Species:="Unknown unknown"]
street[Species=="unknown", Species:="Unknown unknown"]
street[Species=="Purpleleaf plum", Species:="Prunus cerasifera"]
street[Species=="Bradford pear", Species:="Pyrus calleryana"]
street[Species=="Liriondendron tulipifera" | Species=="Liriodendron tulipifera ", Species:="Liriodendron tulipifera"]
street[Species=="Thornless hawthorn", Species:="Crataegus crus-galli"]
street[Species=="Ginko biloba", Species:="Ginkgo biloba"]
street[Species=="Cerceidiphyllum japonicum", Species:="Cercidiphyllum japonicum"]
street[Species=="Katsura", Species :="Cercidiphyllum spp."]
street[Species=="Pyrus cerasifera", Species :="Prunus cerasifera"]
street[Species=="Ulmus american", Species :="Ulmus americana"]
street[Species=="Acer palmatum atropurpureum", Species :="Acer palmatum"]
street[Species=="Crataegus crusgalli var. inermis", Species :="Crataegus crus-galli"]
street[Species=="Sophora japonica 'Regent'", Species :="Sophora japonica"]
street[Species=="Gleditsia triacanthos var. inerm", Species :="Gleditsia triacanthos"]
street[,Species:=(as.character(Species))]
table(street$Species)

genspec <- strsplit(as.character(street$Species), " ")
gen <- unlist(lapply(genspec, "[[", 1))
street[,genus:=gen] ## 40 genera

### make a simplified generic grouping for stats (excel sheet for explanation: genus.simp.definitions.xlsx)
street[,genus.simp:=genus]
street[genus%in%c("Amelanchier", "Crataegus"), genus.simp:="Malus"] ## not calling it "Malinae" because we are not including Pyrus
street[genus%in%c("Quercus", "Betula", "Carpinus", "Carya", "Corylus", "Fagus", "Ostrya"), genus.simp:="Fagales"]
street[genus%in%c("Aesculus", "Koelreuteria", "Acer"), genus.simp:="Sapindaceae"]
street[genus%in%c("Catalpa", "Celtis", "Cercidiphyllum", "Cornus", "Liquidambar", "Syringa",
                  "Halesia", "Liriodendron", "Magnolia", "Morus", "Pinus", "Populus", "Unknown"), genus.simp:="Other"] ## 129 of these, 84 w/o valid genus
street[genus%in%c("Cercis", "Maackia", "Robinia", "Sophora", "Gleditsia"), genus.simp:="Fabaceae"] ## 21 of these
table(street[record.good==1,genus.simp]) ## 15 taxonomic groups
# table(street$genus) ## 39 unique genera
write.csv(street, "processed/boston/street.trees.dbh.csv")

# boxplot(diam.rate~genus, data=street[record.good==1,])
# boxplot(diam.rate~genus.simp, data=street[record.good==1,])
# g <- street[record.good==1, median(diam.rate, na.rm=T), by=genus.simp]
# g[order(genus.simp),]
# street[record.good==1, range(dbh.2006, na.rm=T)] 
# street[record.good==1, range(npp.ann, na.rm=T)] ## 202 records show negative growth -- questionable 2006 dbh record?

# ### figure out a weighted average wood density for the unnassigned taxa given their prevalence in the data set
# ((street[Species=="Quercus rubra", length(Species)]*560)+
#     (street[Species=="Quercus alba", length(Species)]*600)+
#     (street[Species=="Quercus palustris", length(Species)]*580)+
#     (street[Species=="Quercus macrocarpa", length(Species)]*580)+
#     (street[Species=="Pyrus calleryana" | Species=="Pyrus", length(Species)]*603)+
#     (street[Species=="Ginkgo biloba", length(Species)]*455)+
#     (street[Species=="Malus", length(Species)]*610)+
#     (street[Species=="Ulmus americana" | Species=="Ulmus", length(Species)]*460))/
#   (street[Species%in%c("Quercus rubra", "Quercus alba","Quercus macrocarpa", 
#                        "Quercus palustris","Pyrus calleryana", "Pyrus", "Ginkgo biloba",
#                        "Malus", "Ulmus americana", "Ulmus"), length(Species)])
# ### weighted avg of 549 kg/m3
# 
# street[Species%in%c("Quercus rubra", "Quercus alba","Quercus macrocarpa", 
#                     "Quercus palustris","Pyrus calleryana", "Pyrus", "Ginkgo biloba",
#                     "Malus", "Ulmus americana", "Ulmus") & record.good==1, length(Species)] ## 329 good records that are not also represented by specific volume equations
# street[!(genus%in%c("Acer", "Zelkova", "Fraxinus", "Pinus", "Platanus", "Gleditsia", "Magnolia", "Liquidambar", "Celtis", "Tilia")) & 
#          record.good==1, length(genus)] ## 425 good records outside of the main genera, for which we have an average density representating 77% of them
# street[genus%in%c("Acer", "Zelkova", "Fraxinus", "Pinus", "Platanus", "Gleditsia", "Magnolia", "Liquidambar", "Celtis", "Tilia") & 
#          record.good==1, length(genus)] ## 2167 records that we can map to an urban-specific volume allometry, 83% of good stems
# street[record.good==1, length(genus)] ## 2592 records total

# hist(street[record.good==1, log(dbh.2006)])
# qqnorm(street[record.good==1, log(dbh.2006)])
# qqline(street[record.good==1, log(dbh.2006)])
# qqnorm(street[record.good==1, dbh.2006])
# qqline(street[record.good==1, (dbh.2006)])
# muf <- ecdf(street[record.good==1, dbh.2006])
# muf(40)
# muf(50)
# muf(60)

# ### a question: why is there so much periodicity in the dbh.2006 measures? Is this a function of diameter measured to nearest 1 inch?
# a <- street$dbh.2006
# a <- a[order(a)]
# a <- a[!is.na(a)]
# a <- a/2.54 ## convert to inches
# length(a) #3492
# length(a[a%%1==0]) ##3027
# 3027/3492 ## 87% are whole inches
# unique(a[a%%1==0]) ## basically every inch reading up to 60"
#####


### stem growth~dbh
#####
library(data.table)
street <- as.data.table(read.csv("processed/boston/street.trees.dbh.csv"))
par(mfrow=c(1,1), mar=c(4,4,3,1))
street[record.good==1,] ## 2592 records total
# dip <- as.data.frame(table(street[record.good==1, Species]))
# write.csv(dip, "docs/street.species.csv")

### basic summaries
street[record.good==1, quantile(dbh.2006, probs=c(0.05, 0.5, 0.95), na.rm=T)]
street[record.good==1, quantile(diam.rate, probs=c(0.05, 0.5, 0.95), na.rm=T)]

### now the growth modeling
### basics: how does delta diameter vary?
plot(street[record.good==1, dbh.2006], street[record.good==1, delta.diam])
summary(lm(delta.diam~dbh.2006, data=street[record.good==1,])) ## so about a 0.1 cm decline in delta per cm of diameter
plot(street[record.good==1, dbh.2006], street[record.good==1, diam.rate])
summary(lm(diam.rate~dbh.2006, data=street[record.good==1,])) ## 1cm/yr, subtract 0.1 cm/yr for each 10cm increase in size
hm <- (lm(diam.rate~poly(dbh.2006, degree=4), data=street[record.good==1,]))
points(street[record.good==1, dbh.2006], predict(hm), col="red", pch=16, cex=0.2)
summary(hm) ## low predictive, RSE 0.63, R2 0.11
street[,diam.rate.rel:=diam.rate/dbh.2006]
# plot(street[record.good==1, dbh.2006], street[record.good==1, diam.rate.rel]) ## same exponential-looking curve
### hyperbolic but maybe just because we are taking something nearly constant (annual dbh increment) and dividing by x

hm <- (lm(diam.rate~poly(dbh.2006, degree=3), data=street[record.good==1,]))
points(street[record.good==1, dbh.2006], predict(hm), col="purple", pch=16, cex=0.2)
summary(hm) ## low predictive, RSE 0.63, R2 0.11
hm <- (lm(diam.rate~poly(dbh.2006, degree=2), data=street[record.good==1,]))
points(street[record.good==1, dbh.2006], predict(hm), col="blue", pch=16, cex=0.2)
s.hm <- summary(hm) ## low predictive, RSE 0.63, R2 0.11, this is as good as we are likely to get
s.hm

## estimate proper raw coefficients
plot(street[record.good==1, dbh.2006], street[record.good==1, diam.rate])
hm <- (lm(diam.rate~poly(dbh.2006, degree=2, raw=T), data=street[record.good==1,]))
points(street[record.good==1, dbh.2006], predict(hm), col="blue", pch=16, cex=0.2)
s.hm <- summary(hm) ## low predictive, RSE 0.63, R2 0.11
points(street[record.good==1, dbh.2006], 
       (hm$coefficients[1]-(1.96*s.hm$coefficients[1,2]))+((hm$coefficients[2]-(1.96*s.hm$coefficients[2,2]))*street[record.good==1, dbh.2006]^1)+((hm$coefficients[3]-(1.96*s.hm$coefficients[3,2]))*street[record.good==1, dbh.2006]^2),
       col="green")
points(street[record.good==1, dbh.2006], 
       (hm$coefficients[1]+(1.96*s.hm$coefficients[1,2]))+((hm$coefficients[2]+(1.96*s.hm$coefficients[2,2]))*street[record.good==1, dbh.2006]^1)+((hm$coefficients[3]+(1.96*s.hm$coefficients[3,2]))*street[record.good==1, dbh.2006]^2),
       col="red") ## reasonable-ish, let's use this one
# plot(hm)
# table(street[record.good==1, genus])
mod.street.dbhdelta <- hm
save(mod.street.dbhdelta, file="processed/mod.street.dbhdelta.sav")

# ### contrast: FIA data
# live <- as.data.table(read.csv("processed/fia.live.stem.dbh.growth.csv"))
# live[,delta.diam:=DIAM_T1-DIAM_T0]
# plot(live$DIAM_T0, live$delta.diam) ## generally lower diameter gain
# summary(lm(delta.diam~DIAM_T0, data=live)) ## about a 0.03 cm decline in delta per cm of dbh (i.e. less sesitive to size overall)
# live[,diam.rate:=delta.diam/4.8] ## annualized
# plot(live[DIAM_T0>0, DIAM_T0], live[DIAM_T0>0, diam.rate], pch=15, cex=0.3, ylim=c(-1, 2))
# summary(lm(diam.rate~DIAM_T0, data=live[DIAM_T0>0 & STATUS==1,])) ## dumb model (data not clean) but about 0.08cm/yr, add 0.01 cm/yr for each 10cm increment
# live[,diam.rate.rel:=diam.rate/DIAM_T0]
# plot(live$DIAM_T0, live$diam.rate.rel)


## compare street and FIA scatters and lms
par(mfrow=c(1,2))
hm1 <- (lm(diam.rate~poly(dbh.2006, degree=1, raw=T), data=street[record.good==1,]))
hm2 <- (lm(diam.rate~poly(dbh.2006, degree=2, raw=T), data=street[record.good==1,])) ## I think this is the guy
hm3 <- (lm(diam.rate~poly(dbh.2006, degree=3, raw=T), data=street[record.good==1,]))
summary(hm1); summary(hm2); summary(hm3); hist(street[record.good==1, diam.rate])
AIC(hm1); AIC(hm2); AIC(hm3)
plot(street[record.good==1, dbh.2006], hm$residuals); hist(hm$residuals) ## perfect. wow.
plot(street[record.good==1, dbh.2006], street[record.good==1, diam.rate], 
     col="salmon", pch=15, cex=0.8, ylim=c(-1, 3), main="Street trees", xlab="DBH (cm)", ylab="DBH change (cm/yr)")
points(street[record.good==1, dbh.2006], predict(hm1), col="red", cex=0.2)
points(street[record.good==1, dbh.2006], predict(hm2), col="blue", cex=0.2)

### perhaps we can account for species differences?
boxplot(diam.rate~genus, data=street[record.good==1,])
boxplot(dbh.2006~genus, data=street[record.good==1,])
boxplot(diam.rate~genus.simp, data=street[record.good==1,])
boxplot(dbh.2006~genus.simp, data=street[record.good==1,])

#####
#### mixed effects accounting for taxa differences
library(lme4)
hm2.me <- lmer(diam.rate~poly(dbh.2006, degree=2, raw=T)+
                 (1|genus.simp), REML=F, data=street[record.good==1,])
summary(hm2.me)

hm2.me2 <- lmer(diam.rate~poly(dbh.2006, degree=2, raw=T)+
                  (1+dbh.2006|genus.simp), REML=F, data=street[record.good==1,])
summary(hm2.me2) ## coefficients a lot like the previous fixed model and random intercepts model
plot(street[record.good==1, dbh.2006], street[record.good==1, diam.rate], col="black", cex=0.5)
points(street[record.good==1 & genus=="Zelkova", dbh.2006], street[record.good==1 & genus=="Zelkova", diam.rate], col="red", cex=0.9, pch=16)
points(street[record.good==1, dbh.2006], predict(hm2.me2), col=street[record.good==1, genus.simp], cex=0.5, pch=16)
lines(1:120, coef(summary(hm2.me2))[1]+(coef(summary(hm2.me2))[2]*(1:120))+(coef(summary(hm2.me2))[3]*(1:120)^2), col="blue")

hm1.me2 <-lmer(diam.rate~poly(dbh.2006, degree=1, raw=T)+
                 (1+dbh.2006|genus.simp), REML=F, data=street[record.good==1,])
anova(hm1.me2, hm2.me2) ## helps to have the polynomial term

hm0.me2 <- lmer(diam.rate~1+
                  (1+dbh.2006|genus.simp), REML=F, data=street[record.good==1,])
anova(hm0.me2, hm2.me2) ## better to have a dbh term than not

## different RE formulation
hm2.me3 <-  lmer(diam.rate~poly(dbh.2006, degree=2, raw=T)+
                             (dbh.2006|genus.simp), REML=F, 
                 data=street[record.good==1,]) ## same results
hm1.me3 <-  lmer(diam.rate~poly(dbh.2006, degree=1, raw=T)+
                   (dbh.2006|genus.simp), REML=F, 
                 data=street[record.good==1,])
hm3.me3 <- lmer(diam.rate~poly(dbh.2006, degree=3, raw=T)+
                      (dbh.2006|genus.simp), REML=F, 
                data=street[record.good==1,]) ## kills sig of poly terms
hm0.me3 <- lmer(diam.rate~1+
                  (dbh.2006|genus.simp), REML=F, 
                data=street[record.good==1,])
anova(hm2.me3, hm3.me3) ## NS 3rd order
anova(hm1.me3, hm2.me3) ## 2nd order significant ***
anova(hm0.me3, hm1.me3)

mod.street.dbhdelta.me <- hm2.me3
save(mod.street.dbhdelta.me, file="processed/mod.street.dbhdelta.me.sav")
save(mod.street.dbhdelta.me, file="processed/mod.street.final.sav") ## for use in results writeup
load("processed/mod.street.dbhdelta.me.sav")
plot(mod.street.dbhdelta.me)
hist(residuals(mod.street.dbhdelta.me)) ## looks good enough to me
#####


##### PART 2: STREET TREE SIMULATOR
######
### Simulator version history
## V1: Big trees processed separately from small trees (<20k kg), using a modified sampling distribution from street.dbh
## V2: Sampling window on DBH was moving/shrinking as a function of cell biomass (never fully processed)
## V3: Sample weighting of dbh distribution was adjusted towards large end if simulations failed
## V4: tolerance on matching biomass distribution was adjusted to a static threshold (addresses consistent undershoot in large biomass cells)
## V5: urban-specific allometrics on the front end
## V6: uses the correct canopy 1m map aggregated to 30m, based on Raciti's biomass>0 coverage at 1m (NOT dataverse boston canopy tif)

# can <- raster("processed/boston/bos.can.redux30m.tif")
# can.dat <- as.data.table(as.data.frame(can))
# badcan <- raster("processed/boston/bos.can30m.tif")
# badcan.dat <- as.data.table(as.data.frame(badcan))
# plot(badcan.dat$bos.can30m, can.dat$bos.can.redux30m, pch=15, cex=0.6)
# abline(a=0, b=1, col="red") ## almost always the biomass>0 canopy is less than the badcan coverage at 30m pixel scale

###### Approach 5: Reconfigure street tree analysis for more effective search
## Look at cell biomass sample smartly from dbh data, bounded to stop excessive density or BA to reach the biomass total
## read in biomass and canopy data
biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m kg-biomass to 30m pixel
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- crop(biom, aoi) ## biomass was slightly buffered, need to clip to match canopy fraction raster
biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[,aoi:=as.vector(getValues(aoi))]
can <- raster("processed/boston/bos.can.redux30m.tif")
can.dat <- as.data.table(as.data.frame(can))
biom.dat <- cbind(biom.dat, can.dat)
biom.dat[,index:=1:dim(biom.dat)[1]] ## master pixel index for stack of biom/aoi/can, 354068 pix
names(biom.dat)[3] <- "bos.can30m"
# p <- ecdf(biom.dat[!is.na(bos.biom30m) & bos.biom30m>10,bos.biom30m]) ## cdf of cell biomass, vast bulk of # is below 20k
# j <- ecdf(biom.dat[!is.na(bos.biom30m) & bos.biom30m>10,bos.can30m])

## prep street tree data and biomass data for processing
street <- as.data.table(read.csv("processed/boston/street.trees.dbh.csv"))
street.allo <- read.csv("docs/street.biometrics.csv")
street[, biom.2006.urb:=street.allo[match(street[,genus], street.allo$genus, nomatch=8), "b0"]*(street[,dbh.2006]^street.allo[match(street[,genus], street.allo$genus, nomatch=8), "b1"])*street.allo[match(street[,genus], street.allo$genus, nomatch=8), "dens"]]

ba.pred <- function(x){(x/2)^2*pi*0.0001} ## get BA per stump from dbh
clean <- street[record.good==1,] # get a good street tree set ready
clean[, ba:=ba.pred(dbh.2006)]
setkey(clean, biom.2006.urb) # 2592 records in final selection, dbh range 5-112 cm
clean <- clean[order(clean$dbh.2006, decreasing = F),]
clean[,rank:=seq(from=1, to=dim(clean)[1])] ## all the trees have a fixed size rank now (used in adjusting sampling weights)

### set up data
runme <- biom.dat[!is.na(bos.biom30m) & bos.biom30m>10 & !is.na(aoi) & !is.na(bos.can30m) & bos.can30m>0,] ## 108k

### To parallelize this process (for the lower biomass pixels), the script checks for already-written chunks of results, and then tries to produce the next chunk
### calling the script multiple times will result in multiple successive chunks of pixels being run at the same time

## parallel process: check to see if any results have been written to disk, if not queue up the next chunk (or fill in gaps)
## set parameters to divide the file up
chunk.size=2000 ## how many pixel to handle per job ## 10000 is probably too big, if you do this again go for smaller chunks
file.num=ceiling(dim(runme)[1]/chunk.size)
pieces.list <- seq(chunk.size, by=chunk.size, length.out=file.num) ## how the subfiles will be processed
vers <- 6 ## set version for different model runs here

## check existing npp files, find next file to write
check <- list.files("processed/boston/biom_street")
check <- check[grep(check, pattern=paste("ann.npp.street.v", vers, sep=""))] ### version label here
already <- sub(".*weighted\\.", "", check)
already <- as.numeric(sub("\\..*", "", already))
notyet <- pieces.list[!(pieces.list%in%already)]

## if any chunks are not fully processed, next check to see if they're being worked on currently
## the one weakness of this (besides requiring manual launch of each chunk)
## is that if a script times out or otherwise aborts before successfully completing after the step below
## it can't release the unfinished job from being "in process" and will keep the other scripts from running that chunk
if(length(notyet)!=0){ 
  if(!file.exists(paste("processed/boston/biom_street/atwork", vers, "csv", sep="."))){ ## if it isn't there start a file of who is working now
    l <- data.frame(at.work=integer())
    write.csv(l, file=paste("processed/boston/biom_street/atwork", vers, "csv", sep="."))
  }
  atwork <- read.csv(paste("processed/boston/biom_street/atwork", vers, "csv", sep="."))
  atwork <- atwork$at.work
  ## set target for what isn't completed and isn't being worked on
  y <- as.numeric(min(notyet[!(notyet%in%atwork)])) 
  
  if(length(y)==0 | !is.finite(y)){stop("all pixels currently finished or in process")}
  print(paste("going to work on chunk", y)) 
  ## update the at.work file to warn other instances
  atwork <- c(atwork, y)
  l <- data.frame(at.work=atwork)
  write.csv(l, file=paste("processed/boston/biom_street/atwork", vers, "csv", sep="."))
  
  runme.x <- runme[(y-(chunk.size-1)):y,] ## que up the next chunk
  if(y>(dim(runme)[1])){ ## if you're at the end of the file, only grab up to the last row
    runme.x <- runme[(y-(chunk.size-1)):dim(runme)[1],]
  }
}else{stop("all pixels already processed")}

## set up containers
cage.num.trees <- list() ## number of trees per cell, record every successful simulation (vector of 100 per cell)
cage.ann.npp <- list() ## aggregated npp per cell, record every successful simulation (vector of 100 per cell)
cage.dbh <- list() ## the dbh of individual trees chosen in each successful simulation (a vector for each simulation, 100 vectors per cell)
cage.genus <- list() ## the genus of individual trees chosen in each successful simulation (a vector for each simulation, 100 vectors per cell)
index.track <- rep(999999, dim(runme.x)[1]) ## track of order of unique pixel IDs simulated
biom.track <- rep(999999, dim(runme.x)[1]) ## track of biomass in each pixel
cage.wts <- list() ## the sampling weight used in each successful simulation (vector of 100 per cell)
cage.biom.sim <- list() ## the aggregate biomass estimated in each successful simulation (vector of 100 per cell)
attempts.track <- rep(999999, dim(runme.x)[1]) ## the total number of simulations made (success+fail) for each cell
proc.track <- rep(999999, dim(runme.x)[1]) ## the exit status of each cell (1=success, 0 = fail)
cage.spp <- list() ## keep track of species of each selected tree

## for dynamic weighting of sampling of dbh records
incr=300 ## incrememnt to change weighting, how many failures before reaching max weighting shift
k=200 ## constant, arbitrary top of prob weights, higher-->steeper change in sampling weight
D=k/incr ## drop increment, sensitive to number of tries to make while adjusting the sampling weights

## loop each row = pixel
for(t in 1:dim(runme.x)[1]){
  ann.npp <- numeric()
  num.trees <- integer()
  biom.sim.track <- numeric()
  wts.track <- integer()
  cage.genus[[t]] <- list() ### if you can't get a succseful sim, the list entry for that pixel is empty
  cage.dbh[[t]] <- list()
  cage.spp[[t]] <- list()
  
  x <- 0 ## count successful simulations
  q <- 0 ## count attempts since last successful sample
  g <- 0 ## weight adjustment start
  z <- 0 ## total number of sample iterations
  jeez <- 0 ## if it's grinding through a slow one
  
  ### approach below uses moving/shrinking window on the street tree records as biomass starts to get bigger
  ## has the disadvantage of supressing the number of trees at higher biomass (start to oversample very large trees), would need some tuning to get the right sampling window
  # cell.q <- p(runme[t, bos.biom30m]) ## this is the percentile of the cell biomass viz. rest of data
  # spread <- min((1-cell.q), 0.5) ## interval to sample in dbh records, cap it off if the spread gets bigger than 50%
  # window.l <- max(0, (cell.q-spread))
  # window.h <- cell.q+spread
  # if((spread)==0.5){window.h <- 1} ## basic scheme is that for biomass>median, restrict to (rarer) higher window, but if biom<median, sample whole window
  # dbh.lims <- quantile(clean$dbh.2006, probs=c(window.l, window.h)) ## window of dbh to select based on cell biomass
  # grasp <- clean[biom.2006<runme[t, bos.biom30m] &dbh.2006<=dbh.lims[2] & dbh.2006>=dbh.lims[1], .(dbh.2006, biom.2006, ba)] ## this is a street tree record restricted based on cell biomass
  
  dens.lim <- ceiling(0.106*runme.x[t,aoi*bos.can30m]) ## maximum number to select, probably only restrictive in low-canopy cells
  tol <- 100 ## maximum tolerance for simulation deviation from cell biomass, in kg
  tol.window <- c(max(c((runme.x[t, bos.biom30m]-tol), (0.9*runme.x[t, bos.biom30m]))), 
                  min(c((runme.x[t, bos.biom30m]+tol), (1.1*runme.x[t, bos.biom30m]))))
  ### or take a "groomed" record and sample it with adjusting weights if the algorithm fails on density
  grasp <- clean[biom.2006.urb<(tol.window[2]), .(dbh.2006, biom.2006.urb, ba, rank, genus, Species)] ## exclude sampling in this simulation any single trees too big to fit cell biomass
  
  ### sample street trees and try to fill the cells
  while(x<100 & q<1000){ ## select 100 workable samples, or quit after q attempts with no success
    # ## grab a sample of grasp (moving/shrinking, no weighting)
    # samp <- grasp[sample(dim(grasp)[1], size=dens.lim, replace=T),] ## sample grasp with replacement only up to density limit
    
    ## OR: figure out the dynamic weighting to use based on previous failures
    wts=(k-(g*D))+(((g*D)/(dim(grasp)[1]-1))*((grasp$rank)-1))
    # sample grasp with replacement only up to density limit, with specified weights based on iterative reweighting
    samp <- grasp[sample(dim(grasp)[1], size=dens.lim, replace=F, prob = wts),] 
    w=samp[1, biom.2006.urb] ## keep cummulative tally of biomass
    d=1 # keep track of the number of trees
    while(w<tol.window[1] & d<dens.lim){ ## keep adding trees until you just get over the target biomass or run out of records
      w=w+samp[d+1, biom.2006.urb]
      d=d+1
    }
    ### if this is too many trees too tightly packed (over 40m2/ha BA) or not enough biomass in the sample, readjust weights
    if(((samp[1:d, sum(ba)]/runme.x[t, aoi*bos.can30m])*1E4)>40 | w<tol.window[1]){ ## we are counting area of "forest" as the area covered in canopy (NOT pervious, NOT raw ground area)
      if(g<incr){
        g=g+1 ## readjust the sample weights
        if(g==300){print(paste("pixel", runme.x[t, index], "weights at maximum"))}
        q=0 ## reset attempt timeout clock
      }
    }
    if(((samp[1:d, sum(ba)]/runme.x[t, aoi*bos.can30m])*1E4)<40 & w<tol.window[2] & w>tol.window[1]){ ## if the BA density is low enough & got the biomass simulated properly
      x <- x+1 ## record successful sample
#       ann.npp <- c(ann.npp, sum(samp[1:d, biom.2006.urb]*exp((mod.biom.rel$coefficients[2]*log(samp[1:d, dbh.2006]))+mod.biom.rel$coefficients[1])))
      ann.npp <- c(ann.npp, 99999) ## don't calculate NPP in the simulator round
      num.trees <- c(num.trees, d)
      biom.sim.track <- c(biom.sim.track, w) ## simulated biomass in this sample
      cage.dbh[[t]][[x]] <- samp[1:d, dbh.2006] ## which dbhs did you select
      cage.genus[[t]][[x]] <- samp[1:d, genus] # running list of which genera you select
      cage.spp[[t]][[x]] <- samp[1:d, Species]
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
  if(x==0){
    ann.npp <- NA 
    num.trees <- NA
    biom.sim.track <- NA
    wts.track <- NA
    }
  ## store results, record exit status
  cage.ann.npp[[t]] <- ann.npp
  cage.num.trees[[t]] <- num.trees
  cage.wts[[t]] <- wts.track
  biom.track[t] <- runme.x[t, bos.biom30m]
  index.track[t] <- runme.x[t, index]
  cage.biom.sim[[t]] <- biom.sim.track
  attempts.track[t] <- z
  ## exit status
  if(q<1000){
    proc.track[t] <- 1
    print(paste("finished pixel", runme.x[t, index], "=", t))
  } else{
    proc.track[t] <- 0 ## record as incomplete processing
    print(paste("pixel", runme.x[t, index], "failed =", t))
  }
}

## when complete dump everything back into the save file
save(cage.ann.npp, file=paste("processed/boston/biom_street/ann.npp.street.v", vers, ".weighted.", y, ".sav", sep=""))
save(cage.num.trees, file=paste("processed/boston/biom_street/num.trees.street.v", vers, ".weighted.", y, ".sav", sep=""))
save(cage.dbh, file=paste("processed/boston/biom_street/dbh.street.v", vers, ".weighted.", y, ".sav", sep=""))
save(cage.genus, file=paste("processed/boston/biom_street/genus.street.v", vers, ".weighted.", y, ".sav", sep=""))
save(biom.track, file=paste("processed/boston/biom_street/biom.track.street.v", vers, ".weighted.", y, ".sav", sep=""))
save(index.track, file=paste("processed/boston/biom_street/index.track.street.v", vers, ".weighted.", y, ".sav", sep=""))
save(proc.track, file=paste("processed/boston/biom_street/proc.track.street.v", vers, ".weighted.", y, ".sav", sep=""))
save(cage.biom.sim, file=paste("processed/boston/biom_street/biom.sim.street.v", vers, ".weighted.", y, ".sav", sep=""))
save(attempts.track, file=paste("processed/boston/biom_street/attempts.track.street.v", vers, ".weighted.", y, ".sav", sep=""))
save(cage.wts, file=paste("processed/boston/biom_street/wts.street.v", vers, ".weighted.", y, ".sav", sep=""))

## update the at.work file to release this job
check <- list.files("processed/boston/biom_street")
check <- check[grep(check, pattern=paste("ann.npp.street.v", vers, sep=""))] ### version label here
already <- sub(".*weighted\\.", "", check)
already <- as.numeric(sub("\\..*", "", already))
notyet <- pieces.list[!(pieces.list%in%already)]

## update atwork to remove any in-progress files that now have finished files on disk
atwork <- read.csv(paste("processed/boston/biom_street/atwork", vers, "csv", sep="."))
atwork <- atwork$at.work
atwork <- atwork[atwork%in%notyet] ## clear job records that have completed files on disk
l <- data.frame(at.work=atwork)
write.csv(l, file=paste("processed/boston/biom_street/atwork", vers, "csv", sep="."))
######


##### PART 3: RECONSTRUCTION OF SIMULATOR RESULTS
######
### version 1 (not labeled) used static biomass growth equation
### version 2 (V42) resimmed each pixel x100 using the mixed effects dbh increment model and Jenkins/Chognacky allometrics to estimate npp
### version 3 (V43) resimmed each pixel x100 using the mixed effect dbh increment and specific allometrics for dbh-->volume-->biomass
### version 4 (v5) resimmed x100 with urban specific allometrics based on the V5 biomass simulator that also used urban-specific allometrics (instead of default hardwood) to estimate pixel biomass
### Version 5? (v6) 100x resims per pixel using the biomass>0 canopy raster for the 30m canopy fraction

### urban-specific allometrics
## McPherson, E.G., N.S. van Doorn and P.J. Peper. 2016. Urban tree database and allometric equations. Gen. Tech. Rep. PSW-GTR-253. Albany, CA: U.S. Department of Agriculture, Forest Service, Pacific Southwest Research Station. 86 p
street.allo <- read.csv("docs/street.biometrics.csv") ## AG wood vol (m3) as b0*DBH(cm)^b1; multiply by density to get kg-biomass

## to get basal area and biomass of stems --- static biomass growth prediction
# ba <- function(x){(((((x/2)^2)*pi)/900))} ## this is in m2/ha (correct for 900m2 pixel footprint)
ba <- function(x){(((x/2)^2)*pi)/1E4} ### to find JUST the BA of a tree based on dbh (in m2)

#### THIS CODE CHUNK does the growth calculations using the saved object files and exports and R objects. 
#### THE NEXT CODE CHUNK reloads the R objects and reconstructs a map file with all the different iterations of NPP
### for each pixel, get median estimated npp, median # trees
container <- data.frame()
med.dbh.rec <- numeric() ## keep a running tally of the dbh of every tree in the nearest-to-median-npp sample in each pixel
dbh.dump <- list() ## a place to put actual dbh samples from targeted retrievals for deeper analysis (e.g. what do the median retrievals look like?)
load("processed/mod.street.dbhdelta.me.sav") ## polynomial mixed effect with species

### line up a list of randomly selected coefficients for the operatie dbh~growth model
## mod Feb 11 2019 -- make it 1000 long
b0.rand <- rnorm(1000, mean=coef(summary(mod.street.dbhdelta.me))[1,1], sd=coef(summary(mod.street.dbhdelta.me))[1,2])
b1.rand <- rnorm(1000, mean=coef(summary(mod.street.dbhdelta.me))[2,1], sd=coef(summary(mod.street.dbhdelta.me))[2,2])
b2.rand <-rnorm(1000, mean=coef(summary(mod.street.dbhdelta.me))[3,1], sd=coef(summary(mod.street.dbhdelta.me))[3,2])

# x <- 1:100
# plot(x, mean(b0.rand)+(mean(b1.rand)*x)+(mean(b2.rand)*(x^2))) ## looks like a believable growth curve

#### import results objects pulled from parallel processing on the cluster (chunks of 10k pixels)
vers <- 6 ## which simulator model run are we picking at
obj.dump <- list.files("processed/boston/biom_street/")
npp.dump <- obj.dump[grep(obj.dump, pattern = paste("ann.npp.street.v", vers, ".weighted*", sep=""))] ## use npp empty file as a way of organizing the work
npp.dump.chunks <- sub('.*weighted.', '', npp.dump)
npp.dump.chunks <- sub("\\.sav.*", "", npp.dump.chunks) ## later version get named .sav

##### process each pixel chunk; for each pixel estimate npp in a randomly selected dbh collection 100x using b0/b1/b2 values from the vectors above
npp.random <- list() ## container for npp simulations, will contain pixel-n elements per chunk with ~100 values in each
i=1
for(c in 1:length(npp.dump)){
  load(paste("processed/boston/biom_street/", npp.dump[c], sep=""))
  tmp.npp <- sapply(cage.ann.npp, FUN=median)
  tmp.num.sims <- sapply(cage.ann.npp, FUN=length) ## number of successful simulations recorded
  tmp.sim.incomp <- rep(0, length(tmp.num.sims))
  tmp.sim.incomp[tmp.num.sims<100] <- 1 ## vector tracking the incomplete simulations

  ## grab the corresponding other dump files
  load(paste("processed/boston/biom_street/num.trees.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "num.trees" object
  tmp.num <- as.integer(sapply(cage.num.trees, FUN=median))
  tmp.num[tmp.num==9999] <- NA ## filter the NA flags
  load(paste("processed/boston/biom_street/index.track.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "index.track" object
  tmp.index <- index.track
  load(paste("processed/boston/biom_street/biom.track.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "biom.track" object
  tmp.biom <- biom.track
  load(paste("processed/boston/biom_street/proc.track.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "proc.track" object
  tmp.proc <- proc.track
  load(paste("processed/boston/biom_street/biom.sim.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "proc.track" object
  tmp.biom.sim <- sapply(cage.biom.sim, FUN=median)
  load(paste("processed/boston/biom_street/attempts.track.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "proc.track" object
  tmp.attempts <- attempts.track
  load(paste("processed/boston/biom_street/wts.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "cage.wts" object
  tmp.wts <- sapply(cage.wts, FUN=max)
  load(paste("processed/boston/biom_street/genus.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "cage.wts" object

  ### get the dbh collection for every simulation in every pixel
  load(paste("processed/boston/biom_street/dbh.street.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep="")) ## comes in as "dbh.stret.small" object
  ba.grand <- rep(9999, length(cage.dbh)) ## BA values for every retreival per cell
  dbh.grand <- rep(9999, length(cage.dbh)) ## dbh values for the retreival nearest median NPP
  print(paste("patience! Doing pixel-level calculations on chunk", npp.dump.chunks[c]))
  for(b in 1:length(cage.dbh)){ ## b is an individual pixel, with 100 separately simulated dbh samples
    if(length(cage.ann.npp[[b]])<5){  ## if too few NPP retreivals were made for this pixel
      ba.grand[b] <- NA
      dbh.grand[b] <- NA
      dbh.dump[[b]] <- NA
    } else{
      ba.grand[b] <- median(sapply(sapply(cage.dbh[[b]], FUN=ba), FUN=sum)) ## median of the total ba (=m2) for each sample, need to convert based on (canopy) area of pixel
      dbh.grand[b] <- median(unlist(cage.dbh[[b]])) ## grand median dbh of all trees selected for all samples
      ### these look at the summary stats for a particular retreival in each cell (ex. the retrieval nearest the median npp)
      dev <- abs(cage.ann.npp[[b]]-tmp.npp[b]) ## deviance of individual NPP retreivals from the median NPP for this cell
      rrr <- which(dev==min(dev)) ## which retreival is the closest to the median
      dbh.dump[[b]] <- cage.dbh[[b]][[rrr[1]]] ## track the tree sample nearest to every npp median in every pixel (just take the first instance, fuck it)

      #### new hotness: throw varying model parameters at a random dbh collection
      #### also newer hotness: match dbh samples via associated genus to urban-specific volumetric growth equations
      ## randomly select a dbh collection, and successively apply all 1000 combinations of dbh growth coefficients
      npp.random[[b]] <- numeric() ## this is where we are storing the vectors (100 long) of npp estimates for each pixel
      for(p in 1:1000){  ## individual entries in cage.dbh are pixels, which contain vectors of dbh collections
          dbh.rand <- sample(1:length(cage.dbh[[b]]), size=1) ## randomly select a dbh collection
          tmp.dbh0 <- cage.dbh[[b]][[dbh.rand]] ### dbh collection at time 0
          tmp.genus <- cage.genus[[b]][[dbh.rand]] ## genera of dbh colletion

          ### fast match the genus to the specific allometric if possible, or divert to general equation at row 8
          tmp.biom0 <- street.allo[match(tmp.genus, street.allo$genus, nomatch=8), "b0"]*(tmp.dbh0^street.allo[match(tmp.genus, street.allo$genus, nomatch=8), "b1"])*street.allo[match(tmp.genus, street.allo$genus, nomatch=8), "dens"]
          tmp.dbh1 <- tmp.dbh0+(b0.rand[p]+(b1.rand[p]*tmp.dbh0)+(b2.rand[p]*(tmp.dbh0^2))) ## grow the dbh to time 1
          tmp.biom1 <- street.allo[match(tmp.genus, street.allo$genus, nomatch=8), "b0"]*(tmp.dbh1^street.allo[match(tmp.genus, street.allo$genus, nomatch=8), "b1"])*street.allo[match(tmp.genus, street.allo$genus, nomatch=8), "dens"]
          npp.random[[b]] <- c(npp.random[[b]], sum(tmp.biom1-tmp.biom0)) ## store the change in biomass as an npp sample
        } ## end npp random estimator
      } ## end else check for number of retrievals
    if(b%%1000==0){print(b)}
      i=i+1 ## keep track of number of pixels processed
    } ## end pixel loop n=b

  ### save out the chunk of pixels with the results of all random simulations
  save(npp.random, file=paste("processed/boston/biom_street/results/npp.random.v", vers, ".weighted.", npp.dump.chunks[c], ".sav", sep=""))

  ## bind (static) results in container
  container <- rbind(container,
                     cbind(tmp.index, tmp.biom, ## pixel index and map biomass 
                           tmp.num, dbh.grand, ba.grand, tmp.biom.sim, ## median npp and tree number for simulator run
                           tmp.wts, tmp.num.sims, tmp.sim.incomp, tmp.attempts, tmp.proc, ### metrics for how well the simulator performed
                           median(npp.random[[b]])) ## median npp for pixel w randomly selected dbh + model error
  )
} ## end chunk loop

## Collect simulation diagnostics and export
names(container) <- c("pix.ID", "biom.kg",
                       "med.tree.num.sim", "med.dbh.sim", "med.ba.sim", "med.biom.sim", ### simulator medians
                      "max.wts", "num.sims", "sim.incomp", "attempts", "proc.status", ## simulator health indicators
                      "med.npp.rand")
write.csv(container, paste("processed/boston/biom_street/results/streettrees.sim.v", vers, ".diagnostics.csv", sep=""))

hist(container$biom.kg) ## just the cell biomass of the pixels simulated
hist(container$med.tree.num.sim); summary(container$med.tree.num.sim) # med 6 per pixel, up to 38 (slight decrease from using the earlier version with the more generous canopy -- there's been some stricture on putting in trees)
hist(container$med.dbh.sim); summary(container$med.dbh.sim) ## the median dbh of every sim, most sims have median dbh of about 25, range from below 20 to up to 40
hist(container$med.ba.sim) ## UP TO 3.5, @900m2/100% canopy = 39 m2/ha -- good so we don't get any higher than andy's edge plots
hist(container$max.wts) ## most get it done before 50 wt, but a substantial minority are maxing out
hist(container$num.sims) ## most get all 100 sims; a minority get less than 20, nothing gets 20-100
table(container$sim.incomp) ## 98k get a simulation; 9.6k get partial/failed
hist(container$attempts); summary(container$attempts) ## most only have to try a few hundred times, but a few go up to 400k times
hist(container$med.npp.rand); summary(container$med.npp.rand) ## median about 100 kg/pix, max is lower at 538 kg (previous sim maxes out at 685 kg/pix)


gawlie <- data.frame(cbind(container$pix.ID, container$num.sims, container$sim.incomp))
colnames(gawlie) <- c("pix.ID", "num.sims", "sim.incomp")
table(gawlie$sim.incomp) ## v6 is 98264 succeeded, 9622 failed (v5 was 5831 failed, 102055 worked)
write.csv(gawlie, paste("processed/boston/biom_street/results/streettrees.sim.v", vers, ".status.csv", sep=""))


###
## prepare data files with each random pixel realization as a column (within-pixel spread is across rows)
vers=6
npp.random.list <- list.files("processed/boston/biom_street/results")
npp.random.list <- npp.random.list[grep(npp.random.list, pattern = paste("npp.random.v", vers, ".weighted.*", sep=""))] ## V5 is 1)mixed model for DBH change; 2) urban-specific allometrics 3) consistent coefficients across each map realization;  4) uses biomass simulator V5 that also used urban-specific allometrics
npp.random.chunks <- sub('.*weighted.', '', npp.random.list)
npp.random.chunks <- sub("\\.sav.*", "", npp.random.chunks) ## later version get named .sav

sim.status <- as.data.table(read.csv(paste("processed/boston/biom_street/results/streettrees.sim.v", vers, ".status.csv", sep="")))

### associate each set of npp retreivals to its pixel ID and write to disk
for(c in 1:length(npp.random.chunks)){
  load(paste("processed/boston/biom_street/results/", npp.random.list[c], sep="")) ## npp.random
  load(paste("processed/boston/biom_street/index.track.street.v", vers, ".weighted.", npp.random.chunks[c], ".sav", sep="")) ##  "index.track"
  tmp.npp <- index.track
  tmp.npp <- cbind(tmp.npp, matrix(nrow=length(index.track), ncol=1000, NA))
  for(i in 1:length(index.track)){ ## row-wise reconstruction
    if(sim.status[pix.ID==index.track[i],num.sims]<100){ ## for incomplete simulations
      tmp.npp[i,2:1001] <- c(npp.random[[i]], rep(NA, 1000-length(npp.random[[i]]))) ## fill in any remaining holes in the npp retreival
    }else{
      tmp.npp[i,2:1001] <- npp.random[[i]]
    }
    # print(paste("pixel", i))
  }
  write.csv(tmp.npp, file = paste("processed/boston/biom_street/results/street.npp.random.v", vers, ".estimates.", npp.random.chunks[c], ".csv", sep=""))
  print(paste("just wrote chunk", npp.random.chunks[c], "npp results to disk"))
}

### load up results matrices and bind to single matrix
npp.results.list <- list.files("processed/boston/biom_street/results/")
npp.results.list <- npp.results.list[grep(pattern = paste0("street.npp.random.v", vers), npp.results.list)]
gosh <- data.frame()
for(c in 1:length(npp.results.list)){
  tmp <- read.csv(paste0("processed/boston/biom_street/results/", npp.results.list[c]))
  gosh <- rbind(gosh, tmp)
  print(paste('chunk', c, "loaded"))
}

## merge npp simulation results in with map pixel data
library(data.table)
library(raster)
biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m kg-biomass to 30m pixel
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- crop(biom, aoi) ## biomass was slightly buffered, need to clip to match canopy fraction raster
biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[,aoi:=as.vector(getValues(aoi))]
can <- raster("processed/boston/bos.can.redux30m.tif")
can <- crop(can, aoi)
biom.dat[,can:=getValues(can)]
isa <- raster("processed/boston/bos.isa30m.tif")
isa <- crop(isa, aoi)
biom.dat[,isa:=getValues(isa)]
lulc <- raster("processed/boston/bos.lulc30m.lumped.tif")
lulc <- crop(lulc, aoi)
biom.dat[,lulc:=getValues(lulc)]
biom.dat[,pix.ID:=1:dim(biom.dat)[1]]
names(biom.dat)[1] <- "biom"

gosh <- gosh[,-1]
names(gosh) <- c("pix.ID", paste("npp.street.random.iter.", 1:1000, ".kg", sep=""))
map <- merge(x=biom.dat, y=gosh, by="pix.ID", all.x=T, all.y=T)
write.csv(map, paste0("processed/results/street/streettrees.npp.simulator.v", vers, ".results.random.csv"))

### Exploratory of what's in the random street tree results iterations
npp.street <- as.data.table(read.csv("processed/results/street/streettrees.npp.simulator.v6.results.random.csv"))
sum.na <- function(x){sum(x, na.rm=T)}
street.tot <- apply(npp.street[, 8:1007], MARGIN=2, FUN=sum.na)
hist(street.tot/2000)
median(street.tot/2000) ## 9.9 ktC for whole-map sum (missing some values due to more things not being simmed; was in v5 11.1k tC across all map
nonfor.tot <- apply(npp.street[aoi>800 & lulc!=1 & biom<20000, 8:1007], MARGIN=2, FUN=sum.na)
summary(nonfor.tot/2000); hist(nonfor.tot/2000) ## about 7.5 ktC successfully simulated for the nonforest/nonbig pixels
for.tot <- apply(npp.street[aoi>800 & lulc==1 | aoi>800 & biom>=20000, 8:1007], MARGIN=2, FUN=sum.na) 
summary(for.tot/2000)### about 2.3 simulated in forest/large-biomass pixels

## how many sims do we have in the different pixel classes?
dim(npp.street[lulc!=1 & aoi>800,]) ## 125k
dim(npp.street[lulc==1 & aoi>800,]) ## 11.3k
na.detect <- function(x){sum(is.na(x))}
nonfor.sims <- apply(npp.street[aoi>800 & lulc!=1, 8:1007], MARGIN=1, FUN=na.detect)
table(nonfor.sims) ### 95352 nonforest pix have 1000 successful retrievals, 29759 are not retrieved
for.sims <- apply(npp.street[aoi>800 & lulc==1, 8:1007], MARGIN=1, FUN=na.detect)
table(for.sims) ## 11248 forest pixels with no sim values, 22 that have em (???)
bigbiom.sims <- apply(npp.street[lulc!=1 & biom>=20000, 8:1007], MARGIN=2, FUN=na.detect)
dim(npp.street[lulc!=1 & biom>=20000 & aoi>800, 8:1007]) ### 2.7k of these
View(npp.street[lulc!=1 & biom>=20000 & aoi>800, 1:20]) ## OK so there seem to be large-biomass simulations aplenty
View(npp.street[lulc==1 & aoi>800, 1:20])
dim(npp.street[lulc==1 & aoi>800, 1:20]) ## 11.3k forest pixels a lot of which seem to have simulation runs recorded

quantile(nonfor.tot/2000, probs=c(0.05, 0.5, 0.95)) ## ##v6 4.4, 7.5 10.8k tC for simulated street trees, old v5 was 4.8-12.2k tC for simulated street trees



### export the tifs
npp.street <- fread("processed/results/street/streettrees.npp.simulator.v6.results.random.csv")
median.na <- function(x){median(x, na.rm=T)}
med.street <- apply(as.matrix(npp.street[,8:1007]), MARGIN=1, FUN=median.na)
r <- biom
r <- setValues(r, med.street)
writeRaster(r, "processed/results/street/bos.street.V6.npp.med.tif", format="GTiff", overwrite=T)

# for(j in 1:4){
#   r <- biom
#   r <- setValues(r, map[[(j+5)]])
#   writeRaster(r, paste("processed/boston/", names(map)[5+j], ".v1.tif", sep=""),
#               format="GTiff", overwrite=T)
# }


# ## brief exploratory
# ### do simulations track cell biomass well enough?
# vers=4
# container <- read.csv(paste("processed/streettrees.npp.simulator.v", vers, "3.results.random.csv", sep=""))
# container <- as.data.table(container)
# 
# ### how well did the median biomass simulation per cell get to measured biomass?
# summary(container$bos.biom30m)
# summary(container$biom.kg)
# summary(container$med.biom.all)
# plot(container$biom.kg, container$med.biom.all) ## mostly in a tight range -- no successful sims over ~35k
# abline(a=0, b=1) ## the simulated biomass is very near the cell biomass
# container[,biom.dev:=med.biom.all-biom.kg]
# plot(container$biom.kg, container$biom.dev) ## most within +/- 50kg till right at 30k, expands to 100kg
# 
# ### how prevalent are simulation failures, where are they
# sum(container$proc.status==0, na.rm=T)/dim(container[proc.status%in%c(0,1),])[1] # 3.46% failure to simulate
# sum(container$biom.kg>30000, na.rm=T)/dim(container[is.finite(bos.biom30m),])[1] # 1.35% are >30k -- so we are getting some but not all of the largest 1%
# container[biom.kg>25000 & proc.status==0, length(biom.kg)]/container[biom.kg>25000, length(biom.kg)] ## 66% of pixels >25k failed
# container[biom.kg>25000 & proc.status==0, length(biom.kg)]/container[proc.status==0, length(biom.kg)] ## 73% of failures are over 25k
# # map[proc.status==0, median(bos.can30m, na.rm=T)] ## median 99.8% canopy
# # hist(map[proc.status==0, (bos.can30m)]) ## most heavily covered except for a small number that must be low biomass
# # hist(map[proc.status==0, bos.biom30m]) ### bimodal, very low biomass or very high biomass
# ### might need to just toss anything above ~20k kg -- call it "much more like a forest than a piece of city"
# 
# ### how is productivity organized in space
# plot(container$biom.kg, container$med.ann.npp.all) # tent shape, nearly linear increase up to ~25k then steep decline (sampling large old trees?)
# abline(v=22000) ## peaks right at 22k then steep decline
# # plot(map$bos.can30m, map$med.ann.npp.all) ## everything increases but with inreasing variance up to full canopy where it can be all over
# # hist(map[bos.biom30m>10, med.ann.npp.all])
# # hist(map[bos.biom30m>10 & bos.biom30m<22000, med.ann.npp.all]) ## how different are our results if we chuck the weird ass high biomass pixels?
# 
# ## NPP rate, city wide
# container[aoi>800 & bos.biom30m>10, npp.MgC.ha:=(((med.ann.npp.all/1000)*(1/2))/aoi)*1E4]
# container[aoi>800, range(npp.MgC.ha, na.rm=T)] ## 0.04-5.2 MgC/ha/yr
# hist(container[aoi>800, npp.MgC.ha])
# container[bos.biom30m>22000, length(bos.biom30m)]/container[bos.biom30m>10, length(bos.biom30m)] ##5.8% of pix are above 22k
# container[is.finite(med.ann.npp.all) & aoi>800 & bos.biom30m>22000, length(med.ann.npp.all)]/container[aoi>800 & bos.biom30m>22000, length(med.ann.npp.all)] ## we got a value for 66% of the large pixels
# container[aoi>800, sum(med.ann.npp.all/(2*1000), na.rm=T)] #13161 MgC/yr
# container[aoi>800, sum(med.ann.npp.all/(2*1000), na.rm=T)]/(container[, sum(aoi, na.rm=T)]/1E4) ### 1.06 MgC/ha/yr
# ### how many large pixels did not simulate?
# container[aoi>800 & bos.biom30m>20000 & proc.status==0, length(bos.biom30m)] #2.8k
# container[aoi>800 & bos.biom30m>20000 & proc.status==0, length(bos.biom30m)]/container[aoi>800 & is.finite(bos.biom30m), length(bos.biom30m)]
# 
# ### excluding large pixels
# container[aoi>800 & bos.biom30m<22000, sum(med.ann.npp.all/(2*1000), na.rm=T)] #11789 MgC/yr
# container[aoi>800 & bos.biom30m<22000, sum(med.ann.npp.all/(2*1000), na.rm=T)]/(container[, sum(aoi, na.rm=T)]/1E4) ### 0.95 MgC/ha/yr
# 
# ## average small and large pixel productivity
# container[aoi>800 & bos.biom30m<22000, median(npp.MgC.ha, na.rm=T)] ## 0.94 MgC/ha/yr for small pixels
# container[aoi>800 & bos.biom30m>22000, median(npp.MgC.ha, na.rm=T)] ## ### 3.65 MgC/ha/yr
# ## contrast to Andy-forest-based estimate, #13.8k tC, 1.1 tC/ha/yr mean estimate (range was ~7k-22k)
######

