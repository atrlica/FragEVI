library(raster)
library(data.table)
library(ggplot2)


### Part A: Analysis of DBH data from Boston-region FIA plots, provided by Morreale
###########
### Here's the logic of looking at only live stems: 
## 1) The 1m biomass map shows how much is alive at time X, and we are using it to predict NPP up to time Y (1 year)
## 2) The Andy forest calc relies on predicting forest biomass gain *by area* using growth rate~dbh *per stem* for a selection of individual trees
## 2b) The growth rates per stem were determined *from trees that lived* -- no dbh in the forest area was measured for dead wood, no growth mortality rate was estimated for the stand
## 2c) So the Andy growth rate (kg/ha/yr) is estimated as if all biomass survives through the growth period
## 3) Same logic applies to the street tree records: Growth rate *per stem* is only measurable in trees that survived
## 3b) Ian has some model projections on mortality rates, but on an areal basis we only can start with the biomass we have (had, 2006/7) and project forward as if it all survives
## 4) The FIA equation approach is meant apply to a *forest area* directly, and may presumably incorporate mortalities and recruitment along the way
## 5) Alternatively, we can recalculate the areal biomass growth rate as a sum of stem growth~dbh, and determine the growth~dbh relationship from individual stem records
## 5b) To be comparable to results based on the other data sets, the growth~dbh relationship should be estimated based on *the trees that lived* and not on a generalized resample of all trees (as the FIA measures everything, living or dead)


#### Read in FIA data queried by Luca and process 
#####
spp.allo <- read.csv("data/FIA/spp_allometrics.csv")  ## lookup table of allometric equations for the most common representatives in the sample
# live <- read.csv("data/FIA/MA_Tree_Data_ID_NOMORT_SUBID.csv")
# live <- as.data.table(read.csv("data/FIA/MA_Tree_Data_ALLMEASURES.csv"))
live <- as.data.table(read.csv("data/FIA/MA_Tree_Data_STATUS.csv")) ## all stem measurements with status
live[,1] <- NULL
names(live)[3] <- c("TreeID")
names(live)[1] <- c("PlotID")
names(live)[2] <- c("SubID")
spec <- read.csv("data/FIA/REF_SPECIES.csv") ## map of species # to spp
live <- merge(x=live, y=spec[,c("SPCD", "GENUS", "SPECIES")], by.x="SPECIES_CD", by.y="SPCD", all.x=T, all.y=F)
live$GENUS <- as.character(live$GENUS)
live$GENUS <- as.factor(live$GENUS)
live$GENUS.num <- as.numeric(live$GENUS)


### calculate species specific allometrics and biomass changes
live[,spp:=paste(substr(GENUS, 1,1), ".", SPECIES, sep="")]
live <- merge(x=live, y=spp.allo[,c("spp", "b0", "b1")], by="spp", all.x=T)
### if can't find a specific allometric, apply eastern hardwood default
live[is.na(b0), b0:=(-2.48)]
live[is.na(b1), b1:=2.4835]
## calculate biomass at the two measurement periods
biom.pred2 <- function(b0, b1, x){exp(b0+(b1*log(x)))}
live[,biom0.spp:=biom.pred2(b0, b1, DIAM_T0)]
live[,biom1.spp:=biom.pred2(b0, b1, DIAM_T1)]

# Previous approach -- using eastern hardwood defaults for all allometrics
# ### using eastern hardwood defaults, but a lot of the record are Pinus or Tsuga. Might be good to go ahead and use the correct generic equation
# b0 <- -2.48
# b1 <- 2.4835 ## these are eastern hardwood defaults
# biom.pred <- function(x){exp(b0+(b1*log(x)))}
# live$biom0 <- biom.pred(live$DIAM_T0)
# live$biom1 <- biom.pred(live$DIAM_T1)
# live$biom.delt <- live$biom1-live$biom0
# live$biom.rel <- (live$biom.delt/4.8)/live$biom0 ## annualized relative growth increment
# live$dbh.delt <- live$DIAM_T1-live$DIAM_T0
# summary(live$biom.delt) ## some living trees are losing a lot of biomass

# ### how much correction are we seeing by adjusting the equations?
# plot(live[,biom0], live[,biom0.spp])
# abline(a=0, b=1) ## most are slightly more biomass than predicted by E hardwood defaults
# plot(live[, biom1], live[, biom1.spp])
# abline(a=0, b=1) ## same, slightly higher biomass
# summary(live$DIAM_T0) ### most of these aren't too far out of spec for the allometric equations to use

## class as hard or soft wood
live[,type:="H"]
live[spp%in%c("P.strobus", "P.resinosa", "T.canadensis", "A.balsamea"), type:="S"]
live[,type:=as.factor(type)]
live[,biom.delt.spp:=biom1.spp-biom0.spp]

### figure out time lags between measurements of individual stems
### first fix the unique tree ID part -- it can be confused by similar sequence of plot/sub/tree
live[,uniqueID:=paste(PlotID, SubID, TreeID, sep=".")]
live[,UNIQUETREE:=NULL] ## kill the confusing tree ID

### sort out the history of measurement for each stem
arboles <- live[, .(min(YEAR), ## first appearance
                  length(unique(DIAM_T0))), by=uniqueID] ## number of measures
names(arboles) <- c("uniqueID", "first.instance", "num.measures")
arboles[, PlotID:=sub("\\..*","", uniqueID)]
### 8545 first appearance records, total record is 17898, so most have only been measured twice
hist(arboles$first.instance) ## lot of plot establishment in 2003-2006, then a handful of new ones each year
hist(arboles[num.measures==1, first.instance]) ## 2148 trees only measured once, some go back to 2003
summary(live[uniqueID%in%arboles[num.measures==1 & first.instance<2009, uniqueID], DIAM_T0]) ## everything measured '03-08 that disappears has no prior measure, no diam time series
live[uniqueID%in%arboles[num.measures==1, uniqueID] & DIAM_T0>0,] ## 36 stems have only one entry but have prior data (first instance is not recorded here). Not worth bothering, can't tell when first Diam measurement happened
arboles <- arboles[num.measures>1,] ## 6397 individual stems have at least two measurements -- ie. a growth series
live <- live[uniqueID%in%arboles[,uniqueID],] ## 15750 stems with more than one appearance

### figure out the lag in years since the previous measurement
live[, lag:=0] ## initialize, first instance gets lag 0
## figure out the lag time between successive measurements
for(u in 1:dim(arboles)[1]){
  tt <- live[uniqueID==arboles[u, uniqueID],]
  tt <- tt[order(YEAR),]
  tt[2:dim(tt)[1], lag:=diff(tt[,YEAR])]
  for(d in 2:length(tt$YEAR)){
    live[uniqueID==arboles[u, uniqueID] & YEAR==tt[d, YEAR], lag:=tt[d,lag]] ## apply the lags where they belong
  }
  print(paste("getting lag for tree", u, sep=" "))
}

## make sure the lag calculation went well
View(live[lag==0 & order(uniqueID),])
orph <- live[lag==0 & DIAM_T0>0, uniqueID] ## 13 orphans with prior records but no first-instance record
live <- live[!(uniqueID%in%orph), ] ### kill the orphans, down to 15724
live[, length(unique(uniqueID))] ## 6384 stems with DBH histories
length(live[,unique(PlotID)])## 249 unique plots
a <- live[PlotID==231,]
View(a[order(TreeID),])
a[,sum(length(DIAM_T0)), by=YEAR]

### get time series of dbh change
live[lag>0, delta.diam:=DIAM_T1-DIAM_T0] ## 9340 second or more measures
live[lag>0, diam.rate:=delta.diam/lag]
live[lag>0, growth.ann.rel:=(biom.delt.spp/biom0.spp)/lag]

#### Additional grooming: Get rid of records from partially forested subplots
### if number of LIVE stems in subplot is ever <5 in any sample year, cull all records from this subplot
a <- live[STATUS==1,length(DIAM_T0), by=.(PlotID, SubID, YEAR)]
a[,compID:=paste(PlotID, SubID, sep=".")]
kill.me <- a[V1<5, compID] ## 275 unique subplots (out of total 996) drop below 5 live stems at some point in measurement history
live[,compID:=paste(PlotID, SubID, sep=".")]
live <- live[!(compID%in%kill.me),] ## 13005 records, 7743 second or more measures
# j <- live[,length(unique(YEAR)), by=PlotID] ## how many full revisits happen in each plot
# hist(j$V1) ## all plots are visited at least 2 times, or 3
# View(live[PlotID==6,])
# live[PlotID==6, length(unique(uniqueID)), by=YEAR] ## so we get 24 stems in this plot in 3 measurement events
write.csv(live, "processed/fia.live.stem.dbh.growth.csv")
#####

### Modeling growth~dbh (individual stem level) ## this is purely out of curiosity & need for a comparable figure panel
#####
library(lme4)
live <- as.data.table(read.csv("processed/fia.live.stem.dbh.growth.csv"))  ##  too-sparse subplots removed
live$taxa <- live[,paste(GENUS, SPECIES)]
### new hotness: dbh increment modeling
plot(live$DIAM_T0, live$delta.diam)
plot(live$DIAM_T0, live$diam.rate); mean(live[lag>0, mean(diam.rate)]) ### these boys grow about 0.2 cm/yr
lo <- (lm(diam.rate~DIAM_T0, data=live[lag>0,])) ## very low avg. rate, 0.1 cm/yr, 0.07 gain per 10cm --> contrast andy forest intercept is 0.5 cm/yr
points(live[lag>0,DIAM_T0], predict(lo), col="red", cex=0.3, pch=15)
hm <- (lm(diam.rate~poly(DIAM_T0, degree=2), data=live))
points(live[is.finite(diam.rate),DIAM_T0], predict(hm), col="blue", pch=16, cex=0.2)
summary(hm) ## low predictive, RSE 0.24, R2 0.08, not an improvement from the linear model
live[,diam.rate.rel:=diam.rate/DIAM_T0]
plot(live$DIAM_T0, live$diam.rate) ## same exponential-looking curve -- you're inducing this structure by taking a basically constant value and dividing by x
plot(live[lag>0 & STATUS==1 & diam.rate>0, log(DIAM_T0)], live[lag>0 & STATUS==1 & diam.rate>0, diam.rate])
yyy <- lmer(diam.rate~DIAM_T0 +
       (1|PlotID) +
       (1|GENUS)+
       (1|YEAR), data=live[lag>0 & STATUS ==1,], REML=F)
yyy.null <- lmer(diam.rate~1 +
                   (1|PlotID) +
                   (1|GENUS)+
                   (1|YEAR), data=live[lag>0 & STATUS ==1,], REML=F)

### this model fails to converge ... dodgy estimate quality?
# zzz <- lmer(diam.rate~DIAM_T0 +
#               (1+DIAM_T0|PlotID) +
#               (1+DIAM_T0|GENUS)+
#               (1+DIAM_T0|YEAR), data=live[lag>0 & STATUS ==1,], REML=F)


anova(yyy.null, yyy) ## super significant

save(yyy, file="processed/mod.fia.final.sav")
#####

### a lot of ancillary tests on how stem growth rate varies with time, taxon, etc.
#####
# ## log model of BIOMASS growth rate, growth>0, no hard/soft designation
# summary(live$growth.ann.rel) ## 0.9-3.4%, (-13-53%) -- so looks lower generally than Andy or Street trees
# mod.fia.stem.log <- lm(log(growth.ann.rel)~log(DIAM_T0), data=live[growth.ann.rel>0])
# m1 <- summary(mod.fia.stem.log) #R2 0.03, signficant
# mod.fia.stem.type.log <- lm(log(growth.ann.rel)~log(DIAM_T0)*type, data=live[growth.ann.rel>0])
# m2 <- summary(mod.fia.stem.type.log) # R2 0.04, type H/S not significant
# m3 <- lm(log(growth.ann.rel)~log(DIAM_T0)+type, data=live[growth.ann.rel>0]) ## it IS significant as a categorical factor!
# summary(m3)
# plot(log(live$DIAM_T0), log(live$growth.ann.rel), col=as.numeric(live$type)+1) #S=2, H=1
# abline(a=m2$coefficients[1], b=m2$coefficients[2], lty=2, col="gray55")
# abline(a=m1$coefficients[1], b=m1$coefficients[2], lty=1, col="black")
# ## work in untransformed space (exponential model)
# mod1.nls <- nls(growth.ann.rel ~ exp(a + b * log(DIAM_T0)), data=live, start=list(a=0, b=0)) ### OK THIS is the real exponential non-linear model that can handle the negatives
# mm <- summary(mod1.nls) ## all factors significant
# mod2.nls <- nls(growth.ann.rel~exp(a+b*log(DIAM_T0)+as.numeric(type=="S")*c), data=live, start=list(a=0, b=0, c=0)) ## tiny upward correction for softwoods
# summary(mod2.nls)
# 
# ## an early draft figure for stem diameter growth
# col.type <- c("green", "forestgreen")
# plot(live$DIAM_T0, live$growth.ann.rel, col=col.type[as.numeric(live$type)], pch=15, cex=0.3,
#      xlab="Stem diameter (cm)", ylab="Growth rate (kg/kg)", main="FIA stems")
# test <- seq(live[,min(DIAM_T0)], live[,max(DIAM_T0)], length.out=100)
# lines(test, exp(mm$coefficients[1]+(mm$coefficients[2]*log(test))), col="black", lwd=3)
# legend(x=60, y=0.4, bty="n", legend = c("Hardwood", "Softwood"), fill=c("green", "forestgreen"))
# 
# ## OK there you have it: a common model for growth~dbh, exponential, but low predictive power (lots of spread)
# 
# ### what is mean growth rate at different intervals?
# live[biom.delt.spp>0 & DIAM_T0<20, .(length(DIAM_T0), 
#                                              median(growth.ann.rel, na.rm=T),
#                                              mean(growth.ann.rel, na.rm=T),
#                                              sd(growth.ann.rel, na.rm=T)*1.96)]
# live[biom.delt.spp>0 & DIAM_T0>20 & DIAM_T0<30, .(length(DIAM_T0), 
#                                              median(growth.ann.rel, na.rm=T),
#                                              mean(growth.ann.rel, na.rm=T),
#                                              sd(growth.ann.rel, na.rm=T)*1.96)]
# live[biom.delt.spp>0 & DIAM_T0>30 & DIAM_T0<40, .(length(DIAM_T0), 
#                                              median(growth.ann.rel, na.rm=T),
#                                              mean(growth.ann.rel, na.rm=T),
#                                              sd(growth.ann.rel, na.rm=T)*1.96)]
# live[biom.delt.spp>0 & DIAM_T0>40 & DIAM_T0<50, .(length(DIAM_T0), 
#                                              median(growth.ann.rel, na.rm=T),
#                                              mean(growth.ann.rel, na.rm=T),
#                                              sd(growth.ann.rel, na.rm=T)*1.96)]
# live[biom.delt.spp>0 & DIAM_T0>50 , .(length(DIAM_T0), 
#                                              median(growth.ann.rel, na.rm=T),
#                                              mean(growth.ann.rel, na.rm=T),
#                                              sd(growth.ann.rel, na.rm=T)*1.96)]
# 
# ### so in general there is a lot of variability but a growth rate of about 2ish percent is fine for any diameter -- the model adds litle predictive power to this
# ggplot(live[growth.ann.rel>0,], aes((DIAM_T0), (growth.ann.rel))) + geom_bin2d(bins=80) ## most values are low growth across the diameter range
# ggplot(live[growth.ann.rel>0,], aes(log(DIAM_T0), log(growth.ann.rel))) + geom_bin2d(bins=80)
# 
# ## Any obvious taxonomic split?
# plot(log(live$DIAM_T0), log(live$growth.ann.rel), col=live$GENUS.num) ## no
# plot(log(live$DIAM_T0), log(live$growth.ann.rel), col=live$type) ## no
# table(live[growth.ann.rel>0, GENUS]) ## mostly acer pinus quercus tsuga betula
# live[GENUS%in%c("Quercus", "Pinus", "Betula", "Tsuga", "Acer"), length(DIAM_T0)]/dim(live)[1] #91% in five genera
# ## this is my motivation for applying spp specific allometrics: Large # of Tsuga and Pinus predicted with default hardwood equation
# 
# ## what is stem denstiy in each plot?
# fack <- live[lag>0,.(length(DIAM_T0), sum(biom0.spp, na.rm=T)), by=.(PlotID, YEAR)]
# hist(fack$V1); summary(fack$V1) ## down to the min 5 per subplot, up to 61
# plot(fack$V2, fack$V1) ## generally higher biomass == higher stem count

# ### why not make a nice plot?
# par(mfrow=c(1,1), mar=c(4,4,3,1))
# main.xlim <- c(4, 100)
# main.ylim <- c(-0.15, 1)
# ## FIA stem growth~dbh
# col.type <- c("green", "forestgreen")
# plot(live$DIAM_T0, live$growth.ann.rel, col=col.type[as.numeric(live$type)],
#      pch=15, cex=0.4, xlim=main.xlim, ylim=main.ylim,
#      xlab="Stem DBH (cm)", ylab="Relative growth (kg/kg)", main="FIA stem growth~DBH")
# test <- seq(from=main.xlim[1], to=main.xlim[2], length.out=100)
# lines(test, exp(mm$coefficients[1]+(mm$coefficients[2]*log(test))),
#       col="black", lwd=3, lty=2)
# legend(x=30, y=0.6, bty="n", legend = c("Hardwood", "Softwood"), fill=c("green", "forestgreen"))
# abline(v=live[, median(DIAM_T0)])
# abline(h=live[, median(growth.ann.rel)])
#####

#### Plot-level assessment of growth~biomass.density
#####
### equivalent plot-level assessment in live-only data ## 331 plots have at least one subplot dense enough to bother with
library(data.table)
live <- as.data.table(read.csv("processed/fia.live.stem.dbh.growth.csv")) ## this is data groomed to remove artifacts and subplots with <5 stems
live[DIAM_T0>0 & STATUS==1, .(quantile(DIAM_T0, probs=c(0.05, 0.5, 0.95), na.rm=T), 
                              quantile(diam.rate, probs=c(0.05, 0.5, 0.95), na.rm=T))]
unique(live[DIAM_T0>0 & STATUS==1, unique(YEAR)])
### attempting to track the fate of individual stems (didn't work that well)
#####
# ba.pred <- function(x){(x/2)^2*pi*0.0001} ## get BA per stump from dbh
# ## 296 individual plots
# live[YEAR<2007, unique(lag)] ## nothing has been measured prior to this window
# vets <- live[YEAR<2007, unique(PlotID)] ## 201 plots started before 2007
# arboles[PlotID %in% vets, mean(num.measures)] ## around 2
# arboles[PlotID %in% vets, unique(num.measures)] ## 2-3
# arboles[PlotID %in% vets, ] ## 4521 trees in this set of plots
# # arboles[PlotID %in% vets & num.measures==1,] ## 1315 trees in essentially the whole set of plots, so some trees get missed no matter what
# arboles[PlotID %in% vets & num.measures==3, unique(PlotID)] ### 2286 trees in about half the plots
# arboles[PlotID %in% vets & num.measures>=2, unique(PlotID)] ### 178 of 201 plots have 2-3 measurement events recorded

## 1) figure out t0 for each plot (earliest record, plots started in the first several years)
## 2) Pick out the subplots in t0 that are too sparse and ditch -- in current and subsequent samples
## 3) Calculate total biomass in t0 and number of stems
## 4) Find the next sampling event, calculate total biomass and number of stems
### 4b) Track the frequency of death/removal since last measurement, and toss if too high

# for(p in 1:vets){
#   events <- live[PlotID==vets[p], unique(YEAR)]
#   if(length(events)>1){ ## no sense processing this plot if it has only been visited once
#     a <- live[PlotID==vets[p] & YEAR==min(events), sum(biom1.spp)] ## this is 2006 initial biomass
#     b <- live[PlotID==vets[p] & YEAR==2011, sum(biom1.spp)] ## this is 2006 initial biomass
#     b-a ## 987 kg growth
#     ## compare to inline biomass delta calculation (don't trust how well the query is lining up successive tree measurements)
#     live[PlotID==vets[p] & YEAR==2011, sum(biom.delt.spp)] ## 1025 kg growth
#     live[PlotID==vets[p] & YEAR==min(events), ] ## 14 stems in 2006
#     live[PlotID==vets[p] & YEAR==2011,] ## well well well, 20 stems in 2014
#     live[PlotID==vets[p] & YEAR==min(events), unique(uniqueID)] ## well well well, 20 stems in 2014
#     live[PlotID==vets[p] & YEAR==2011, unique(uniqueID)] ## well well well, 20 stems in 2014
#   }
# }
### there is evidently trees dropping out, and new trees that are growing into the qualifying range over the intervals 
### OK these records may not be solid enough to track tree populations down to individual over time
### INSTEAD: Bulk DBH (biomass) change from t0 (earliest instance) to t1 (next instance)
### remove/account for sparse subplots
### affix the records-based lag for annualizing change
### this approach doesn't pay attention to the fate of individual stems -- these may come and go, but we know how much biomass is there at each sampling event and how many subplots we've included


#####

###
### groom data to get plot-basis biomass change over time
#####
### get the total biomass at t0 and t1 + # of subplots in each sampling event past the first instance
live[lag>0, prev.sample:=YEAR-lag] ## when was plot sampled prior to this sample?
dim(live[lag>0 & STATUS==1,]) ## 7104 LIVE stem measurements (not dead/removed) that have a first-instance already
plot.HW <- live[lag>0 & type=="H" & STATUS==1, .(sum(biom0.spp), sum(biom1.spp)), by=.(PlotID, YEAR)] ## isolate biomass by plot for HW species only
plot.sum <- live[lag>0 & STATUS==1,   ### filter for lag>0 so you only get things on their 2nd or higher sampling event
                 .(sum(biom0.spp),
                   sum(biom1.spp),
                   length(unique(SubID)),
                   mean(prev.sample)), by=.(PlotID, YEAR)]
plot.HW[,yr.Plot:=paste(YEAR, PlotID, sep=".")]
plot.sum[,yr.Plot:=paste(YEAR, PlotID, sep=".")]
plot.sum <- merge(plot.sum, plot.HW[,.(yr.Plot, V1, V2)], by="yr.Plot")
names(plot.sum) <- c("yr.Plot", "PlotID", "YEAR", "biom0.sum", "biom1.sum", "num.subplots", "prev.sample", "biom0.HW", "biom1.HW")
plot.sum[,HWfrac:=biom0.HW/biom0.sum]
subplot.area <- (7.3152^2)*pi
hist(plot.sum$HWfrac) ## most are high HW fraction
hist(plot.sum$num.subplots) ## about evenly split 1-4
dim(plot.sum) ## 319 plots

## biomass gain by AREA
plot.sum[,biom.delt.ann:=((biom1.sum-biom0.sum)/(YEAR-prev.sample))] ## absolute annual bulk biomass change 
plot.sum[,biom.delt.ann.HW:=((biom1.HW-biom0.HW)/(YEAR-prev.sample))] ## absolute annual HW biomass change
plot.sum[,biom.delt.ann.MgC.ha:=((biom.delt.ann/2000)/(subplot.area*num.subplots))*1E4] ## TOTAL gain by area, between 0-5 MgC/ha/yr
plot.sum[,biom.delt.ann.HW.MgC.ha:=((biom.delt.ann.HW/2000)/(subplot.area*num.subplots))*1E4] ## HW gain by area, 0-5 MgC/ha/yr, much more constricted

## biomass gain by prior BIOMASS present
plot.sum[,biom.delt.ann.rel:=biom.delt.ann/biom0.sum] ### annualized total biomass change vs. starting total biomass 
plot.sum[,biom.delt.ann.rel.HW:=biom.delt.ann.HW/biom0.sum] ## annualized HW biomass change vs. starting total biomass
plot.sum[,biom0.MgC.ha:=((biom0.sum/2000)/(subplot.area*num.subplots))*1E4] ## the biomass density avg across all included subplots
plot.sum[,biom0.MgC.ha.HW:=((biom0.HW/2000)/(subplot.area*num.subplots))*1E4]
plot.sum[,biom.delt.ann.rel.HW.HW:=biom.delt.ann.HW/biom0.HW] ## biomass gain relative only to the HW present

## basic exploratory -- what does the data and relationship look like?
length(unique(plot.sum$PlotID)) ## 216 unique plots out of 333 entries, some with more than 1 interval
# mult <- (plot.sum[duplicated(PlotID), PlotID])
# aaaa <- plot.sum[PlotID %in% mult,] ## 107 plots appear twice as successive samples
# View(aaaa[order(PlotID),])
plot.sum[HWfrac>0.25, length(unique(PlotID))] ## 200 plots with 1) at least 1 subplot maintaining 5 or more LIVE stems across measurement period and 2) at least 25% HW biomass at start of measurement window
hist(plot.sum[HWfrac>0.25, biom.delt.ann.rel.HW])
hist(plot.sum[HWfrac>0.25, biom.delt.ann.rel.HW.HW])
write.csv(plot.sum, "processed/fia.live.plot.groomedV2.csv")
#####

## mixed model plot-area growth rate
#####
live.plot <- as.data.table(read.csv("processed/fia.live.plot.groomedV2.csv"))
summary(live.plot[,YEAR-prev.sample])
live.plot[HWfrac>0.25, .(quantile(biom0.MgC.ha, probs=c(0.05, 0.5, 0.95), na.rm=T),
                         quantile(biom.delt.ann.rel.HW, probs=c(0.05, 0.5, 0.95), na.rm=T))]

### modeling growth as MgC/ha/yr(rate)~MgC/ha (density)
plot(live.plot$biom0.MgC.ha, live.plot$biom.delt.ann) ## pretty linear, low slope
points(live.plot[HWfrac>0.25, biom0.MgC.ha], live.plot[HWfrac>0.25, biom.delt.ann.HW], pch=16, col="blue") ## basically flat! ~100 kg/subplot, eliminates high fliers bc pine dominated
points(live.plot[HWfrac<0.25, biom0.MgC.ha], live.plot[HWfrac<0.25, biom.delt.ann], pch=3, col="red") ## some pine dominated plots are v high density and annual absolute gain
live.plot[HWfrac>0.25, mean(biom.delt.ann.HW.MgC.ha)] ## avg. 1.6 MgC/ha/yr -- matches COLE expectations
plot(live.plot[HWfrac>0.25, biom0.MgC.ha], live.plot[HWfrac>0.25, biom.delt.ann.HW.MgC.ha])
check <- lm(biom.delt.ann.HW.MgC.ha~biom0.MgC.ha, data=live.plot[HWfrac>0.25]) ## RSE 0.88
points(live.plot[HWfrac>0.25, biom0.MgC.ha], predict(check), cex=0.4, col="red")
check2 <- lm(biom.delt.ann.HW.MgC.ha~poly(biom0.MgC.ha, degree=2, raw=T), data=live.plot[HWfrac>0.25])
points(live.plot[HWfrac>0.25, biom0.MgC.ha], predict(check2), cex=0.4, col="blue") ## nothing, not significant model
sigma(check) ## 0.88 MgC/ha/yr, R2 0.05

## how different is the scatter based on how we calculate the biomass density or relative biomass gain rate?
plot(live.plot[HWfrac>0.25, biom0.MgC.ha], 
     live.plot[HWfrac>0.25, biom.delt.ann.rel.HW]) ## HW gain rel to TOTAL, vs total biomass
points(live.plot[HWfrac>0.25, biom0.MgC.ha], 
       live.plot[HWfrac>0.25, biom.delt.ann.rel.HW.HW], pch=13, col="blue") ## HW gain rel to HW, vs total biomass
points(live.plot[HWfrac>0.25, biom0.MgC.ha.HW],
       live.plot[HWfrac>0.25, biom.delt.ann.rel.HW.HW], pch=14, col="red") ## HW gain rel to HW, vs HW biomass
summary(live.plot[HWfrac>0.25, biom.delt.ann.rel.HW])
summary(live.plot[HWfrac>0.25, biom.delt.ann.rel.HW.HW]) ## pretty similar

## I think total biomass density is the right metric, and the two metrics of rate are pretty similar

### absolute biomass gain rate
check.me <- lmer(biom.delt.ann.HW.MgC.ha~biom0.MgC.ha+
                   (1|PlotID)+
                   (1|prev.sample), 
                 data=live.plot[HWfrac>0.25,])
summary(check.me) ## nearly identical with the original model, not clear that we can get the repeated part to sort plotID repeated effects
coef(check.me)
sigma(check.me)
fixed.coef <- coef(summary(check.me))
plot(live.plot[HWfrac>0.25, biom0.MgC.ha], live.plot[HWfrac>0.25, biom.delt.ann.HW.MgC.ha])
points(live.plot[HWfrac>0.25, biom0.MgC.ha], predict(check.me), col="blue", pch=12)
b0.dist <- rnorm(1000, mean=fixed.coef[1,1], sd = fixed.coef[1,2])
b1.dist <- rnorm(1000, mean=fixed.coef[2,1], sd = fixed.coef[2,2])
for(i in 1:100){
  b0.r <- rnorm(1, fixed.coef[1,1], fixed.coef[1,2])
  b1.r <- rnorm(1, fixed.coef[2,1], fixed.coef[2,2])
  lines(c(0,250), b0.r+c(0,250)*b1.r, col=i, cex=0.2, pch=16)
} ## OK so this is how it be -- a model with error of linear/low-predictive MgC/ha/yr to MgC/ha relationship


## relative growth rate MgC/MgC (rate)~MgC/ha (density)
library(lme4)
hist(live.plot[,biom0.MgC.ha]) ## looks a lot like urban -- peaks 50-100, long tail up to 300

plot(live.plot$biom0.MgC.ha, live.plot$biom.delt.ann.rel.HW) ## there we are
plot(live.plot$biom0.MgC.ha, live.plot[,log(biom.delt.ann.rel.HW)])
points(plot.sum[HWfrac<0.25, biom0.MgC.ha], plot.sum[HWfrac<0.25, log(biom.delt.ann.rel.HW)], pch=3, col="red") ## low HW = low growth across density -- understory?

### linear model
a <- lmer(biom.delt.ann.rel.HW~biom0.MgC.ha+
       (1|PlotID)+
       (1|prev.sample),
     data=live.plot[HWfrac>0.25], REML=F)
summary(a) ## tiny negative association with density
hist(residuals(a))
## contrast: mean growth in andy edge is 0.064 MgC/MgC, here it is 0.032 MgC/MgC
sigma(a) ## 0.0057
AIC(a)
co <- coef(summary(a))
plot(live.plot[, biom0.MgC.ha], live.plot[, biom.delt.ann.rel.HW])
points(live.plot[, biom0.MgC.ha],
       co[1]+
         (co[2]*live.plot[, biom0.MgC.ha]),
       col="red", pch=16, cex=0.3)

## log xform the response
a.log <- lmer(log(biom.delt.ann.rel.HW)~biom0.MgC.ha+
            (1|PlotID)+
            (1|prev.sample),
          data=live.plot[HWfrac>0.25], REML=F)
summary(a.log) ## tiny negative association with density
summary(lm(log(biom.delt.ann.rel.HW)~biom0.MgC.ha, data=live.plot))
hist(residuals(a.log))
## contrast: mean growth in andy edge is 0.064 MgC/MgC, here it is 0.032 MgC/MgC
sigma(a.log) ## 0.235
AIC(a.log) ##535.8
co.log <- coef(summary(a.log))
plot(live.plot[HWfrac>0.25, biom0.MgC.ha], live.plot[HWfrac>0.25, biom.delt.ann.rel.HW], xlim=c(0,250))
points(live.plot[HWfrac>0.25, biom0.MgC.ha],
       exp(co.log[1]+(co.log[2]*live.plot[HWfrac>0.25, biom0.MgC.ha])),
       col="red", pch=16, cex=0.3)
points(live.plot[HWfrac>0.25, biom0.MgC.ha],
       co[1]+(co[2]*live.plot[HWfrac>0.25, biom0.MgC.ha]),
       col="blue", pch=16, cex=0.3)
a.log.null <- lmer(log(biom.delt.ann.rel.HW)~1+
                (1|PlotID)+
                (1|prev.sample),
              data=live.plot[HWfrac>0.25], REML=F)
anova(a.log.null, a.log, test="Chisq") ## yep significant
coef(summary(a.log))
# ### log xform the predictor
# a.log2 <- lmer(biom.delt.ann.rel.HW~log(biom0.MgC.ha)+
#                (1|PlotID)+
#                (1|prev.sample),
#              data=live.plot[HWfrac>0.25], REML=F)
# AIC(a.log2)
# sigma(a.log2) #0.0059
# summary(a.log2)
# plot(live.plot[, biom0.MgC.ha], live.plot[, biom.delt.ann.rel.HW])
# points(live.plot[, (biom0.MgC.ha)],
#       exp(co[1]+(co[2]*live.plot[, (biom0.MgC.ha)])),
#        col="red", pch=16, cex=0.3)
# points(live.plot[, biom0.MgC.ha],
#        exp(co[1]+(co[2]*live.plot[, biom0.MgC.ha])),
#        col="blue", pch=16, cex=0.3) ## these are the same -- log xform the y or the x


# ## looks bent to me
# b <- lmer(biom.delt.rel.ann.HW~poly(biom0.MgC.ha, degree=2, raw=T)+
#             (1|PlotID)+
#             (1|prev.sample),
#           data=live.plot[,], REML=F)
# summary(b) ## significant coefficients
# sigma(b) ## 0.0055
# AIC(b) ## lower, looking good
# co <- coef(summary(b))
# plot(live.plot[, biom0.MgC.ha], live.plot[, biom.delt.rel.ann.HW]) 
# points(live.plot[, biom0.MgC.ha], 
#        co[1]+
#          (co[2]*live.plot[, biom0.MgC.ha])+
#          (co[3]*(live.plot[, biom0.MgC.ha]^2)),
#        col="red", pch=16, cex=0.3)
# ### filter for low HWfrac
# b.hw <- lmer(biom.delt.rel.ann.HW~poly(biom0.MgC.ha, degree=2, raw=T)+
#             (1|PlotID)+
#             (1|prev.sample),
#           data=live.plot[HWfrac>0.25,], REML=F)
# summary(b.hw) ## significant coefficients
# sigma(b.hw) ## 0.0057
# live.plot[HWfrac>0.25, sd(biom.delt.rel.ann.HW)] ## 0.016 -- so we get some prediction out of the fit
# AIC(b.hw) ## lower, looking good
# co <- coef(summary(b.hw))
# # plot(live.plot[, biom0.MgC.ha], live.plot[, biom.delt.rel.ann.HW]) 
# points(live.plot[, biom0.MgC.ha], ## about the same through 150 MgC/ha, but HW filtered goes higher
#        co[1]+
#          (co[2]*live.plot[, biom0.MgC.ha])+
#          (co[3]*(live.plot[, biom0.MgC.ha]^2)),
#        col="blue", pch=16, cex=0.3)
# ### for this second order poly, let's see what kind of spread we actually get
# plot(live.plot[HWfrac>0.25, biom0.MgC.ha], live.plot[HWfrac>0.25, biom.delt.rel.ann.HW], pch=15, col="black")
# b0.rand <- rnorm(1000, co[1,1], co[1,2])
# b1.rand <- rnorm(1000, co[2,1], co[2,2])
# b2.rand <- rnorm(1000, co[3,1], co[3,2])
# t <- 0:250
# for(i in 1:100){
#   b0.r <- sample(b0.rand, 1)
#   b1.r <- sample(b1.rand, 1)
#   b2.r <- sample(b2.rand, 1)
#   lines(t, b0.r+(b1.r*t)+(b2.r*(t^2)), col=i)
# } ### much less stable than a linear model based on growth per ha
# 
# ## put another kink in it
# c <- lmer(biomHW.gain.rel~poly(biom0.MgC.ha, degree=3, raw=T)+
#             (1|PlotID)+
#             (1|Year),
#           data=dat[PlotID%in%runme,], REML=F)
# summary(c) ## significant coefficients
# sigma(c) ## 0.0062
# AIC(c) ## lower, looking good
# co <- coef(summary(c))
# plot(dat[PlotID%in%runme, biom0.MgC.ha], dat[PlotID%in%runme, biomHW.gain.rel]) 
# points(dat[PlotID%in%runme, biom0.MgC.ha], 
#        co[1]+
#          (co[2]*dat[PlotID%in%runme, biom0.MgC.ha])+
#          (co[3]*(dat[PlotID%in%runme, biom0.MgC.ha]^2))+
#          (co[4]*(dat[PlotID%in%runme, biom0.MgC.ha]^3)),
#        col="red", pch=16, cex=0.3)
# 
# ## do it again
# d <- lmer(biomHW.gain.rel~poly(biom0.MgC.ha, degree=4, raw=T)+
#             (1|PlotID)+
#             (1|Year),
#           data=dat[PlotID%in%runme,], REML=F)
# summary(d) ## significant coefficients
# sigma(d) ## 0.009
# AIC(d) ## lower, looking good
# co <- coef(summary(d))
# plot(dat[PlotID%in%runme, biom0.MgC.ha], dat[PlotID%in%runme, biomHW.gain.rel], ylim=c(-0.01, 0.08)) 
# points(dat[PlotID%in%runme, biom0.MgC.ha], 
#        co[1]+
#          (co[2]*dat[PlotID%in%runme, biom0.MgC.ha])+
#          (co[3]*(dat[PlotID%in%runme, biom0.MgC.ha]^2))+
#          (co[4]*(dat[PlotID%in%runme, biom0.MgC.ha]^3))+
#          (co[5]*(dat[PlotID%in%runme, biom0.MgC.ha]^4)),
#        col="red", pch=16, cex=0.3)
# hist(dat[biomHW.gain.rel>0, biom0.MgC.ha]) ## peak about 100, then long tail up to 300 MgC/ha
# ## I reckon that's alright
# write.csv(dat, "processed/fia.live.plot.groomedV2.csv")
#####




### PART B: BOSTON NPP FROM FIA EMPRICAL GROWTH~BIOMASS ANALYSIS
######
### V2.2: 1) uses species-specific biometrics; 2) models hardwood growth separately from trees in general; 3) uses nls to avoid dumping negative growth records
### V2.3 1) Uses subplot IDs to remove dbh records from subplot sites that are not fully forested; 2) filtered plots that have low hw fraction to determine hw-only growth rate
### V2.4 1) uses all-status tree data, get bulk biomass by plot in each sample event; 2) remove any subplots that do not maintain at >=5 live stems in all sample events; 
### 3) filter plots with <25% HW biomass; 4) simple LMER with random intercepts for starting year and PlotID (some plots are revisited more than once)
### 
##
# live.plot <- as.data.table(read.csv("processed/fia.live.plot.groomed.csv")) ## previous data set not tracking morts/removals or plot/stem multiple records appearances
library(lme4)
library(data.table)
library(raster)
live.plot <- as.data.table(read.csv("processed/fia.live.plot.groomedV2.csv")) ## newest thing
live.plot[HWfrac>0.25,.(min(biom.delt.ann.rel.HW),
                        max(biom.delt.ann.rel.HW))]

## model using biomass gain MgC/MgC logxform
## the below is model a.log above; 297 end-to-end plot records, 200 unique plots
mod.live.plot.final <- lmer(log(biom.delt.ann.rel.HW)~biom0.MgC.ha+
                (1|PlotID)+
                (1|prev.sample),
              data=live.plot[HWfrac>0.25,], REML=F)

plot(live.plot[,biom0.MgC.ha], live.plot[,biom.delt.ann.rel.HW])
points(live.plot[HWfrac<=0.25,biom0.MgC.ha], live.plot[HWfrac<=0.25,biom.delt.ann.rel.HW], pch=12, col="red")


y <- summary(mod.live.plot.final)
save(mod.live.plot.final, file="processed/mod.live.plot.final.sav")

## load the biomass data and reprocess
biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m kg-biomass to 30m pixel
aoi <- raster("processed/boston/bos.aoi30m.tif")
can <- raster("processed/boston/bos.can.redux30m.tif")
isa <- raster("processed/boston/bos.isa30m.tif") ## this is the newly reregistered guy
lulc <- raster("processed/boston/bos.lulc30m.lumped.tif")
biom <- crop(biom, aoi)
can <- crop(can, aoi)
isa <- crop(isa, aoi)
lulc <- crop(lulc,aoi)

biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[, aoi:=as.vector(getValues(aoi))]
biom.dat[, can.frac:=getValues(can)]
biom.dat[, isa.frac:=getValues(isa)]
biom.dat[, lulc:=getValues(lulc)]
biom.dat[, pix.ID:=seq(1, dim(biom.dat)[1])]
names(biom.dat)[1] <- "biom"
dim(biom.dat[!is.na(biom) & biom>10 & aoi>800,]) ## 106659 valid biomass pixels

### Now we have to decide how to calculate the "forest" density of the biomass on the ground
### Three approaches: GROUND (kg-biomass/m2-pixel); "FOREST" (kg-biomass/m2-canopy); "PERV" (kg-biomass/m2-pervious)
## Ground biomass density
biom.dat[, live.MgC.ha.ground:=((biom/2000)/aoi)*(1E4)] ## convert kg-biomass/pixel based on size of pixel (summed by 1m2 in "aoi"), 
biom.dat[biom<10, live.MgC.ha.ground:=0] ## cancel out tiny biomass cells

## Forest biomass density
biom.dat[,live.MgC.ha.forest:=((biom/2000)/(aoi*can.frac))*(1E4)]
biom.dat[biom<10, live.MgC.ha.forest:=0] ## have to manually fix this because of 0 canopy pix
biom.dat[can.frac<0.01, live.MgC.ha.forest:=0]

## Perv biomass density
biom.dat[,live.MgC.ha.perv:=((biom/2000)/(aoi*(1-isa.frac)))*(1E4)]
biom.dat[biom<10, live.MgC.ha.perv:=0] ## have to manually fix this because of isa=1 pix
biom.dat[isa.frac>0.99, live.MgC.ha.perv:=0]

# ### what do our densities look like
# hist(biom.dat[aoi>800 & live.MgC.ha.ground>0, live.MgC.ha.ground]) ## up to 600 Mgbiom/ha = 300 MgC/ha -- a lot of this below the range sampled by the FIA plot data
# summary(biom.dat[aoi>800, live.MgC.ha.ground]) ## mostly below 40 MgC/ha
# hist(biom.dat[aoi>800 & live.MgC.ha.forest>0, live.MgC.ha.forest]) ## 0-300 MgC/ha, nice normal peak ~80
# summary(biom.dat[aoi>800,live.MgC.ha.forest]) ## everything nudged denser
# hist(biom.dat[aoi>800 & isa.frac<0.98 & live.MgC.ha.perv>0, live.MgC.ha.perv]) ## up to 3000 MgC/ha skewed as fuck
# summary(biom.dat[aoi>800, live.MgC.ha.perv]) ## a few very high
# biom.dat[aoi>800 & live.MgC.ha.perv>300,] ### 3k v high, these are the classic artifacts, v. low pervious cover, but low overall biomass
# ### contrast: what is the range sampled in the live.plot FIA data?
# summary(live.plot[HWfrac>0.25,biom0.MgC.ha]) ## don't get into the low range of things! Minimum biomass density is about where most of the pixels fall
# hist(live.plot[HWfrac>0.25,biom0.MgC.ha]) ## looks a lot like the forest spread in my data -- but these are presumably wall to wall trees!

# ## root out any artifacts in density calcs
# biom.dat[aoi>800, length(aoi)] #136667 valid pixels
# biom.dat[aoi>800, length(bos.biom30m)] #136667 valid pixels
# biom.dat[aoi>800 & is.finite(bos.biom30m), length(bos.biom30m)] ##135705 is the number to hit
# biom.dat[!is.na(live.MgC.ha.ground) & aoi>800,] ## 135705 pix
# biom.dat[is.finite(live.MgC.ha.ground) & aoi>800,] ## 135705 pix
# biom.dat[is.finite(live.MgC.ha.perv) & aoi>800,] ## 135714 pix
# ### this is biomass associated with 100% paved pixels, added fix above

## calculate growth factors per cell
### new hotness: build in model estimate uncertainty

### loop and interatively make maps with different model realizations
max.fact <- live.plot[HWfrac>0.25, max(biom.delt.ann.rel.HW)]
min.fact <- live.plot[HWfrac>0.25, min(biom.delt.ann.rel.HW)]
b0.rand <- rnorm(1000, coef(y)[1,1], coef(y)[1,2])
b1.rand <- rnorm(1000, coef(y)[2,1], coef(y)[2,2])

# ## what is the expected general range of gfacts given the distribution of biomass densities?
# gfacts.ground.mn <- exp(mean(b0.rand)+biom.dat[aoi>800 & live.MgC.ha.ground>0, live.MgC.ha.ground]*mean(b1.rand))
# gfacts.ground.mn[gfacts.ground.mn>max.fact] <- max.fact
# hist(gfacts.ground.mn); summary(gfacts.ground.mn) ### skew left, a lot pile up at the high end but none hit the limit
# plot(biom.dat[aoi>800 & live.MgC.ha.ground>0, live.MgC.ha.ground], gfacts.ground.mn) ## looks just like the fixed fit
# plot(live.plot[HWfrac>0.25, biom0.MgC.ha], live.plot[HWfrac>0.25, biom.delt.ann.rel.HW], xlim=c(0,250))
# points(live.plot[HWfrac>0.25, biom0.MgC.ha], exp(mean(b0.rand)+(live.plot[HWfrac>0.25, biom0.MgC.ha]*mean(b1.rand))), col="red", cex=0.4, pch=15)
# 
# gfacts.forest.mn <- exp(mean(b0.rand)+biom.dat[aoi>800 & live.MgC.ha.forest>0, live.MgC.ha.forest]*mean(b1.rand))
# gfacts.forest.mn[gfacts.forest.mn>max.fact] <- max.fact
# hist(gfacts.forest.mn); summary(gfacts.forest.mn) ### a pretty pretty bell
# plot(biom.dat[aoi>800 & live.MgC.ha.forest>0, live.MgC.ha.forest], gfacts.forest.mn) ## looks just like the fixed fit
# plot(live.plot[HWfrac>0.25, biom0.MgC.ha], live.plot[HWfrac>0.25, biom.delt.ann.rel.HW], xlim=c(0,250))
# points(live.plot[HWfrac>0.25, biom0.MgC.ha], exp(mean(b0.rand)+(live.plot[HWfrac>0.25, biom0.MgC.ha]*mean(b1.rand))), col="red", cex=0.4, pch=15)
# 
# ## a comparison of histograms, or: Why it do be like it is sometimes.
# par(mar=c(2,2,2,2), oma=c(0,0,0,0), mfrow=c(2,2))
# hist(biom.dat[aoi>800 & live.MgC.ha.ground>0, live.MgC.ha.ground])
# hist(gfacts.ground.mn)
# hist(biom.dat[aoi>800 & live.MgC.ha.forest>0, live.MgC.ha.forest])
# hist(gfacts.forest.mn)

for(i in 1:1000){
  b0.sel <- b0.rand[i]
  b1.sel <- b1.rand[i]
  dump <- copy(biom.dat)
  ## assign growth factors based on position and biomass density
  dump[,ground.gfact:=exp(b0.sel+(b1.sel*live.MgC.ha.ground))] ## predict growth factors ground
  dump[(live.MgC.ha.ground)<0.5, ground.gfact:=0] ## anything without less than about 9 kg biomass/cell gets 0 factor
  dump[,forest.gfact:=exp(b0.sel+(b1.sel*live.MgC.ha.forest))] ## predict growth factors ground
  dump[(live.MgC.ha.forest)<0.5, forest.gfact:=0]
  dump[can.frac<0.01, forest.gfact:=0]
  dump[,perv.gfact:=exp(b0.sel+(b1.sel*live.MgC.ha.perv))] ## predict growth factors ground
  dump[(live.MgC.ha.perv)<0.5, perv.gfact:=0]
  dump[isa.frac>0.99, perv.gfact:=0]
  
  ### all of these are limited to the minimum, but can reach too high maxima 
  dump[ground.gfact>max.fact, ground.gfact:=max.fact]
  dump[forest.gfact>max.fact, forest.gfact:=max.fact]
  dump[perv.gfact>max.fact, perv.gfact:=max.fact]
  
    ## calculate NPP and store
  ### biomass productivity factors
  dump[,npp.kg.hw.ground:=biom*ground.gfact] ## apply ground growth factor to whole pixel area
  dump[,npp.kg.hw.forest:=biom*forest.gfact] ## apply growth factor to FOREST AREA
  dump[,npp.kg.hw.perv:=biom*perv.gfact] ## apply growth factor to PERVIOUS AREA
  names(dump)[13:15] <- paste0("fia.npp.", c("ground", "forest", "perv"), ".iter", i, ".csv")

  fwrite(dump[,13], paste0("processed/results/fia/ground/fia.npp.ground.iter", i, ".csv"), na = NA)
  fwrite(dump[,14], paste0("processed/results/fia/forest/fia.npp.forest.iter", i, ".csv"), na = NA)
  fwrite(dump[,15], paste0("processed/results/fia/perv/fia.npp.perv.iter", i, ".csv"), na = NA)
  print(paste("FIA NPP iteration",i))
}
dim(dump[biom>10 & !is.na(biom) & aoi>800,]) ## 106659 valid pix still

## rebuild the frames
### Version history: 
### V24 = static model 
### V3 = 2nd poly to predict relative growth rate, with model noise;
### V4 = linear model to predict area growth rate (MgC/ha/yr) with model noise
### V5 = log-xform linear model to predict relative growth rate, with model noise
### V6 = biomass>0 canopy map, 1000 iterations, log-xform linear model to predict relative growth rate, with model noise

tmp <- matrix(nrow=dim(biom.dat)[1], ncol=1009)
tmp[,1:9] <- unlist(biom.dat)
for(d in 1:1000){
  a <- fread(paste0("processed/results/fia/ground/fia.npp.ground.iter", d, ".csv"))
  tmp[,d+9] <- unlist(a)
  print(paste("attached iteration", d, "ground"))
}
tmp <- as.data.frame(tmp)
colnames(tmp)[1:9] <- names(biom.dat)
colnames(tmp)[10:1009] <- paste0("fia.npp.ground.iter", 1:1000)
fwrite(tmp, "processed/results/fia/npp.FIA.empirV6.ground.csv", na=NA)

tmp <- matrix(nrow=dim(biom.dat)[1], ncol=1009)
tmp[,1:9] <- unlist(biom.dat)
for(d in 1:1000){
  a <- fread(paste0("processed/results/fia/forest/fia.npp.forest.iter", d, ".csv"))
  tmp[,d+9] <- unlist(a)
  print(paste("attached iteration", d, "forest"))
}
tmp <- as.data.frame(tmp)
colnames(tmp)[1:9] <- names(biom.dat)
colnames(tmp)[10:1009] <- paste0("fia.npp.forest.iter", 1:1000)
fwrite(tmp, "processed/results/fia/npp.FIA.empirV6.forest.csv", na=NA)

tmp <- matrix(nrow=dim(biom.dat)[1], ncol=1009)
tmp[,1:9] <- unlist(biom.dat)
for(d in 1:1000){
  a <- fread(paste0("processed/results/fia/perv/fia.npp.perv.iter", d, ".csv"))
  tmp[,d+9] <- unlist(a)
  print(paste("attached iteration", d, "perv"))
}
tmp <- as.data.frame(tmp)
colnames(tmp)[1:9] <- names(biom.dat)
colnames(tmp)[10:1009] <- paste0("fia.npp.perv.iter", 1:1000)
fwrite(tmp, "processed/results/fia/npp.FIA.empirV6.perv.csv", na=NA)


### what is the general distribution of predicted growth factors for the forest FIA runs?
# tmp.fia <- fread("processed/results/fia/npp.FIA.empirV6.forest.csv")
# dim(tmp.fia); names(tmp.fia) ## 354068 total pixels
# dim(tmp.fia[!is.na(biom) & aoi>800,]) ## 135705 biomass pixels
# dim(tmp.fia[aoi>800 & !is.na(fia.npp.forest.iter1) & !is.na(biom),]) ## 135705 valid retrievals
# plot(tmp.fia[,biom], tmp.fia[,fia.npp.forest.iter1]) ## a bell curve, nice
# median.na <- function(x){median(x, na.rm=T)}
# npp.med <- apply(as.matrix(tmp.fia[,10:1009]), MARGIN=1, FUN=median.na)
# hist(npp.med) ## seems believable
# length(npp.med)## 354068, so it's the full raster
# aa <- npp.med-tmp.fia$fia.npp.forest.iter51 ## this appears to truly be the median
# summary(aa) 
# sum(npp.med, na.rm=T)/2000 ## 4.8 ktC

### write out the median rasters
tmp <- fread("processed/results/fia/npp.FIA.empirV6.ground.csv")
median.na <- function(x){median(x, na.rm=T)}
npp.med <- apply(as.matrix(tmp[,10:1009]), MARGIN=1, FUN=median.na)
r <- raster(aoi)
r <- setValues(r, values = npp.med)
r <- mask(r, mask=aoi)
writeRaster(r, filename="processed/results/fia/npp.FIA.empirV6.ground.median.tif", format="GTiff", overwrite=T)
sum(npp.med, na.rm=T)/2000 ## 7.0 ktC

tmp <- fread("processed/results/fia/npp.FIA.empirV6.forest.csv")
median.na <- function(x){median(x, na.rm=T)}
npp.med <- apply(as.matrix(tmp[,10:1009]), MARGIN=1, FUN=median.na)
r <- raster(aoi)
r <- setValues(r, values = npp.med)
r <- mask(r, mask=aoi)
writeRaster(r, filename="processed/results/fia/npp.FIA.empirV6.forest.median.tif", format="GTiff", overwrite=T)
sum(npp.med, na.rm=T)/2000 ## sum 4.8 ktC

tmp <- fread("processed/results/fia/npp.FIA.empirV6.perv.csv")
median.na <- function(x){median(x, na.rm=T)}
npp.med <- apply(as.matrix(tmp[,10:1009]), MARGIN=1, FUN=median.na)
r <- raster(aoi)
r <- setValues(r, values = npp.med)
r <- mask(r, mask=aoi)
writeRaster(r, filename="processed/results/fia/npp.FIA.empirV6.perv.median.tif", format="GTiff", overwrite=T)
sum(npp.med, na.rm=T)/2000 ## sum 5.13 ktC



#### Supplemental analysis
### following section contains:
### 1) Analysis of weird artifact seen in stem-level growth~dbh plots (log-log) -- why parallel lines?
### 2) Initial Boston NPP estimation based on published biomass-->Age-->C.aquisition equations provided by FIA COLE search
#####

### 1) Why do stem growth~dbh plots have distinct parallel lines?
# #### what is the deal with these scatter plots -- parallel lines, artifact bug hunt
# dat <- read.csv("data/FIA/MA_Tree_Data_ID_NOMORT_SUBID.csv")
# dat <- as.data.table(dat)
# plot(dat$DIAM_T0, log(dat$biom.delt))
# plot(log(dat$DIAM_T0), log(dat$biom.delt)) ## linear but only in the ones that are not falling (non negative)
# plot(log(dat$DIAM_T0), log(dat$biom.delt),
#      col=dat$GENUS.num) ## a lot of one particular genus
# plot(log(dat$DIAM_T0), log(dat$biom.delt),
#      col=dat$SPECIES_CD) ## no clear plot ID in lower parts
# plot(dat$DIAM_T0, dat$biom.delt, col=dat$GENUS) ## some are diving
# dat[biom.delt<(-500),] ## only 2 lose more than 500kg
# View(dat[biom.delt<exp(1) & biom.delt>0,]) ## tiny growth are larger pinus growing by 0.1 inch
# plot(dat[log(biom.delt)<1.45 & log(DIAM_T0)<3, log(DIAM_T0)], dat[log(biom.delt)<1.45 & log(DIAM_T0)<3, log(biom.delt)])
# dat[,diam.delt:=DIAM_T1-DIAM_T0]

# ### here's what's happening: The low straight-line changes are where they are detecting a 0.1" increase in circum. (minimum increase on dbh tape) -- locks things into a single growth line
# table(dat[log(biom.delt)<1.5 & log(DIAM_T0)<3, GENUS]) ## dominated by a handfull of genera
# table(dat[log(biom.delt)<1.5 & log(DIAM_T0)<3, PlotID]) ## these low growth small ones are everywhere
# hist(dat[log(biom.delt)<1.5 & log(DIAM_T0)<3, DIAM_T0]) ## fairly even dbh dist between 12 and 20cm
# plts <- unique(dat$PlotID)
# for(i in 1:length(plts)){ ### in general, a positive trend of greater delta.biom with greater dbh to start
#   plot(dat[PlotID==plts[i], log(DIAM_T0)], 
#        dat[PlotID==plts[i],log(biom.delt)], 
#        col=dat[PlotID==plts[i],GENUS.num], main=paste("plot", plts[i]))
# }
# ### hang on: a lot of these are on the same slope
# hist(dat[dbh.delt>0, dbh.delt], breaks=50)
# dat[dbh.delt>2,]
# dat[,dbh.delt.incr:=dbh.delt/0.254] ## 0.245 is the dbh tape increment
# dat[dbh.delt.incr>0 & dbh.delt.incr<3,]
# plot(dat[,log(DIAM_T0)], dat[,log(biom.rel)], col=as.factor(dat[,dbh.delt.incr]))
# View(dat[log(DIAM_T0)<2.8&log(biom.rel)<(-4.1),]) ## a huge proportion of these are 0 growth
# dat[biom.rel<=0,] ### 2k records show 0 to negative growth
# plot(dat[log(DIAM_T0)<2.8&log(biom.rel)<(-4.1), log(DIAM_T0)], 
#      dat[log(DIAM_T0)<2.8&log(biom.rel)<(-4.1), log(biom.rel)],
#      col=as.factor(dat[,dbh.delt.incr]))
# r <- dat[log(DIAM_T0)<2.55&log(biom.rel)<(-4.1),]
# View(r[order(r$dbh.delt.incr),])
# dat[, dbh.units0:=DIAM_T0/0.254] ## convert dbh to tape increment units
# dat[, dbh.units1:=DIAM_T1/0.254]
# plot(dat$dbh.units0, dat$dbh.units1) ### these are def integers
# plot(dat[dbh.units1-dbh.units0>0, dbh.units0], dat[dbh.units1-dbh.units0>0, dbh.units1])
# abline(a=0,b=1) ## positive growth
# plot(dat[dbh.units1-dbh.units0<0, dbh.units0], dat[dbh.units1-dbh.units0<0, dbh.units1])
# abline(a=0,b=1) ## negative growth
# plot(dat[dbh.units1-dbh.units0==0, dbh.units0], dat[dbh.units1-dbh.units0==0, dbh.units1])
# abline(a=0,b=1) ## zero growth
# 
# dat[,dbh.unit.delt:=dbh.units1-dbh.units0]
# dat[log(DIAM_T0)<2.55&log(biom.rel)<(-4.1), unique(dbh.unit.delt)] ## for the low lines it's all 0 or 1 dbh increment
# t.col <- c("blue", "purple", "red")
# plot(dat[log(DIAM_T0)<2.8&log(biom.rel)<(-4.1), log(DIAM_T0)], ### everything is 1 dbh increment growth
#      dat[log(DIAM_T0)<2.8&log(biom.rel)<(-4.1), log(biom.rel)],
#      col=t.col[dat[log(DIAM_T0)<2.8&log(biom.rel)<(-4.1), dbh.unit.delt]])
# 
# View(dat[log(DIAM_T0)<2.8&log(biom.rel)<(-4.1) & log(biom.rel)>(-4.2), ])
# View(dat[log(DIAM_T0)<2.8&log(biom.rel)<(-4.2), ])
# plot(dat[DIAM_T0<exp(2.8) & biom.rel<exp(-4.2), DIAM_T0],
#      dat[DIAM_T0<exp(2.8) & biom.rel<exp(-4.2), biom.rel],
#      col=t.col[abs(dat[DIAM_T0<exp(2.8) & biom.rel<exp(-4.2), dbh.unit.delt])])
# 
# t.col <- c("blue", "red", "black", "green", "orange", "yellow", rep("salmon",100))
# plot(dat[, DIAM_T0],
#      dat[, biom.rel],
#      col=t.col[abs(dat[, dbh.unit.delt])])
# plot(dat[biom.rel>0, DIAM_T0],
#      dat[biom.rel>0, biom.rel],
#      col=t.col[abs(dat[biom.rel>0, dbh.unit.delt])])
# plot(dat[dbh.unit.delt==1, DIAM_T0],
#      dat[dbh.unit.delt==1, biom.rel], ylim=c(-0.05, 0.05), col="green")
# points(dat[dbh.unit.delt==0, DIAM_T0],
#      dat[dbh.unit.delt==0, biom.rel], col="black") ## look like they're on top of each other
# points(dat[dbh.unit.delt==2, DIAM_T0],
#        dat[dbh.unit.delt==2, biom.rel], col="blue") ## look like they're on top of each other
# points(dat[dbh.unit.delt==3, DIAM_T0],
#        dat[dbh.unit.delt==3, biom.rel], col="purple") ## look like they're on top of each other
# points(dat[dbh.unit.delt==-1, DIAM_T0],
#        dat[dbh.unit.delt==-1, biom.rel], col="orange") ## look like they're on top of each other
# points(dat[dbh.unit.delt==-2, DIAM_T0],
#        dat[dbh.unit.delt==-2, biom.rel], col="gold") ## look like they're on top of each other
# points(dat[dbh.unit.delt==-3, DIAM_T0],
#        dat[dbh.unit.delt==-3, biom.rel], col="pink") ## look like they're on top of each other
## congratulations, you've proved that the dbh~biomass allometric is exponential
## basically every linear line you get is a new click on the line: Higher => jumped more dbh units
## but in the log transform you are cutting off all the values that are negative dbh units

# plot(dat[dbh.unit.delt%in%seq(-4,4), DIAM_T0],
#      dat[dbh.unit.delt%in%seq(-4,4), biom.rel], col=dat[dbh.unit.delt%in%seq(-4,4), dbh.unit.delt+5])
# 
# plot(dat[dbh.unit.delt%in%seq(-10,10), DIAM_T0],
#      dat[dbh.unit.delt%in%seq(-10,10), biom.rel], col=dat[dbh.unit.delt%in%seq(-10,10), dbh.unit.delt+11])
# 
# dat[dbh.unit.delt%in%seq(-10,10), ] ## 8018 records of 9402 (85%)
# ### log transforming it gives you this one (byproduct is it elminates the negative growth records entirely)
# plot(dat[dbh.unit.delt%in%seq(1,10), log(DIAM_T0)],
#      dat[dbh.unit.delt%in%seq(1,10), log(biom.rel)], col=dat[dbh.unit.delt%in%seq(1,10), dbh.unit.delt])
# plot(log(dat$DIAM_T0), log(dat$biom.rel)) ## yep, each separate line is a different dbh increment class
# plot(dat[dbh.unit.delt%in%seq(1,50), log(biom0)], ## same feckin story
#      dat[dbh.unit.delt%in%seq(1,50), log(biom.rel)], col=dat[dbh.unit.delt%in%seq(1,50), dbh.unit.delt])
# ## so counting up from the bottom, each line represents 1, 2, 3, etc dbh increment gain, across a range of starting dbh (and bigger dbh class gives *relatively* less biomass gain at a given dbh increment class, but absolutely more)


### SUPPLEMENTAL 2: FIA, EQUATION BASIS BIOMASS-->AGE-->GROWTH
#####
###
## FIA V1: Equation for wood volume~time (proxy for stand biomass density)
a=123.67
b=0.04
gogo <- 120
AGE=seq(0,gogo)
y = a*(1-exp(-b*AGE))^3
plot(AGE, y, ylab="stand MgC/ha")
z = diff(y) ## yearly growth increment, absolute kg
z.rel <- z/y[2:(gogo+1)] ## yearly growth increment, % of previous year biomass
plot(y[2:(gogo+1)], z,  xlab="biomass MgC/ha", ylab="absolute growth rate, MgC/ha/yr") ## zero growth above ~125 MgC/ha
points(y[2:(gogo+1)], z.rel, pch=14, cex=0.5, col="red") #hyperbolic
plot(AGE[2:(gogo+1)], z, xlab="Stand age, yr", ylab="absolute growth rate, MgC/ha/yr") ## this is the gain curve over site age, near zero >200 yrs
### above 250 yrs gain is <1 kgC/ha/yr
### what is equivalent age of stand if we know live biomass?
log(1-(y/a)^(1/3))/(-b) ## predicts 40 years
tC <- seq(0,120)
plot(tC, log(1-(tC/a)^(1/3))/(-b), xlab="live biomass, tC/ha", ylab="Site age, yr")

###
### standing live wood C storage in local forest that resemble (?) what we'd have in Boston:
### i.e. N.red oak; red maple/oak; mixed upland hwoods; red maple uplands: Range is  94.7-105.1 MgC/ha
## read in summaries of C stock change in Q.alba/Q.rubra/Carya and Q.rubra forest
## both forest types include values for reforestation and afforestation conditions
# ne.summ <- read.csv("docs/FIA_CTMANHRI.csv")
# plot(ne.summ$Age.yrs, ne.summ$live.tree.tCha)
# 
# mix.oak.ref <- read.csv("docs/FIA_QUALQURUCATO_Reforest.csv")
# plot(mix.oak.ref$Age.yrs, mix.oak.ref$live.tree.c.inc)
# mix.oak.aff <- read.csv("docs/FIA_QUALQURUCATO_Afforest.csv")
# plot(mix.oak.aff$Age.yrs, mix.oak.aff$live.tree.c.inc)
# 
# red.oak.ref <- read.csv("docs/FIA_QURU_Reforest.csv")
# plot(red.oak.ref$Age.yrs, red.oak.ref$live.tree.c.inc)
# red.oak.aff <- read.csv("docs/FIA_QURU_Afforest.csv")
# plot(red.oak.aff$Age.yrs, red.oak.aff$live.tree.c.inc)
# 
# fia.summ <- as.data.frame(cbind(ne.summ$Age.yrs, ne.summ$live.tree.c.inc, mix.oak.ref$live.tree.c.inc,
#                                 mix.oak.aff$live.tree.c.inc, red.oak.ref$live.tree.c.inc, red.oak.aff$live.tree.c.inc))
# colnames(fia.summ) <- c("Age", "NE.total", "mix.oak.ref", "mix.oak.aff", "red.oak.ref", "red.oak.aff")
# plot(fia.summ$Age, fia.summ$NE.total, pch=15, col="black")
# points(fia.summ$Age, fia.summ$mix.oak.aff, pch=16, col="lightblue")
# points(fia.summ$Age, fia.summ$mix.oak.ref, pch=16, col="blue")
# points(fia.summ$Age, fia.summ$red.oak.aff, pch=17, col="pink")
# points(fia.summ$Age, fia.summ$red.oak.ref, pch=17, col="red")
### the only thing that changes between reforestation and afforestation are values for forest floor and soil C

# ## what is relationship between standing live biomass-C and C increment
# plot(ne.summ$mean.vol.m3, ne.summ$live.tree.c.inc)
# plot(ne.summ$live.tree.tCha, ne.summ$live.tree.c.inc)
# plot(ne.summ$Age.yrs, ne.summ$live.tree.tCha)
# plot(ne.summ$Age.yrs, ne.summ$mean.vol.m3) ## basic sigmoid 0 to max at 100

###
### prototype process for using FIA aggregate data to figure npp from Raciti biomass
# 1) determine MgC/ha in 30m plot, normalize by appropriate area factor (raw ground/canopy/pervious) (i.e. there is X many ha of forest there with Y much living carbon in place)
# 2) Figure out the next-year tC/ha in the plot
# 3) apply this incrememnt to the area fraction in question

biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m kg-biomass to 30m pixel
aoi <- raster("processed/boston/bos.aoi30m.tif")
can <- raster("processed/boston/bos.can30m.tif")
isa <- raster("processed/boston/bos.isa.rereg30m.tif")
biom <- crop(biom, aoi)
can <- crop(can, aoi)
isa <- extend(isa, aoi)
biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[,aoi:=as.vector(getValues(aoi))]
can.dat <- as.data.table(as.data.frame(can))
biom.dat[, can.frac:=can.dat$bos.can30m]
isa.dat <- as.data.table(as.data.frame(isa))
biom.dat[, isa.frac:=isa.dat$bos.isa.rereg30m]
biom.dat[,pix.ID:=seq(1, dim(biom.dat)[1])]

### live MgC by area of GROUND in each pixel
biom.dat[, live.MgC.ha.ground:=(bos.biom30m/aoi)*(1/2)*(1E-03)*(1E4)] ## convert kg-biomass/pixel based on size of pixel (summed by 1m2 in "aoi"), kgC:kgbiomass, Mg:kg, m2:ha
biom.dat[aoi>800, range(live.MgC.ha.ground, na.rm=T)] ## up to 284 MgC/ha, as we'd expect from what we saw in Raciti
hist(biom.dat[aoi>800,live.MgC.ha.ground]) #v skewed, most below 50 MgC/ha
biom.MgC.ha.ground <- raster(biom)
biom.MgC.ha.ground <- setValues(biom.MgC.ha.ground, biom.dat[,live.MgC.ha.ground])

### ALTERNATIVE BIOMASS DENSITIES: correct biomass figures for the amount of canopy or pervious cover present per cell
### i.e. we assume FIA is measuring trees in "forest" with essentially continuous canopy, so that differences in tC/ha as a function of age are purely a function of tree growth and not differences in tree %coverage

### live MgC/ha for "forest" fraction in each pixel
biom.dat[,live.MgC.ha.forest:=(bos.biom30m/(aoi*can.frac))*(1/2)*(1E-03)*(1E4)]
biom.dat[bos.biom30m==0, live.MgC.ha.forest:=0] ## have to manually fix this because of 0 canopy pix
# biom.dat[can.frac<=0.01, live.MgC.ha.forest:=0]
range(biom.dat[aoi>800,live.MgC.ha.forest],  na.rm=T) ## 0 - 284
hist(biom.dat[aoi>800,live.MgC.ha.forest]) ## correcting for canopy cover, more mid-rage values
biom.dat[live.MgC.ha.forest<100, length(live.MgC.ha.forest)]/dim(biom.dat[!is.na(can.frac),])[1] ## 84% of pixels are below 100 MgC/ha
biom.MgC.ha.forest <- raster(biom)
biom.MgC.ha.forest <- setValues(biom.MgC.ha.forest, biom.dat[,live.MgC.ha.forest])
plot(biom.MgC.ha.forest)

## correct for pervious cover
biom.dat[,live.MgC.ha.perv:=(bos.biom30m/(aoi*(1-isa.frac)))*(1/2)*(1E-03)*(1E4)]
biom.dat[bos.biom30m==0, live.MgC.ha.perv:=0] ## have to manually fix this because of isa=1 pix

# biom.dat[isa.frac>0.99, live.MgC.ha.perv:=0]
range(biom.dat[aoi>800 & isa.frac<0.98,live.MgC.ha.perv],  na.rm=T) ## 0 - 3890
hist(biom.dat[aoi>800 & isa.frac<0.98,live.MgC.ha.perv]) ## a small number of very extreme values
biom.dat[live.MgC.ha.perv<100, length(live.MgC.ha.perv)]/dim(biom.dat[!is.na(can.frac),])[1] ## 75% of pixels are below 100 MgC/ha
biom.MgC.ha.perv <- raster(biom)
biom.MgC.ha.perv <- setValues(biom.MgC.ha.perv, biom.dat[,live.MgC.ha.perv])
plot(biom.MgC.ha.perv)

### get delta figures
biom.dat[, delta.C.perv:=live.MgC.ha.perv-(live.MgC.ha.ground)]
biom.dat[, delta.C.forest:=live.MgC.ha.forest-(live.MgC.ha.ground)]
plot(biom.dat[isa.frac<0.9, isa.frac], biom.dat[isa.frac<0.9, delta.C.perv]) ## deviation using pervious correction gets higher with greater impervious fraction
plot(biom.dat[can.frac>0.07, can.frac], biom.dat[can.frac>0.07, delta.C.forest]) ## as canopy nears 100%, NPP estimates converge on raw area

### figure out forest "age" for the cells (using coefficients for NE total)
## age based on ground area
a=123.67
b=0.04
biom.dat[,age.ground:=log(1-(live.MgC.ha.ground/a)^(1/3))/(-b)] ## some of these will produce infinities with too dense biomass
# biom.dat[age.ground>120, age.ground:=120] ## fix the divergent ones to just "old, not growing"
# biom.dat[!is.finite(age.ground), age.ground:=120] ## again, fix the ones that got fucked to "old, not growing"
# biom.dat[age.ground>250, age.ground:=NA] ## don't cancel the high ages -- need to see them in order to fix them in post-process
biom.dat[is.na(aoi), age.ground:=NA] # cancel places out of bounds
# biom.dat[bos.biom30m<=10, age.ground:=NA] ## cancel places with no biomass
biom.dat[is.na(bos.biom30m), age.ground:=NA]

ground.age <- raster(biom)
ground.age <- setValues(ground.age, biom.dat[,age.ground])
plot(ground.age)
hist(biom.dat[,age.ground]) ## got a lot of "old" ones, indicating high density of biomass

## age based on canopy area
biom.dat[,age.forest:=log(1-(live.MgC.ha.forest/a)^(1/3))/(-b)] ## some of these will produce infinities with too dense biomass
# biom.dat[age.forest>120, age.forest:=120] ## fix the divergent ones to just "old, not growing"
# biom.dat[!is.finite(age.forest), age.forest:=120] ## again, fix the ones that got fucked to "old, not growing"
# biom.dat[age.forest>250, age.forest:=NA] ## cancel ages that are unreliably retrieved
biom.dat[is.na(aoi), age.forest:=NA] # cancel places out of bounds
# biom.dat[bos.biom30m<10, age.forest:=NA] ## cancel places with no biomass
biom.dat[is.na(bos.biom30m), age.forest:=NA]

forest.age <- raster(biom)
forest.age <- setValues(forest.age, biom.dat[,age.forest])
plot(forest.age)
hist(biom.dat[,age.forest]) ## many more old forest, peak has moved older

## age based on pervious area
biom.dat[,age.perv:=log(1-(live.MgC.ha.perv/a)^(1/3))/(-b)] ## some of these will produce infinities with too dense biomass
# biom.dat[age.perv>120, age.perv:=120] ## fix the divergent ones to just "old, not growing"
# biom.dat[!is.finite(age.perv), age.perv:=120] ## again, fix the ones that got fucked to "old, not growing"
# biom.dat[age.perv>250, age.perv:=NA] ## cancel ages that are unreliably retrieved
biom.dat[is.na(aoi), age.perv:=NA] # cancel places out of bounds
# biom.dat[bos.biom30m<=10, age.perv:=NA] ## cancel places with no biomass
biom.dat[is.na(bos.biom30m), age.perv:=NA]

perv.age <- raster(biom)
perv.age <- setValues(perv.age, biom.dat[,age.perv])
plot(perv.age)
hist(biom.dat[,age.perv])

### compare
biom.dat[,delta.age.perv:=age.perv-age.ground]
biom.dat[,delta.age.forest:=age.forest-age.ground]
plot(biom.dat$bos.biom30m, biom.dat$delta.age.forest)
plot(biom.dat$bos.biom30m, biom.dat$delta.age.perv)

## frequency distributions of different methods
par(mfrow=c(3,1), mar=c(4,4,2,1))
hist(biom.dat$age.ground,  main="Forest age, unadjusted", xlab="Age, yrs", xlim=c(0, 300), breaks=40)
hist(biom.dat$age.forest,  main="Forest age, canopy area adj.", xlab="Age, yrs", xlim=c(0, 300), breaks=40)
hist(biom.dat$age.perv,  main="Forest age, pervious area adj.", xlab="Age, yrs", xlim=c(0, 300), breaks=40)
biom.dat[is.finite(age.ground), length(age.ground)]/biom.dat[aoi>800, length(aoi)] ## 97% retrieval
biom.dat[is.finite(age.forest), length(age.forest)]/biom.dat[aoi>800, length(aoi)] ## 71% retrieval
biom.dat[is.finite(age.perv), length(age.perv)]/biom.dat[aoi>800, length(aoi)] ## 65% retreival
### so you have different total numbers with the different methods -- beware then when you are comparing freq distributions (i.e. gap fill all you can)

### calculate the npp for each method
### coefficients for growth equation
a=123.67
b=0.04

### I don't like the idea of treating the age=NA pix as if they were unknown -- they are NA because they don't retrieve good age, which means they are "old" and npp=0
## figure out next annual increment possible for the "forest" average present in each cell, based on projected "age" and corrected for area
biom.dat[,npp.ann.ground:=((a*(1-(exp(-b*(age.ground+1))))^3)-(a*(1-(exp(-b*(age.ground))))^3))*(1E-4)*aoi] ## by ground area
biom.dat[,npp.ann.forest:=((a*(1-(exp(-b*(age.forest+1))))^3)-(a*(1-(exp(-b*(age.forest))))^3))*(1E-4)*(aoi*can.frac)] ## by canopy area
biom.dat[,npp.ann.perv:=((a*(1-(exp(-b*(age.perv+1))))^3)-(a*(1-(exp(-b*(age.perv))))^3))*(1E-4)*(aoi*(1-isa.frac))] ## by pervious area
#### convert these to kg-biomass gain rather than MgC gain
biom.dat[, npp.ann.ground:=npp.ann.ground*1000*2]
biom.dat[, npp.ann.forest:=npp.ann.forest*1000*2]
biom.dat[, npp.ann.perv:=npp.ann.perv*1000*2]

hist(biom.dat[,npp.ann.ground]) ## OK
hist(biom.dat[,npp.ann.forest])
hist(biom.dat[,npp.ann.perv])

### clean up artifacts
summary(biom.dat$npp.ann.ground)
summary(biom.dat$npp.ann.forest)
summary(biom.dat$npp.ann.perv) ## a lot more NAs in the forest and perv
biom.dat[is.finite(live.MgC.ha.ground) & !is.finite(age.ground) & aoi>800, min(live.MgC.ha.ground)] ### anything 123.7 and above fails to retrieve age+npp
biom.dat[is.finite(live.MgC.ha.forest) & !is.finite(age.forest) & aoi>800, min(live.MgC.ha.forest)] ### anything 123.7 and above fails to retrieve age
biom.dat[is.finite(live.MgC.ha.perv) & !is.finite(age.perv) & aoi>800, min(live.MgC.ha.perv)] ### anything 123.7 and above fails to retrieve age

par(mfrow=c(3,1))
plot(biom.dat$live.MgC.ha.ground, biom.dat$npp.ann.ground, xlim=c(0,200))
plot(biom.dat$live.MgC.ha.forest, biom.dat$npp.ann.forest, xlim=c(0,200))
plot(biom.dat$live.MgC.ha.perv, biom.dat$npp.ann.perv, xlim=c(0,200)) ## different, but all cut out ~123 MgC/ha

### assign 0 npp to all super-high biomass cells
biom.dat[live.MgC.ha.ground>123.6, npp.ann.ground:=0]
biom.dat[live.MgC.ha.forest>123.6, npp.ann.forest:=0]
biom.dat[live.MgC.ha.perv>123.6, npp.ann.perv:=0]
# summary(biom.dat$npp.ann.ground)
# summary(biom.dat$npp.ann.forest)
# summary(biom.dat$npp.ann.perv) ## a handful of extra NAs in perv

# View(biom.dat[is.finite(npp.ann.ground) & !is.finite(npp.ann.perv),]) ## all partial pix with NA isa, fine
# biom.dat[aoi>800 & is.na(npp.ann.ground),] #962 non retreivs, all missing biomass
# biom.dat[aoi>800 & is.na(npp.ann.forest),] #962 non retreivs, all missing biomass
# biom.dat[aoi>800 & is.na(npp.ann.perv),] #972 non retreivs, all missing biomass
# View(biom.dat[aoi>800 & is.na(npp.ann.perv) & !is.na(npp.ann.ground),]) #972 non retreivs, all missing biomass

## fix for all biomass==0
biom.dat[bos.biom30m==0, npp.ann.ground:=0]
biom.dat[bos.biom30m==0, npp.ann.forest:=0]
biom.dat[bos.biom30m==0, npp.ann.perv:=0] ## good enough, have retrievals for almost everything consistently

### look at some plots 
plot(biom.dat$npp.ann.ground, biom.dat$npp.ann.forest)
plot(biom.dat$npp.ann.ground, biom.dat$npp.ann.perv) ## nothing ever exceeds the ground figure
par(mfrow=c(3,1))
hist(biom.dat$npp.ann.ground, main="NPP, raw area", xlab="MgC/pix/yr", breaks=40)
hist(biom.dat$npp.ann.forest, main="NPP, canopy area", xlab="MgC/pix/yr", breaks=40)
hist(biom.dat$npp.ann.perv, main="NPP, pervious area", xlab="MgC/pix/yr", breaks=40)

### aggregated stats
biom.dat[,sum(npp.ann.ground, na.rm=T)]/(2*1000) #13.8k tC/yr by raw ground area
((biom.dat[,sum(npp.ann.ground, na.rm=T)]/(2*1000))/biom.dat[,sum(aoi, na.rm=T)])*1E4 ### 1.1 MgC/ha/yr raw ground area
biom.dat[,sum(npp.ann.forest, na.rm=T)]/(2*1000) #4.6k tC/yr by canopy area
((biom.dat[,sum(npp.ann.forest, na.rm=T)]/(2*1000))/biom.dat[,sum(aoi, na.rm=T)])*1E4 ### 0.37 MgC/ha/yr canopy area area
biom.dat[,sum(npp.ann.perv, na.rm=T)]/(2*1000) #5.4k tC/yr by pervious area
((biom.dat[,sum(npp.ann.perv, na.rm=T)]/(2*1000))/biom.dat[,sum(aoi, na.rm=T)])*1E4 ### 0.43 MgC/ha/yr pervious area
## contrast 10.3-8.9 = 1.4 NEP for boston region in Hardiman

### age distribution
hist(biom.dat$age.ground)
hist(biom.dat$age.forest)
hist(biom.dat$age.perv)
biom.dat[, median(age.ground, na.rm = T)] ##20.2
biom.dat[, median(age.forest, na.rm = T)] ##39.7
biom.dat[, median(age.perv, na.rm = T)] ##37.3

write.csv(biom.dat, "processed/npp.FIA.v3.csv")

### A SLIGHT TWEAK (not a big or systematic effect apparently)
# ### Applying different FIA coefficients for different forest types to estimate 30m annual NPP (MgC/yr)
# 
# library(data.table)
# library(raster)
# 
# ## cleaned up code, handles different growth parameters for different forests (note "Afforestation" and "Reforestation" values are same viz. live biomass growth rates)
# ## initialize parameters for different forest types that ?? resemble tree species distributions in Boston
# for.type <- c("NEdefault","Mixoak", "Redoak")
# a <- c(123.67, 130.81, 123.09)
# b <- c(0.04, 0.03, 0.04)
# 
# ## test limits of the live biomass~stand age function
# for(f in 1:length(for.type)){
#   x=seq(0,120)
#   liveC=a[f]*(1-(exp(-b[f]*x)))^3
#   plot(x,liveC, main=for.type[f], ylab="live tC/ha", xlab="stand age")
#   x <- seq(0, 120) ## inverse: model age vs. live biomass
#   st.age=log(1-(x/a[f])^(1/3))/(-b[f]) ##
#   plot(x, st.age, main=for.type[f], ylab="live tC/ha", xlab="stand age")
# }
# diff(st.age) ## lagged differences --> yearly increment in C gain
# diff(liveC)
# ## conservatively, none of the models is particularly stable over 100 yrs stand age
# 
# biom.dat[, delta.npp.forest:=npp.ann.forest-npp.ann.ground]
# biom.dat[, delta.npp.perv:=npp.ann.perv-npp.ann.ground]
# 
# ## package up some summary rasters
# biom.dat.r <- biom.dat
# for(g in 5:19){
#   r <- raster(biom)
#   r <- setValues(r, biom.dat.r[[g]])
#   writeRaster(r, filename=paste("processed/boston/fia/fia", colnames(biom.dat)[g], "tif", sep="."),
#               format="GTiff", overwrite=T)
# }
#####