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


#### Read in FIA data queried by Luca.
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
# View(live[lag==0 & order(uniqueID),])
# orph <- live[lag==0 & DIAM_T0>0, uniqueID] ## 13 orphans with prior records but no first-instance record
# live <- live[!(uniqueID%in%orph), ] ### kill the orphans, down to 15724
# live[, length(unique(uniqueID))] ## 6384 stems with DBH histories 
# length(live[,unique(PlotID)])## 249 unique plots
# a <- live[PlotID==231,]
# View(a[order(TreeID),])

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


### Modeling growth~dbh (individual stem level) ## this is purely out of curiosity & need for a comparable figure panel
#####
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

zzz <- lmer(diam.rate~DIAM_T0 +
              (1+DIAM_T0|PlotID) +
              (1+DIAM_T0|GENUS)+
              (1+DIAM_T0|YEAR), data=live[lag>0 & STATUS ==1,], REML=F)
yyy.null <- lmer(diam.rate~1 +
              (1|PlotID) +
              (1|GENUS)+
              (1|YEAR), data=live[lag>0 & STATUS ==1,], REML=F)

anova(yyy.null, yyy)

save(yyy, file="processed/mod.fia.final.sav")

# 
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
### equivalent plot-level assessment in live-only data ## 331 plots have at least one subplot dense enough to bother with
live <- as.data.table(read.csv("processed/fia.live.stem.dbh.growth.csv")) ## this is data groomed to remove artifacts and subplots with <5 stems

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
## remove/account for sparse subplots
## affix the records-based lag for annualizing change
### this approach doesn't pay attention to the fate of individual stems -- these may come and go, but we know how much biomass is there at each sampling event and how many subplots we've included

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
# hist(plot.sum$HWfrac) ## most are high HW fraction
# hist(plot.sum$num.subplots) ## about evenly split 1-4

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


### for some reason I thought this was the right way to look at the data but the above seems to show it doesn't have to be this complicated.# 
# # subs <- live[, length(DIAM_T0), by=.(PlotID, SubID, YEAR)]
# # hist(subs$V1)
# # subs[V1<=5,] ## 823 subplot samples out of 2214 less than 5
# # subs[,max(V1)]
# plot.track <- integer()
# subs.track <- integer()
# events.track <- integer()
# year0.track <- integer()
# year1.track <- integer()
# year2.track <- integer()
# biom0HW.track <- numeric()
# biom0SW.track <- numeric()
# biom1HW.track <- numeric()
# biom1SW.track <- numeric()
# biom2HW.track <- numeric()
# biom2SW.track <- numeric()
# stems0.track <- integer()
# stemsBAD1.track <- integer()
# stemsBAD2.track <- integer()
# 
# ## do this intermediate processing step for every plot for up to two resampling events
# vets <- live[YEAR<2007, unique(PlotID)] ## 201 plots started before 2007
# ### OR: Run everything, culling for shit too young to be remeasured, and account for sampling window with starting year
# vets <- live[, unique(PlotID)] ## ~300 plots
# for(p in 1:length(vets)){
#   events <- live[PlotID==vets[p], unique(YEAR)] ## when measured
#   events <- events[order(events)] ## these have to be chronological
#   if(length(events)>1){ ## do not process if visited only once
#     plot.track <- c(plot.track, vets[p]) ## record this plot as having data
#     events.track <- c(events.track, length(events)) ## flag how many measurement events there are for this plot
#     subs <- live[PlotID==vets[p] & YEAR==min(events), length(DIAM_T1), by=SubID]
#     keep <- subs[V1>=5, SubID] ## which subplots have at least 5 trees
#     year0.track <- c(year0.track, min(events)) ## track year0 for this plot
#     subs.track <- c(subs.track, length(keep)) ## track how many subs you're watching to start
#     biom0HW.track <- c(biom0HW.track, live[PlotID==vets[p] & YEAR==min(events) & SubID%in%keep & type=="H", sum(biom1.spp)])
#     biom0SW.track <- c(biom0SW.track, live[PlotID==vets[p] & YEAR==min(events) & SubID%in%keep & type=="S", sum(biom1.spp)])
#     stems0.track <- c(stems0.track, live[PlotID==vets[p] & YEAR==min(events) & SubID%in%keep, length(DIAM_T0)]) ## how many stems initially present?
#     for(i in 2:length(events)){
#       assign(paste("year", (i-1), ".track", sep=""), c(get(paste("year", (i-1), ".track", sep="")), events[i])) ## track year0 for this plot
#       assign(paste("biom", (i-1), "HW.track", sep=""), c(get(paste("biom", (i-1), "HW.track", sep="")), live[PlotID==vets[p] & YEAR==events[i] & SubID%in%keep & type=="H", sum(biom1.spp)])) 
#       assign(paste("biom", (i-1), "SW.track", sep=""), c(get(paste("biom", (i-1), "SW.track", sep="")), live[PlotID==vets[p] & YEAR==events[i] & SubID%in%keep & type=="S", sum(biom1.spp)]))
#       assign(paste("stemsBAD", (i-1), ".track", sep=""), c(get(paste("stemsBAD", (i-1), ".track", sep="")), live[PlotID==vets[p] & YEAR==events[i] & SubID%in%keep & STATUS%in%c(2,3), length(DIAM_T0)])) 
#     }
#     # live[PlotID==vets[p] & YEAR==events[i] & SubID%in%keep,]
#   } 
# }
# ## find a plot that gets measured a lot and has morts/removals
# lots <- arboles[num.measures==3, PlotID]
# live[PlotID %in% lots & STATUS%in%c(2,3),]
# View(live[PlotID==60 & STATUS%in%c(2,3),]) ## a bunch of small birches died
# ### now constitute this into a long-form biomass difference tracker
# subplot.area <- (7.3152^2)*pi
# dat <- data.frame(Year=year0.track)
# dat$PlotID <- plot.track
# dat$HWfrac <- biom0HW.track/(biom0HW.track+biom0SW.track)
# dat$stems0 <- stems0.track
# dat$GONE1 <- stemsBAD1.track
# dat$tot.biom0 <- (biom0HW.track+biom0SW.track)
# dat$biom0.MgC.ha <- ((dat$tot.biom0/2000)/(subs.track*subplot.area))*1E4
# dat$delta.biomHW <- (biom1HW.track-biom0HW.track)
# ## recall that in andy forest we used a model for *annual* dbh gain applied to the stems to estimate biomass gain
# ## but here, we just have biomass at two lagged times, separated by a variable lag -- we have to annualize by lag
# dat$delta.biomHW.ann <- dat$delta.biomHW/(year1.track-year0.track)
# dat$biomHW.gain.rel <- dat$delta.biomHW.ann/dat$tot.biom0
# dat <- as.data.table(dat)
# plot(dat$biom0.MgC.ha, dat$biomHW.gain.rel)
# points(dat[HWfrac<0.25, biom0.MgC.ha], dat[HWfrac<0.25, biomHW.gain.rel], pch=16, col="green")
# points(dat[GONE1/stems0>0.25, biom0.MgC.ha], dat[GONE1/stems0>0.25, biomHW.gain.rel], pch=16, col="red")
# # dat[biomHW.gain.rel<(-0.10),] ### INCLUDING the morts etc. there is just one plot that loses >10% biomass
# # View(live[PlotID==136,]) ### this plot gets logged!
# # live[PlotID==136, length(DIAM_T1), by=YEAR] ## 34 trees in 2004 to just 10 in 2013, and we watch the same subplots and don't restrict the stems we allow in
# hist(dat$HWfrac) ## only a fraction are low HWfrac
# rocked.tab <- dat[GONE1/stems0>0.25, PlotID]
# 
# 
# ### load up the second measurement round (119 additional plot measurements)
# dat2 <- as.data.table(data.frame(Year=year1.track))
# dat2[,PlotID:=plot.track]
# dat2[,subs:=subs.track]
# dat2[,events:=events.track]
# dat2[,biom0HW:=biom1HW.track] ## treat time 1 as time 0 in the repeated measurement
# dat2[,tot.biom0:=(biom1HW.track+biom1SW.track)]
# dat2[,HWfrac:=biom0HW.track/(biom0HW.track+biom0SW.track)] ## original time0 HW fraction
# dat2[,stems0:=stems0.track]
# dat2[,biom0.MgC.ha:=((tot.biom0/2000)/(subs.track*subplot.area))*1E4]
# dat2 <- dat2[events.track==3,] ## now we match length with year2.track
# dat2[,GONE1:=stemsBAD2.track]
# dat2[,delta.biomHW:=(biom2HW.track-biom0HW)]
# dat2[,delta.biomHW.ann:=delta.biomHW/(year2.track-Year)]
# dat2[,biomHW.gain.rel:=delta.biomHW.ann/tot.biom0]
# ### stems marked removed in one year do not appear in next but morts do?
# stems1 <- dat[PlotID%in%dat2[,unique(PlotID)], stems0-GONE1]
# dat2[,stems0:=stems1] ## recalibrate your baseline of stem numbers
# plot(dat2$biom0.MgC.ha, dat2$biomHW.gain.rel) ## same pattern, some serious losses
# points(dat2[GONE1/stems0>0.25, biom0.MgC.ha], dat2[GONE1/stems0>0.25, biomHW.gain.rel], col="red", pch=16)
# dat2$subs <- NULL
# dat2$events <- NULL
# dat2$biom0HW <- NULL
# dat2 <- setcolorder(dat2, names(dat))
# 
# dat <- rbind(dat, dat2) ## rocking 375 records now
# dat[,rocked:=0]
# dat[PlotID%in%rocked.tab, rocked:=1] ## flag plots that showed high stem loss early on, 41 plots early on
# plot(dat$biom0.MgC.ha, dat$biomHW.gain.rel) ## hyperbolic but mostly flat, say 20 serious biomass loss plots
# points(dat[HWfrac<0.25, biom0.MgC.ha], dat[HWfrac<0.25, biomHW.gain.rel], pch=16, col="green")
# points(dat[rocked==1, biom0.MgC.ha], dat[rocked==1, biomHW.gain.rel], pch=16, col="red")
# points(dat[GONE1/stems0>0.25, biom0.MgC.ha], dat[GONE1/stems0>0.25, biomHW.gain.rel], pch=16, col="purple")
# dat[,PlotID:=as.factor(PlotID)]
# dat[,Year:=as.factor(Year)]
# plot(dat[!PlotID %in% c(47,283,183,5), biom0.MgC.ha], dat[!PlotID %in% c(47,283,183,5), biomHW.gain.rel]) ## hyperbolic but mostly flat, say 20 serious biomass loss plots
# points(dat[!PlotID %in% c(47,283,183,5) & rocked==1, biom0.MgC.ha], dat[!PlotID %in% c(47,283,183,5) & rocked==1, biomHW.gain.rel], col="red", pch=16) ## hyperbolic but mostly flat, say 20 serious biomass loss plots
# dat[rocked==0 & biomHW.gain.rel<(-0.05),PlotID] ## these are the weird plots that get stupid in later sample events
# 
# runme <- unique(dat[rocked==0 & !PlotID%in%c(47,283,183,5) & HWfrac>=0.25, PlotID]) ## 184 unique plots


## OK, figure a model and coefficient error terms
live.plot <- as.data.table(read.csv("processed/fia.live.plot.groomedV2.csv"))
summary(live.plot[,YEAR-prev.sample])

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
plot(live.plot[HWfrac>0.25, biom0.MgC.ha], live.plot[HWfrac>0.25, biom.delt.ann.rel.HW])
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

#######
### previous analysis not paying attention to stem or plot repeat visits
live.plot <- live[,
                  .(length(unique(SubID)),
                    sum(ba.pred(DIAM_T0)),
                    sum(biom.delt.spp, na.rm=T), 
                    sum(biom0.spp, na.rm=T),
                    length(DIAM_T0)),
                  by=PlotID]
# dim(live.plot) ## 221 plots with at least one valid subplot
names(live.plot) <- c("PlotID", "num.subplots", "total.BA", "biom.growth", "total.biom0.kg", "num.stems")
plot(live.plot[,num.stems], live.plot[,total.biom0.kg]) ## linear but gets unconstrained above ~20 stems/plot
# cor(live.plot[, total.biom0.kg], live.plot[, num.stems]) ## rho 0.60
# live.plot[total.biom0.kg<5000,] ## mostly single-subplot plots (ie plot is mostly non-forest)

## plot-level growth
subplot.area <- (7.3152^2)*pi ## subplots are 24ft in radius
live.plot[,total.biom0.MgC.ha:=(total.biom0.kg/2000)/((num.subplots*subplot.area)/1E4)] ## biomass density in aggregated area covered by the fully forested subplots
live.plot[,biom.growth.ann:=biom.growth/4.8]
live.plot[,biom.growth.ann.rel:=(biom.growth.ann)/total.biom0.kg]
# hist(live.plot$biom.growth.ann.rel) ## most plots below 5%, a few are wildly productive, a couple decline
# hist(live.plot$BA.plot) ## ok BA is looking like forest mostly -- 20-40
# summary(live.plot$biom.growth.ann.rel) ## -1% to 14%

# ## plot level density
live.plot[,stems.ha:=(num.stems/(subplot.area*num.subplots))*1E4]
live.plot[,BA.plot:=total.BA/(num.subplots*subplot.area/1E4)]
# hist(live.plot$BA.plot)## peak below 40m2/ha, a few up to 80
# summary(live.plot$ba.ha)
# hist(live.plot$stems.ha) # 400-500 peak
# hist(live.plot[, total.biom0.kg], breaks=10) ## peak 5-10k kg
# hist((live.plot[, total.biom0.kg]/675)*1E4/2000, breaks=10) ## i.e. peaking 50-100 MgC/ha, up to 200
# summary(live.plot[, total.biom0.kg/(subplot.area*num.subplots)*1E4]/2000) ## 90 MgC/ha, 34-307 -- so a forest
# hist(live.plot$num.stems)
# hist(live.plot$num.subplots) ## pretty even distribution, about 1/4 in each class from 1 to 4 subplots available



# ## basic growth~biomass modeling
# plot(live.plot$total.biom0.kg, live.plot$biom.growth.ann.rel) ## exponential negative, the super productive ones were very low biomass to begin
# plot(live.plot$total.biom0.kg, live.plot$biom.growth.ann) ## linear-ish, more biomass --> more growth (no bottoming out)
# plot(live.plot$total.biom0.MgC.ha, live.plot$biom.growth.ann.rel) ## OK more happily exponential
# plot(live.plot[,log(total.biom0.MgC.ha)], live.plot[,log(biom.growth.ann.rel)]) ## here's our regular log-log linear relationship, pretty shotgunned
# 
# l <- lm(log(biom.growth.ann.rel)~log(total.biom0.Mg.ha), data=live.plot) ### R2 0.17, significant
# plot(live.plot[,log(total.biom0.Mg.ha)], live.plot[,log(biom.growth.ann.rel)])
# s <- summary(l)
# abline(a=s$coefficients[1], b=s$coefficients[2], col="red") ## ok ho hum
# ## can do this as a proper expnential nls to also incorporate the negative growth plots


## look at specific growth rates in the hardwood/softwood population separate from total growth in all trees
## spp in these names indicates that species-specific allometries were applied to the underlying stem data
hwood <- live[type=="H", .(sum(biom1.spp-biom0.spp, na.rm=T), 
                           sum(biom0.spp, na.rm=T)), 
              by=PlotID]
names(hwood)[2:3] <- c("biom.growth.spp.hw", "total.biom0.spp.kg.hw")
swood <- live[type=="S", .(sum(biom1.spp-biom0.spp, na.rm=T), 
                           sum(biom0.spp, na.rm=T)),
              by=PlotID]
names(swood)[2:3] <- c("biom.growth.spp.sw", "total.biom0.spp.kg.sw")

live.plot <- merge(live.plot, hwood, by="PlotID")
live.plot <- merge(live.plot, swood, by="PlotID")
live.plot[,biom.growth.ann.hw:=(biom.growth.spp.hw/4.8)/total.biom0.spp.kg.hw]
live.plot[,biom.growth.ann.sw:=(biom.growth.spp.sw/4.8)/total.biom0.spp.kg.sw]
live.plot[,hw.frac:=total.biom0.spp.kg.hw/total.biom0.kg] ## get fraction of hardwood in total biomass
write.csv(live.plot, "processed/fia.live.plot.groomed.csv")

live.plot <- as.data.table(read.csv("processed/fia.live.plot.groomed.csv"))
# ### here follows a lot of hemming and hawing until the final model, excluding low hw.frac, is discovered.
# ## hard wood and soft wood growth in stands, identifying stands that a low in either
# plot(live.plot[,total.biom0.MgC.ha], live.plot[,biom.growth.ann.hw]) ## ding dong goofball with big outliers, all low hw
# points(live.plot[hw.frac<0.25, total.biom0.MgC.ha], live.plot[hw.frac<0.25, biom.growth.ann.hw], cex=1.5, col="red", pch=14)
# 
# plot(live.plot[,total.biom0.MgC.ha], live.plot[,biom.growth.ann.sw]) ## low softwood frac is the norm, but more well behaved
# points(live.plot[hw.frac>0.75, total.biom0.MgC.ha], live.plot[hw.frac>0.75, biom.growth.ann.sw], cex=1.5, col="green", pch=14)
# 
# ### comparative models in hardwoods and softwoods in isolation
# hw <- lm(live.plot[biom.growth.ann.hw>0,log(biom.growth.ann.hw)]~live.plot[biom.growth.ann.hw>0,log(total.biom0.MgC.ha)])
# summary(hw) ## barely significant (p 0.07), R2 0.01
# sw <- lm(live.plot[biom.growth.ann.sw>0,log(biom.growth.ann.sw)]~live.plot[biom.growth.ann.sw>0,log(total.biom0.MgC.ha)])
# summary(sw) ## significant, R2 0.12
# 
# ### model without respect to hard/softwood
# plot(live.plot[,log(total.biom0.MgC.ha)], live.plot[,log(biom.growth.ann.rel)])
# points(live.plot[hw.frac<0.25,log(total.biom0.MgC.ha)], live.plot[hw.frac<0.25,log(biom.growth.ann.rel)], col="red", pch=14, cex=1.6)
# tot.mod <- lm(live.plot[biom.growth.ann.rel>0,log(biom.growth.ann.rel)]~live.plot[biom.growth.ann.rel>0,log(total.biom0.MgC.ha)])
# summary(tot.mod) ### ok does like R2 0.18 and significant
# 
# plot(live.plot[,log(total.biom0.MgC.ha)], live.plot[,log(biom.growth.ann.rel)], 
#      main="FIA growth~biomass (plot level)", ylab="Relative growth rate (Mg/Mg/ha)", xlab="Plot biomass density (Mg/ha)",
#      xaxt="n", yaxt="n")
# axis(side = 1, at = c(4.5, 5, 5.5, 6, 6.5), labels = round(exp(c(4.5, 5, 5.5, 6, 6.5)), digits=0))
# axis(side=2, at=seq(-5.5, -2.5, by=0.5), labels=round(exp(seq(-5.5, -2.5, by=0.5)), digits=3))
# abline(a=tot.mod$coefficients[1], b=tot.mod$coefficients[2], lwd=1.4, col="black")
# # points(live.plot[,log(((total.biom0.spp.kg.hw/1000)/(num.subplots*subplot.area))*1E4)],
# #        live.plot[,log(biom.growth.ann.hw)], pch=15, cex=0.6, col="red")
# points(live.plot[,log(total.biom0.Mg.ha)],
#        live.plot[,log(biom.growth.ann.hw)], pch=15, cex=0.6, col="red")
# abline(a=hw$coefficients[1], b=hw$coefficients[2], col="red")
# # points(live.plot[,log(((total.biom0.spp.kg.sw/1000)/(num.subplots*subplot.area))*1E4)],
# #        live.plot[,log(biom.growth.ann.sw)], pch=17, cex=0.6, col="blue")
# points(live.plot[,log(total.biom0.Mg.ha)],
#        live.plot[,log(biom.growth.ann.sw)], pch=17, cex=0.6, col="blue")
# abline(a=sw$coefficients[1], b=sw$coefficients[2], col="blue")
# legend(x=4.3, y=-5, legend = c("All", "Hardwoods", "Softwoods"), 
#        fill = c("black", "red", "blue"))
# summary(live.plot$biom.growth.ann.hw) ## hardwood growth rate just doesn't vary much across range of TOTAL biomass density: 2-3% growth (0.4 to 15%)
# 
# 
# #### let's look at it in untransformed space
# tot.mod <- lm(live.plot[biom.growth.ann.rel>0,log(biom.growth.ann.rel)]~live.plot[biom.growth.ann.rel>0,log(total.biom0.Mg.ha)])
# summary(tot.mod) ### R2 0.18 and significant
# hw.mod <- lm(live.plot[biom.growth.ann.hw>0,log(biom.growth.ann.hw)]~live.plot[biom.growth.ann.hw>0,log(total.biom0.Mg.ha)])
# summary(hw.mod) ### R2 0.01 and only sig at p<0.10
# hw.mod$coefficients
# tot.mod$coefficients ## slope in hw.mod is very low, not particularly high in tot.mod
# 
# ## untransformed
# hw.mod.exp <- nls(biom.growth.ann.hw ~ exp(a + b * log(total.biom0.Mg.ha)),
#                   data=live.plot, start=list(a=0, b=0))
# q <- summary(hw.mod.exp) ## b is not even close to significant
# exp(q$coefficients[1]) ## ie. about 3% boyo
# 
# tot.mod.exp <- nls(biom.growth.ann.rel ~ exp(a + b * log(total.biom0.Mg.ha)),
#                    data=live.plot, start=list(a=0, b=0)) ### not terrifically differfent from final HW-frac filtered model below
# r <- summary(tot.mod.exp) ## definitely significant
# 
# ## visual comparison of the plot-level models
# plot(live.plot[,(total.biom0.Mg.ha)], 
#      live.plot[,(biom.growth.ann.rel)], main="FIA growth~biomass (plot-level)",
#      xlab="Biomass Density (Mg/ha)", ylab="Relative growth (Mg/Mg/ha)")
# points(live.plot[,(total.biom0.Mg.ha)], 
#        live.plot[,(biom.growth.ann.hw)], col="red", pch=17, cex=0.6)
# test <- seq(min(live.plot$total.biom0.Mg.ha), max(live.plot$total.biom0.Mg.ha), length.out=100)
# lines(test, exp(r$coefficients[1]+(r$coefficients[2]*log(test))),
#       lwd=3, lty=2)
# # abline(h=mean(live.plot$biom.growth.ann.hw), col="orangered", lwd=3, lty=4)
# 
# ## why so poor fit on hw model, why some high biomass plots have such high hw growth rate?
# # summary(live.plot$hw.frac) ## most have a majority hardwoods, a few v low
# # hist(live.plot$hw.frac)
# # live.plot[total.biom0.Mg.ha>400, hw.frac] ## high biomass contains the low hw fraction ones
# points(live.plot[hw.frac<0.25, total.biom0.Mg.ha],
#        live.plot[hw.frac<0.25, biom.growth.ann.hw],
#        cex=1.2, lwd=3, pch=5, col="seagreen")

### OOOOOOKKKKKKAAAAYYYYY Let's elminate the handful of weird low-HW plots (could be weird places that favor pines or something, not really how an urban forest do)
hw.mod.exp.filt <- nls(biom.growth.ann.hw ~ exp(a + b * log(total.biom0.MgC.ha)),
                       data=live.plot[hw.frac>0.25,], start=list(a=0, b=0))
y <- summary(hw.mod.exp.filt) ## all significant

## so we will use the low-hw filtered plots to approximate growth in the final NPP calculations
plot(live.plot[,total.biom0.MgC.ha], live.plot[,biom.growth.ann.hw], pch=17, col="red", cex=0.6)
points(live.plot[hw.frac<0.25, total.biom0.MgC.ha], live.plot[hw.frac<0.25, biom.growth.ann.hw], pch=5, col="seagreen", cex=0.8)
test <- seq(min(live.plot$total.biom0.MgC.ha), max(live.plot$total.biom0.MgC.ha, length.out=100))
lines(test, exp(y$coefficients[1]+(y$coefficients[2]*log(test))),
      lwd=3, lty=2, col="black")
legend(x=100, y=0.01, legend=c("Hardwoods", "Low HW frac.", "All"), 
       fill=c("red",  "seagreen", "black"), bty="n")




########
##### PART B: BOSTON NPP FROM FIA EMPRICAL GROWTH~BIOMASS ANALYSIS
######
### V2.2: 1) uses species-specific biometrics; 2) models hardwood growth separately from trees in general; 3) uses nls to avoid dumping negative growth records
### V2.3 1) Uses subplot IDs to remove dbh records from subplot sites that are not fully forested; 2) filtered plots that have low hw fraction to determine hw-only growth rate
### V2.4 1) uses all-status tree data, get bulk biomass by plot in each sample event; 2) remove any subplots that do not maintain at >=5 live stems in all sample events; 
### 3) filter plots with <25% HW biomass; 4) simple LMER with random intercepts for starting year and PlotID (some plots are revisited more than once)
### 
##
# live.plot <- as.data.table(read.csv("processed/fia.live.plot.groomed.csv")) ## previous data set not tracking morts/removals or plot/stem multiple records appearances
library(lme4)
live.plot <- as.data.table(read.csv("processed/fia.live.plot.groomedV2.csv")) ## newest thing
# runme <- unique(live.plot[rocked==0 & !PlotID%in%c(47,283,183,5) & HWfrac>=0.25, PlotID]) ## 184 unique plots

# mod.live.plot.final <- nls(biom.growth.ann.hw ~ exp(a + b * log(total.biom0.MgC.ha)),
#                            data=live.plot[hw.frac>0.25,], start=list(a=0, b=0))
# summary(mod.live.plot.final) ## all significant
# load("processed/mod.live.plot.final.sav") ## direct from npp_anlysis.R
# 
# mod.live.plot.final <- lmer(biomHW.gain.rel~poly(biom0.MgC.ha, degree=2, raw=T)+
#             (1|PlotID)+
#             (1|Year),
#           data=live.plot[PlotID%in%runme,], REML=F)
### 4th degree fits slightly better but leads to unstable growth factor prediction when you vary the coefficients by a bit

# hw.mod.exp.filt <- mod.live.plot.final
# y <- summary(hw.mod.exp.filt) ## residual standard error = standard error of regression = how far off values may be from predicted (vs. R2, which is unreliable)

## model using biomass gain MgC/ha/yr
# mod.live.plot.final <- lmer(biom.delt.ann.HW.MgC.ha~biom0.MgC.ha+
#                               (1|PlotID)+
#                               (1|prev.sample), 
#                             data=live.plot[HWfrac>0.25,])

## model using biomass gain MgC/MgC logxform
## the below is model a.log above
mod.live.plot.final <- lmer(log(biom.delt.ann.rel.HW)~biom0.MgC.ha+
                (1|PlotID)+
                (1|prev.sample),
              data=live.plot[HWfrac>0.25], REML=F)

plot(live.plot[,biom0.MgC.ha], live.plot[,biom.delt.ann.rel.HW])
points(live.plot[HWfrac<=0.25,biom0.MgC.ha], live.plot[HWfrac<=0.25,biom.delt.ann.rel.HW], pch=12, col="red")


y <- summary(mod.live.plot.final)
save(mod.live.plot.final, file="processed/mod.live.plot.final.sav")

## load the biomass data and reprocess
biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m kg-biomass to 30m pixel
aoi <- raster("processed/boston/bos.aoi30m.tif")
can <- raster("processed/boston/bos.can30m.tif")
isa <- raster("processed/boston/bos.isa30m.tif") ## this is the newly reregistered guy
biom <- crop(biom, aoi)
can <- crop(can, aoi)
isa <- extend(isa, aoi)
lulc <- raster("processed/boston/bos.lulc30m.lumped.tif")
biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[,aoi:=as.vector(getValues(aoi))]
can.dat <- as.data.table(as.data.frame(can))
biom.dat[, can.frac:=can.dat$bos.can30m]
isa.dat <- as.data.table(as.data.frame(isa))
biom.dat[, isa.frac:=isa.dat$bos.isa30m]
lulc.dat <- as.data.table(as.data.frame(lulc))
biom.dat[,lulc:=lulc.dat$bos.lulc30m.lumped]
biom.dat[,pix.ID:=seq(1, dim(biom.dat)[1])]


### Now we have to decide how to calculate the "forest" density of the biomass on the ground
### Three approaches: GROUND (kg-biomass/m2-pixel); "FOREST" (kg-biomass/m2-canopy); "PERV" (kg-biomass/m2-pervious)
## Ground biomass density
biom.dat[, live.MgC.ha.ground:=((bos.biom30m/2000)/aoi)*(1E4)] ## convert kg-biomass/pixel based on size of pixel (summed by 1m2 in "aoi"), 
biom.dat[bos.biom30m<10, live.Mgbiom.ha.ground:=0] ## cancel out tiny biomass cells

## Forest biomass density
biom.dat[,live.MgC.ha.forest:=((bos.biom30m/2000)/(aoi*can.frac))*(1E4)]
biom.dat[bos.biom30m<10, live.MgC.ha.forest:=0] ## have to manually fix this because of 0 canopy pix
biom.dat[can.frac<0.01, live.MgC.ha.forest:=0]

## Perv biomass density
biom.dat[,live.MgC.ha.perv:=((bos.biom30m/2000)/(aoi*(1-isa.frac)))*(1E4)]
biom.dat[bos.biom30m<10, live.MgC.ha.perv:=0] ## have to manually fix this because of isa=1 pix
biom.dat[isa.frac>0.99, live.MgC.ha.perv:=0]

# ### what do our densities look like
# hist(biom.dat[aoi>800, live.MgC.ha.ground]) ## up to 600 Mgbiom/ha = 300 MgC/ha -- a lot of this below the range sampled by the FIA plot data
# summary(biom.dat[aoi>800, live.MgC.ha.ground]) ## mostly below 40 MgC/ha
# hist(biom.dat[aoi>800, live.MgC.ha.forest]) ## same, more medium sized
# summary(biom.dat[aoi>800,live.MgC.ha.forest]) ## everything nudged denser
# hist(biom.dat[aoi>800 & isa.frac<0.98, live.MgC.ha.perv]) ## up to 800k kgbiom/ha
# summary(biom.dat[aoi>800, live.MgC.ha.perv]) ## a few very high
# biom.dat[aoi>800 & live.MgC.ha.perv>300,] ### 3k v high, these are the classic artifacts, v. low pervious cover, but low overall biomass
# ### contrast: what is the range sampled in the live.plot FIA data?
# summary(live.plot[,total.biom0.MgC.ha]) ## don't get into the low range of things! Minimum biomass density is about where most of the pixels fall

# ## root out any artifacts in density calcs
# biom.dat[aoi>800, length(aoi)] #136667 valid pixels
# biom.dat[aoi>800, length(bos.biom30m)] #136667 valid pixels
# biom.dat[aoi>800 & is.finite(bos.biom30m), length(bos.biom30m)] ##135705 is the number to hit
# biom.dat[!is.na(npp.kg.hw.ground) & aoi>800, length(npp.kg.hw.ground)] ## 135705 pix
# biom.dat[!is.na(npp.kg.hw.ground) & aoi>800, length(live.kgbiom.ha.ground)] ## 135705 pix
# biom.dat[is.finite(npp.kg.hw.forest) & aoi>800 & is.finite(bos.biom30m), length(npp.kg.hw.forest)] ## 135705 pix
# biom.dat[is.finite(npp.kg.hw.perv) & aoi>800 & is.finite(bos.biom30m), length(npp.kg.hw.perv)] ## 135695, a few have can but no isa
# biom.dat[is.na(npp.kg.hw.perv) & aoi>800 & is.finite(bos.biom30m),] ## 127736 pix
# ### this is biomass associated with 100% paved pixels, added fix above

## calculate growth factors per cell

### new hotness: build in model estimate uncertainty
dump <- copy(biom.dat)

### loop and interatively make maps with different model realizations
sav.ground <- biom.dat[,1:9, with=F]
sav.forest <- biom.dat[,1:9, with=F]
sav.perv <- biom.dat[,1:9, with=F]
# sav.floor <- numeric()
# sav.floor.implemented <- numeric()
max.fact <- live.plot[,max(biom.delt.ann.rel.HW)]
for(i in 1:100){
  b0.sel <- rnorm(1, coef(y)[1,1], coef(y)[1,2])
  b1.sel <- rnorm(1, coef(y)[2,1], coef(y)[2,2])

  ## assign growth factors based on position and biomass density
  dump[,ground.gfact:=exp(b0.sel+(b1.sel*live.MgC.ha.ground))] ## predict growth factors ground
  dump[(live.MgC.ha.ground)<0.5, ground.gfact:=0] ## anything without less than about 9 kg biomass/cell gets 0 factor
  # dump[ground.gfact<gfact.min, ground.gfact:=gfact.min] ## anything high enough to be past the reliable part of the curve gets stet to min
  dump[,forest.gfact:=exp(b0.sel+(b1.sel*live.MgC.ha.forest))] ## predict growth factors ground
  dump[(live.MgC.ha.forest)<0.5, forest.gfact:=0]
  # dump[forest.gfact<gfact.min, forest.gfact:=gfact.min]
  dump[can.frac<0.01, forest.gfact:=0]
  dump[,perv.gfact:=exp(b0.sel+(b1.sel*live.MgC.ha.perv))] ## predict growth factors ground
  dump[(live.MgC.ha.perv)<0.5, perv.gfact:=0]
  # dump[perv.gfact<gfact.min, perv.gfact:=gfact.min]
  dump[isa.frac>0.99, perv.gfact:=0]
  ### all of these are limited to the minimum, but can reach too high maxima 
  dump[ground.gfact>max.fact, ground.gfact:=max.fact]
  dump[forest.gfact>max.fact, forest.gfact:=max.fact]
  dump[perv.gfact>max.fact, perv.gfact:=max.fact]
  
  
  # ### limit the predictions to something reasonable in this polynomial model
  # test=0:300
  # pred <- b0.sel+(b1.sel*test)+(b2.sel*test^2)+(b3.sel*test^3)+(b4.sel*test^4)
  # # plot(test, pred)
  # gfact.min <- min(pred) ### for when you get into the high biomass stuff
  # sav.floor <- c(sav.floor, gfact.min)
  # if(gfact.min<0){gfact.min <- 0} ## sometimes predict large biomass losses, can't abide this
  # sav.floor.implemented <- c(sav.floor.implemented, gfact.min)
  # ## assign growth factors based on position and biomass density
  # dump[,ground.gfact:=b0.sel+(b1.sel*live.MgC.ha.ground)+(b2.sel*live.MgC.ha.ground^2)+(b3.sel*live.MgC.ha.ground^3)+(b4.sel*live.MgC.ha.ground^4)] ## predict growth factors ground
  # dump[(live.MgC.ha.ground)<0.5, ground.gfact:=0] ## anything without less than about 9 kg biomass/cell gets 0 factor
  # dump[ground.gfact<gfact.min, ground.gfact:=gfact.min] ## anything high enough to be past the reliable part of the curve gets stet to min
  # dump[,forest.gfact:=b0.sel+(b1.sel*live.MgC.ha.forest)+(b2.sel*live.MgC.ha.forest^2)+(b3.sel*live.MgC.ha.forest^3)+(b4.sel*live.MgC.ha.forest^4)] ## predict growth factors ground
  # dump[(live.MgC.ha.forest)<0.5, forest.gfact:=0]
  # dump[forest.gfact<gfact.min, forest.gfact:=gfact.min]
  # dump[can.frac<0.01, forest.gfact:=0]
  # dump[,perv.gfact:=b0.sel+(b1.sel*live.MgC.ha.perv)+(b2.sel*live.MgC.ha.perv^2)+(b3.sel*live.MgC.ha.perv^3)+(b4.sel*live.MgC.ha.perv^4)] ## predict growth factors ground
  # dump[(live.MgC.ha.perv)<0.5, perv.gfact:=0]
  # dump[perv.gfact<gfact.min, perv.gfact:=gfact.min]
  # dump[isa.frac>0.99, perv.gfact:=0]
  
  ## calculate NPP and store
  ### biomass productivity factors
  dump[,npp.kg.hw.ground:=bos.biom30m*ground.gfact] ## apply ground growth factor to whole pixel area
  dump[,npp.kg.hw.forest:=bos.biom30m*forest.gfact] ## apply growth factor to FOREST AREA
  dump[,npp.kg.hw.perv:=bos.biom30m*perv.gfact] ## apply growth factor to PERVIOUS AREA
  
  # ### area productivity factors
  # dump[,npp.kg.hw.ground:=(aoi/1E4)*ground.gfact*2000] ## apply ground growth factor to whole pixel area
  # dump[,npp.kg.hw.forest:=((aoi*can.frac)/1E4)*forest.gfact*2000] ## apply growth factor to FOREST AREA
  # dump[,npp.kg.hw.perv:=((aoi*(1-isa.frac))/1E4)*perv.gfact*2000] ## apply growth factor to PERVIOUS AREA

  sav.ground <- cbind(sav.ground, dump[,npp.kg.hw.ground])
  names(sav.ground)[9+i] <- paste("npp.fia.ground.iter.", i, ".kg", sep="")
  sav.forest <- cbind(sav.forest, dump[,npp.kg.hw.forest])
  names(sav.forest)[9+i] <- paste("npp.fia.forest.iter.", i, ".kg", sep="")
  sav.perv <- cbind(sav.perv, dump[,npp.kg.hw.perv])
  names(sav.perv)[9+i] <- paste("npp.fia.perv.iter.", i, ".kg", sep="")
  print(paste("iteration",i))
}
## save them out
### Version history: 
### V24 = static model 
### V3 = 2nd poly to predict relative growth rate, with model noise;
### V4 = linear model to predict area growth rate (MgC/ha/yr) with model noise
### V5 = log-xform linear model to predict relative growth rate, with model noise
write.csv(sav.ground, "npp.FIA.empirV5.ground.csv")
write.csv(sav.forest, "npp.FIA.empirV5.forest.csv")
write.csv(sav.perv, "npp.FIA.empirV5.perv.csv")
# 
# dump[aoi>800, .(sum(npp.kg.hw.ground, na.rm=T)/2000/1000,
#                 sum(npp.kg.hw.forest, na.rm=T)/2000/1000,
#                 sum(npp.kg.hw.perv, na.rm=T)/2000/1000)]
# hist(dump[aoi>800, ground.gfact])
# hist(dump[aoi>800, forest.gfact])
# hist(dump[aoi>800, perv.gfact])


#####
# ### static, using model coefficients
# biom.dat[,ground.gfact:=exp(y$coefficients[1]+y$coefficients[2]*log(live.MgC.ha.ground))]
# biom.dat[,forest.gfact:=exp(y$coefficients[1]+y$coefficients[2]*log(live.MgC.ha.forest))]
# biom.dat[,perv.gfact:=exp(y$coefficients[1]+y$coefficients[2]*log(live.MgC.ha.perv))]

# ### are we getting reasonable mapping of growth factors?
# summary(biom.dat[aoi>800, ground.gfact]) ## up to inf
# biom.dat[aoi>800 & ground.gfact>0.07,] ## 75k unfiltered are above 7% max measured biomass gain rate
# summary(biom.dat[aoi>800, forest.gfact]) ## better, still max at inf
# biom.dat[aoi>800 & forest.gfact>0.07,] ## 34k above 7%
# summary(biom.dat[aoi>800, perv.gfact]) ## a lot of inf
# biom.dat[aoi>800 & perv.gfact>0.07,] ## 54k above 7%

## kill basic artifacts
# # ((10/2000)/900)*1E4 # 10 kg-biom/pix is 0.06 MgC/ha
# biom.dat[aoi>800 & bos.biom30m<10, ground.gfact:=0]
# biom.dat[aoi>800 & bos.biom30m<10, forest.gfact:=0]
# biom.dat[aoi>800 & can.frac<0.01, forest.gfact:=0]
# biom.dat[aoi>800 & bos.biom30m<10, perv.gfact:=0]
# biom.dat[aoi>800 & isa.frac>0.99, perv.gfact:=0]

# summary(biom.dat[aoi>800, ground.gfact]) ## nax 2,8
# biom.dat[aoi>800 & ground.gfact>0.07,] ## 46k unfiltered are above 7% max measured biomass gain rate
# summary(biom.dat[aoi>800, forest.gfact]) ## max 0.68
# biom.dat[aoi>800 & forest.gfact>0.07,] ## only 917 above 7%
# summary(biom.dat[aoi>800, perv.gfact]) ## max 2.8
# biom.dat[aoi>800 & perv.gfact>0.07,] ## 15kk above 7%

# ## old filter set for v23
# biom.dat[bos.biom30m<10 | live.Mgbiom.ha.forest<1, forest.gfact:=0]
# biom.dat[bos.biom30m<10 | live.Mgbiom.ha.perv<1, perv.gfact:=0]
# biom.dat[bos.biom30m<10 | live.Mgbiom.ha.ground<1, ground.gfact:=0]

### New refinement: We will cap the low biomass pixels at the predicted maximum growth rate at minimum measured biomass
cutoff <- exp(y$coefficients[1]+y$coefficients[2]*log(live.plot[hw.frac>0.25, min(total.biom0.MgC.ha)])) ## cut off is 34 MgC/ha = 0.043

## new filter set for v24
# biom.dat[aoi>800 & ground.gfact>cutoff,] ## 66k are valid with too high growth factor!
# hist(biom.dat[aoi>800 & ground.gfact>cutoff, live.MgC.ha.ground]) ## all are below 35, skewed small
biom.dat[ground.gfact>cutoff, ground.gfact:=cutoff]

# biom.dat[aoi>800 & forest.gfact>cutoff,] ## 17k
# hist(biom.dat[aoi>800 & forest.gfact>cutoff, live.MgC.ha.forest]) ## skewed large
biom.dat[forest.gfact>cutoff, forest.gfact:=cutoff]

# biom.dat[aoi>800 & perv.gfact>cutoff,] ## 29k
# hist(biom.dat[aoi>800 & perv.gfact>cutoff, live.MgC.ha.perv]) ## mostly uniform
biom.dat[perv.gfact>cutoff, perv.gfact:=cutoff]

## final cleanup: remove gfacts for NA AOI
biom.dat[is.na(aoi), ground.gfact:=NA]
biom.dat[is.na(aoi), forest.gfact:=NA]
biom.dat[is.na(aoi), perv.gfact:=NA]

## calculate npp from these growth factors
## regression coeff*biom.density(kg/ha, 3 approaches)-->growth factors (kg/kg) per cell
## growth factors * cell biomass (kg | Mg | MgC) --> NPP (kg biomass per cell, or etc.) 
biom.dat[,npp.kg.hw.ground:=bos.biom30m*ground.gfact] ## this is just the growth factor, can be applied to any biomass units
biom.dat[,npp.kg.hw.forest:=bos.biom30m*forest.gfact]
biom.dat[,npp.kg.hw.perv:=bos.biom30m*perv.gfact]

### totals for aoi
biom.dat[aoi>800 & is.finite(isa.frac) & is.finite(can.frac), sum(npp.kg.hw.ground, na.rm=T)]/2000 ## 9.7k MgC/yr v24; v23 12k MgC/yr ground basis
biom.dat[aoi>800 & is.finite(isa.frac) & is.finite(can.frac), sum(npp.kg.hw.forest, na.rm=T)]/2000 ## 7.8k MgC/yr v24; v23 7.8k MgC/yr forest basis
biom.dat[aoi>800 & is.finite(isa.frac) & is.finite(can.frac), sum(npp.kg.hw.perv, na.rm=T)]/2000  ## 7.4k MgC/yr v24; v23 7.8k MgC/yr perv basis
biom.dat[aoi>800 & is.finite(isa.frac) & is.finite(can.frac), ((sum(npp.kg.hw.ground, na.rm=T)/2000)/sum(aoi))*1E4] ## 0.79 MgC/ha/yr v24; v23 0.98 MgC/ha/yr ground basis
biom.dat[aoi>800 & is.finite(isa.frac) & is.finite(can.frac), ((sum(npp.kg.hw.forest, na.rm=T)/2000)/sum(aoi))*1E4] ## 0.63 MgC/ha/yr v24; v23 0.63 MgC/ha/yr forest basis
biom.dat[aoi>800 & is.finite(isa.frac) & is.finite(can.frac), ((sum(npp.kg.hw.perv, na.rm=T)/2000)/sum(aoi))*1E4]  ## 0.60 MgC/ha/yr v24; 0.64 MgC/ha/yr perv basis

biom.dat[aoi>800, .(median(npp.kg.hw.ground, na.rm=T), ## median is ~60-90 kg-biomass/pix
                    median(npp.kg.hw.forest, na.rm=T),
                    median(npp.kg.hw.perv, na.rm=T))]

hist(biom.dat[aoi>800, ((npp.kg.hw.forest/2000)/aoi)*1E4]) ### sombitch, it too craps out at about 3 MgC/ha/yr
hist(biom.dat[aoi>800, ((npp.kg.hw.ground/2000)/aoi)*1E4]) ### 
hist(biom.dat[aoi>800, ((npp.kg.hw.perv/2000)/aoi)*1E4]) ### 
## all of these seem to recapitulate the productivity distributions seen IN THE EQUATION-BASED FIA THING WHATT???

write.csv(biom.dat, "processed/npp.FIA.empirV24.csv")

# biom.dat <- as.data.table(read.csv("processed/npp.FIA.empirV24.csv"))

#######
#### Supplemental analysis
#### following section contains:
#### 1) Analysis of weird artifact seen in stem-level growth~dbh plots (log-log) -- why parallel lines?
#### 2) Initial Boston NPP estimation based on published biomass-->Age-->C.aquisition equations provided by FIA COLE search
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


##### SUPPLEMENTAL 2
##### FIA, EQUATION BASIS BIOMASS-->AGE-->GROWTH
###### FIA V1: Equation for wood volume~time (proxy for stand biomass density)
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

##########
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

# ######################
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






