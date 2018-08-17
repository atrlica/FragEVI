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

spp.allo <- read.csv("data/FIA/spp_allometrics.csv")  ## lookup table of allometric equations for the most common representatives in the sample
live <- read.csv("data/FIA/MA_Tree_Data_ID_NOMORT_SUBID.csv")
live <- as.data.table(live)
names(live)[1] <- c("TreeID")
names(live)[2] <- c("PlotID")
names(live)[3] <- c("SubID")
spec <- read.csv("data/FIA/REF_SPECIES.csv") ## map of species # to spp
live <- merge(x=live, y=spec[,c("SPCD", "GENUS", "SPECIES")], by.x="SPECIES_CD", by.y="SPCD", all.x=T, all.y=F)
live$GENUS <- as.character(live$GENUS)
live$GENUS <- as.factor(live$GENUS)
live$GENUS.num <- as.numeric(live$GENUS)

### calculate species specific allometrics
live[,spp:=paste(substr(GENUS, 1,1), ".", SPECIES, sep="")]
live <- merge(x=live, y=spp.allo[,c("spp", "b0", "b1")], by="spp", all.x=T)
### if can't find a specific allometric, apply eastern hardwood default
live[is.na(b0), b0:=(-2.48)]
live[is.na(b1), b1:=2.4835]
biom.pred2 <- function(b0, b1, x){exp(b0+(b1*log(x)))}
live[,biom0.spp:=biom.pred2(b0, b1, DIAM_T0)]
live[,biom1.spp:=biom.pred2(b0, b1, DIAM_T1)]
## class as hard or soft wood
live[,type:="H"]
live[spp%in%c("P.strobus", "P.resinosa", "T.canadensis", "A.balsamea"), type:="S"]
live[,type:=as.factor(type)]
live[,biom.delt.spp:=biom1.spp-biom0.spp]
live[,growth.ann.rel:=(biom.delt.spp/biom0.spp)/4.8]
# summary(live$growth.ann.rel)

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


#### Additional grooming: Get rid of records from partially forested subplots
## first cull out the subplots that are too sparse (no fewer than 5 trees per subplot)
a <- live[,length(DIAM_T0), by=.(PlotID, SubID)]
a[,compID:=paste(PlotID, SubID, sep=".")]
kill.me <- a[V1<5, compID]
live[,compID:=paste(PlotID, SubID, sep=".")]
live <- live[!(compID%in%kill.me),] ## remove records of dbh if they come from subplots with fewer than 5 trees
write.csv(live, "processed/fia.live.stem.dbh.growth.csv")


##### Modeling growth~dbh (individual stem level)
## log model, growth>0, no hard/soft designation
live <- read.csv("processed/fia.live.stem.dbh.growth.csv")
live <- as.data.table(live) ## 6875 records
summary(live$growth.ann.rel) ## 0.9-3.4%, (-13-53%) -- so looks lower generally than Andy or Street trees
mod.fia.stem.log <- lm(log(growth.ann.rel)~log(DIAM_T0), data=live[growth.ann.rel>0])
m1 <- summary(mod.fia.stem.log) #R2 0.03, signficant
mod.fia.stem.type.log <- lm(log(growth.ann.rel)~log(DIAM_T0)*type, data=live[growth.ann.rel>0])
m2 <- summary(mod.fia.stem.type.log) # R2 0.04, type H/S not significant
plot(log(live$DIAM_T0), log(live$growth.ann.rel), col=as.numeric(live$type)+1) #S=2, H=1
abline(a=m2$coefficients[1], b=m2$coefficients[2], lty=2, col="gray55")
abline(a=m1$coefficients[1], b=m1$coefficients[2], lty=1, col="black")
## work in untransformed space (exponential model)
mod1.nls <- nls(growth.ann.rel ~ exp(a + b * log(DIAM_T0)), data=live, start=list(a=0, b=0)) ### OK THIS is the real exponential non-linear model that can handle the negatives
mm <- summary(mod1.nls) ## all factors significant
col.type <- c("green", "forestgreen")
plot(live$DIAM_T0, live$growth.ann.rel, col=col.type[as.numeric(live$type)], pch=15, cex=0.3,
     xlab="Stem diameter (cm)", ylab="Growth rate (kg/kg)", main="FIA stems")
test <- seq(live[,min(DIAM_T0)], live[,max(DIAM_T0)], length.out=100)
lines(test, exp(mm$coefficients[1]+(mm$coefficients[2]*log(test))), col="black", lwd=3)
legend(x=60, y=0.4, bty="n", legend = c("Hardwood", "Softwood"), fill=c("green", "forestgreen"))

## OK there you have it: a common model for growth~dbh, exponential, but low predictive power (lots of spread)

### what is mean growth rate at different intervals?
live[biom.delt.spp>0 & DIAM_T0<20, .(length(DIAM_T0), 
                                             median(growth.ann.rel, na.rm=T),
                                             mean(growth.ann.rel, na.rm=T),
                                             sd(growth.ann.rel, na.rm=T)*1.96)]
live[biom.delt.spp>0 & DIAM_T0>20 & DIAM_T0<30, .(length(DIAM_T0), 
                                             median(growth.ann.rel, na.rm=T),
                                             mean(growth.ann.rel, na.rm=T),
                                             sd(growth.ann.rel, na.rm=T)*1.96)]
live[biom.delt.spp>0 & DIAM_T0>30 & DIAM_T0<40, .(length(DIAM_T0), 
                                             median(growth.ann.rel, na.rm=T),
                                             mean(growth.ann.rel, na.rm=T),
                                             sd(growth.ann.rel, na.rm=T)*1.96)]
live[biom.delt.spp>0 & DIAM_T0>40 & DIAM_T0<50, .(length(DIAM_T0), 
                                             median(growth.ann.rel, na.rm=T),
                                             mean(growth.ann.rel, na.rm=T),
                                             sd(growth.ann.rel, na.rm=T)*1.96)]
live[biom.delt.spp>0 & DIAM_T0>50 , .(length(DIAM_T0), 
                                             median(growth.ann.rel, na.rm=T),
                                             mean(growth.ann.rel, na.rm=T),
                                             sd(growth.ann.rel, na.rm=T)*1.96)]

### so in general there is a lot of variability but a growth rate of about 2ish percent is fine for any diameter -- the model adds litle predictive power to this
ggplot(live[growth.ann.rel>0,], aes((DIAM_T0), (growth.ann.rel))) + geom_bin2d(bins=80) ## most values are low growth across the diameter range
ggplot(live[growth.ann.rel>0,], aes(log(DIAM_T0), log(growth.ann.rel))) + geom_bin2d(bins=80)

## Any obvious taxonomic split?
plot(log(live$DIAM_T0), log(live$growth.ann.rel), col=live$GENUS.num) ## no
table(live[growth.ann.rel>0, GENUS]) ## mostly acer pinus quercus tsuga betula
live[growth.ann.rel>0 & GENUS%in%c("Quercus", "Pinus", "Betula", "Tsuga", "Acer"), length(DIAM_T0)]/dim(live)[1] #82% in five genera
## this is my motivation for applying spp specific allometrics: Large # of Tsuga and Pinus predicted with default hardwood equation

## what is stem denstiy in each plot?
fack <- live[,.(length(DIAM_T0), sum(biom0.spp, na.rm=T)), by=PlotID]
hist(fack$V1) ## up to ~50 measured stems/plot, weak max below 10 -- these are culled for sparse subplots
plot(fack$V2, fack$V1) ## generally higher biomass == higher stem count

### OK: We know that most of the trees in the FIA plots grow more slowly than the stem records we have for street and Andy forests.

### why not make a nice plot? ## export 600x600
par(mfrow=c(1,1), mar=c(4,4,3,1))
main.xlim <- c(4, 100)
main.ylim <- c(-0.15, 1)
## FIA stem growth~dbh
col.type <- c("green", "forestgreen")
plot(live$DIAM_T0, live$growth.ann.rel, col=col.type[as.numeric(live$type)],
     pch=15, cex=0.4, xlim=main.xlim, ylim=main.ylim,
     xlab="Stem DBH (cm)", ylab="Relative growth (kg/kg)", main="FIA stem growth~DBH")
test <- seq(from=main.xlim[1], to=main.xlim[2], length.out=100)
lines(test, exp(mm$coefficients[1]+(mm$coefficients[2]*log(test))),
      col="black", lwd=3, lty=2)
legend(x=30, y=0.6, bty="n", legend = c("Hardwood", "Softwood"), fill=c("green", "forestgreen"))
abline(v=live[, median(DIAM_T0)])
abline(h=live[, median(growth.ann.rel)])



#####
#### Plot-level assessment of growth~biomass.density
### equivalent plot-level assessment in live-only data ## 331 plots have at least one subplot dense enough to bother with
live <- read.csv("processed/fia.live.stem.dbh.growth.csv")
live <- as.data.table(live)

live.plot <- live[,
                  .(length(unique(SubID)),
                    sum(biom.delt.spp, na.rm=T), 
                    sum(biom0.spp, na.rm=T),
                    length(DIAM_T0)),
                  by=PlotID]
dim(live.plot) ## 221 plots with at least one valid subplot
names(live.plot) <- c("PlotID", "num.subplots", "biom.growth", "total.biom0.kg", "num.stems")
plot(live.plot[,total.biom0.kg], live.plot[,num.stems]) ## linear but gets unconstrained above ~20 stems/plot
cor(live.plot[, total.biom0.kg], live.plot[, num.stems]) ## rho 0.60
## the "fully forested" plots sample 0 to 25k kg/plot 
live.plot[total.biom0.kg<5000,] ## mostly single-subplot plots (ie plot is mostly non-forest)


## plot-level growth
subplot.area <- (7.3152^2)*pi ## subplots are 24ft in radius
live.plot[,total.biom0.Mg.ha:=((total.biom0.kg/1000)/(num.subplots*subplot.area))*1E4] ## biomass density in aggregated area covered by the fully forested subplots
live.plot[,biom.growth.ann:=biom.growth/4.8]
live.plot[,biom.growth.ann.rel:=(biom.growth.ann)/total.biom0.kg]
hist(live.plot$biom.growth.ann.rel) ## most plots below 5%, a few are wildly productive, a couple decline
summary(live.plot$biom.growth.ann.rel) ## -1% to 14%

## plot level density
live.plot[,stems.ha:=(num.stems/(subplot.area*num.subplots))*1E4]
hist(live.plot$stems.ha) # 400-500 peak
hist(live.plot[, total.biom0.kg], breaks=10) ## peak 5-10k kg
hist((live.plot[, total.biom0.kg]/675)*1E4/2000, breaks=10) ## i.e. peaking 50-100 MgC/ha, up to 200
summary(live.plot[, total.biom0.kg/(subplot.area*num.subplots)*1E4]/2000) ## 90 MgC/ha, 34-307 -- so a forest

## growth~biomass modeling
plot(live.plot$total.biom0.kg, live.plot$biom.growth.ann.rel) ## exponential negative, the super productive ones were very low biomass to begin
plot(live.plot$total.biom0.kg, live.plot$biom.growth.ann) ## linear-ish, more biomass --> more growth (no bottoming out)
plot(live.plot$total.biom0.Mg.ha, live.plot$biom.growth.ann.rel) ## OK more happily exponential
plot(live.plot[,log(total.biom0.Mg.ha)], live.plot[,log(biom.growth.ann.rel)]) ## here's our regular log-log linear relationship, pretty shotgunned

l <- lm(log(biom.growth.ann.rel)~log(total.biom0.Mg.ha), data=live.plot) ### R2 0.17, significant
plot(live.plot[,log(total.biom0.Mg.ha)], live.plot[,log(biom.growth.ann.rel)])
s <- summary(l)
abline(a=s$coefficients[1], b=s$coefficients[2], col="red") ## ok ho hum
## can do this as a proper expnential nls to also incorporate the negative growth plots

# ## compare the models of areal biomass growth with raw data and using only hardwoods
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

### comparative models in hardwoods and softwoods in isolation
hw <- lm(live.plot[biom.growth.ann.hw>0,log(biom.growth.ann.hw)]~live.plot[biom.growth.ann.hw>0,log(total.biom0.Mg.ha)])
summary(hw) ## barely significant (p 0.07), R2 0.01
sw <- lm(live.plot[biom.growth.ann.sw>0,log(biom.growth.ann.sw)]~live.plot[biom.growth.ann.sw>0,log(total.biom0.Mg.ha)])
summary(sw) ## significant, R2 0.12

### model without respect to hard/softwood
tot.mod <- lm(live.plot[biom.growth.ann.rel>0,log(biom.growth.ann.rel)]~live.plot[biom.growth.ann.rel>0,log(total.biom0.Mg.ha)])
summary(tot.mod) ### ok does like R2 0.18 and significant

plot(live.plot[,log(total.biom0.Mg.ha)], live.plot[,log(biom.growth.ann.rel)], 
     main="FIA growth~biomass (plot level)", ylab="Relative growth rate (Mg/Mg/ha)", xlab="Plot biomass density (Mg/ha)",
     xaxt="n", yaxt="n")
axis(side = 1, at = c(4.5, 5, 5.5, 6, 6.5), labels = round(exp(c(4.5, 5, 5.5, 6, 6.5)), digits=0))
axis(side=2, at=seq(-5.5, -2.5, by=0.5), labels=round(exp(seq(-5.5, -2.5, by=0.5)), digits=3))
abline(a=tot.mod$coefficients[1], b=tot.mod$coefficients[2], lwd=1.4, col="black")
# points(live.plot[,log(((total.biom0.spp.kg.hw/1000)/(num.subplots*subplot.area))*1E4)],
#        live.plot[,log(biom.growth.ann.hw)], pch=15, cex=0.6, col="red")
points(live.plot[,log(total.biom0.Mg.ha)],
       live.plot[,log(biom.growth.ann.hw)], pch=15, cex=0.6, col="red")
abline(a=hw$coefficients[1], b=hw$coefficients[2], col="red")
# points(live.plot[,log(((total.biom0.spp.kg.sw/1000)/(num.subplots*subplot.area))*1E4)],
#        live.plot[,log(biom.growth.ann.sw)], pch=17, cex=0.6, col="blue")
points(live.plot[,log(total.biom0.Mg.ha)],
       live.plot[,log(biom.growth.ann.sw)], pch=17, cex=0.6, col="blue")
abline(a=sw$coefficients[1], b=sw$coefficients[2], col="blue")
legend(x=4.3, y=-5, legend = c("All", "Hardwoods", "Softwoods"), 
       fill = c("black", "red", "blue"))
summary(live.plot$biom.growth.ann.hw) ## hardwood growth rate just doesn't vary much across range of TOTAL biomass density: 2-3% growth (0.4 to 15%)


#### let's look at it in untransformed space
tot.mod <- lm(live.plot[biom.growth.ann.rel>0,log(biom.growth.ann.rel)]~live.plot[biom.growth.ann.rel>0,log(total.biom0.Mg.ha)])
summary(tot.mod) ### R2 0.18 and significant
hw.mod <- lm(live.plot[biom.growth.ann.hw>0,log(biom.growth.ann.hw)]~live.plot[biom.growth.ann.hw>0,log(total.biom0.Mg.ha)])
summary(hw.mod) ### R2 0.01 and only sig at p<0.10
hw.mod$coefficients
tot.mod$coefficients ## slope in hw.mod is very low, not particularly high in tot.mod

## untransformed
hw.mod.exp <- nls(biom.growth.ann.hw ~ exp(a + b * log(total.biom0.Mg.ha)),
                  data=live.plot, start=list(a=0, b=0))
q <- summary(hw.mod.exp) ## b is not even close to significant
exp(q$coefficients[1]) ## ie. about 3% boyo

tot.mod.exp <- nls(biom.growth.ann.rel ~ exp(a + b * log(total.biom0.Mg.ha)),
                   data=live.plot, start=list(a=0, b=0))
r <- summary(tot.mod.exp) ## definitely significant

## visual comparison of the plot-level models
plot(live.plot[,(total.biom0.Mg.ha)], 
     live.plot[,(biom.growth.ann.rel)], main="FIA growth~biomass (plot-level)",
     xlab="Biomass Density (Mg/ha)", ylab="Relative growth (Mg/Mg/ha)")
points(live.plot[,(total.biom0.Mg.ha)], 
       live.plot[,(biom.growth.ann.hw)], col="red", pch=17, cex=0.6)
test <- seq(min(live.plot$total.biom0.Mg.ha), max(live.plot$total.biom0.Mg.ha), length.out=100)
lines(test, exp(r$coefficients[1]+(r$coefficients[2]*log(test))),
      lwd=3, lty=2)
# abline(h=mean(live.plot$biom.growth.ann.hw), col="orangered", lwd=3, lty=4)

## why so poor fit on hw model, why some high biomass plots have such high hw growth rate?
# summary(live.plot$hw.frac) ## most have a majority hardwoods, a few v low
# hist(live.plot$hw.frac)
# live.plot[total.biom0.Mg.ha>400, hw.frac] ## high biomass contains the low hw fraction ones
points(live.plot[hw.frac<0.25, total.biom0.Mg.ha],
       live.plot[hw.frac<0.25, biom.growth.ann.hw],
       cex=1.2, lwd=3, pch=5, col="seagreen")

### OOOOOOKKKKKKAAAAYYYYY Let's elminate the handful of weird low-HW plots (could be weird places that favor pines or something, not really how an urban forest do)
hw.mod.exp.filt <- nls(biom.growth.ann.hw ~ exp(a + b * log(total.biom0.Mg.ha)),
                       data=live.plot[hw.frac>0.25,], start=list(a=0, b=0))
y <- summary(hw.mod.exp.filt) ## b is not even close to significant
y$coefficients
r$coefficients ## quite comparable
## so we will use the low-hw filtered plots to approximate growth in the final NPP calculations
### show on plot
lines(test, exp(y$coefficients[1]+(y$coefficients[2]*log(test))),
      lwd=3, lty=2, col="red")
legend(x=100, y=0.01, legend=c("Hardwoods", "Low HW frac.", "All"), 
       fill=c("red",  "seagreen", "black"), bty="n")




########
##### PART B: BOSTON NPP FROM FIA EMPRICAL GROWTH~BIOMASS ANALYSIS
######
### V2.2: 1) uses species-specific biometrics; 2) models hardwood growth separately from trees in general; 3) uses nls to avoid dumping negative growth records
### V2.3 1) Uses subplot IDs to remove dbh records from subplot sites that are not fully forested; 2) filtered plots that have low hw fraction to determine hw-only growth rate

## read in data and get plot-level metrics
live <- read.csv("processed/fia.live.stem.dbh.growth.csv")
live <- as.data.table(live)

live.plot <- live[,
                  .(length(unique(SubID)),
                    sum(biom.delt.spp, na.rm=T), 
                    sum(biom0.spp, na.rm=T),
                    length(DIAM_T0)),
                  by=PlotID]
names(live.plot) <- c("PlotID", "num.subplots", "biom.growth", "total.biom0.kg", "num.stems")
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
subplot.area <- (7.3152^2)*pi ## subplots are 24ft in radius
live.plot[,total.biom0.Mg.ha:=((total.biom0.kg/1000)/(num.subplots*subplot.area))*1E4] ## biomass density in aggregated area covered by the fully forested subplots
live.plot[,biom.growth.ann:=biom.growth/4.8]
live.plot[,biom.growth.ann.rel:=(biom.growth.ann)/total.biom0.kg]
live.plot[,hw.frac:=total.biom0.spp.kg.hw/total.biom0.kg]

#######
## final exponential model fit, hardwood growth~biomass
hw.mod.exp.filt <- nls(biom.growth.ann.hw ~ exp(a + b * log(total.biom0.Mg.ha)),
                       data=live.plot[hw.frac>0.25,], start=list(a=0, b=0))
y <- summary(hw.mod.exp.filt) ## residual standard error = standard error of regression = how far off values may be from predicted (vs. R2, which is unreliable)

## load the biomass data and reprocess
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

### Now we have to decide how to calculate the "forest" density of the biomass on the ground
### Three approaches: GROUND (kg-biomass/m2-pixel); "FOREST" (kg-biomass/m2-canopy); "PERV" (kg-biomass/m2-pervious)
### live MgC by area of GROUND in each pixel
## Ground biomass density
biom.dat[, live.Mgbiom.ha.ground:=(bos.biom30m/aoi)*(1E4)/1E3] ## convert kg-biomass/pixel based on size of pixel (summed by 1m2 in "aoi"), 
## Forest biomass density
biom.dat[,live.Mgbiom.ha.forest:=(bos.biom30m/(aoi*can.frac))*(1E4)/1E3]
biom.dat[bos.biom30m==0, live.Mgbiom.ha.forest:=0] ## have to manually fix this because of 0 canopy pix
biom.dat[can.frac<0.01, live.Mgbiom.ha.forest:=0]
## Perv biomass density
biom.dat[,live.Mgbiom.ha.perv:=(bos.biom30m/(aoi*(1-isa.frac)))*(1E4)/1E3]
biom.dat[bos.biom30m==0, live.Mgbiom.ha.perv:=0] ## have to manually fix this because of isa=1 pix
biom.dat[isa.frac>0.99, live.Mgbiom.ha.perv:=0]

### what do our densities look like
hist(biom.dat[aoi>800, live.Mgbiom.ha.ground]) ## up to 600 Mgbiom/ha = 300 MgC/ha -- a lot of this below the range sampled by the FIA plot data
summary(biom.dat[aoi>800, live.Mgbiom.ha.ground])
hist(biom.dat[aoi>800, live.Mgbiom.ha.forest]) ## same, more medium sized
summary(biom.dat[aoi>800,live.Mgbiom.ha.forest]) ## everything nudged denser
hist(biom.dat[aoi>800 & isa.frac<0.98, live.Mgbiom.ha.perv]) ## up to 800k kgbiom/ha
summary(biom.dat[aoi>800, live.Mgbiom.ha.perv])
### contrast: what is the range sampled in the live.plot FIA data?
summary(live.plot[,total.biom0.Mg.ha]) ## don't get into the low range of things! Minimum biomass density is about where most of the pixels fall

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
biom.dat[,ground.gfact:=exp(y$coefficients[1]+y$coefficients[2]*log(live.Mgbiom.ha.ground))]
biom.dat[,forest.gfact:=exp(y$coefficients[1]+y$coefficients[2]*log(live.Mgbiom.ha.forest))]
biom.dat[,perv.gfact:=exp(y$coefficients[1]+y$coefficients[2]*log(live.Mgbiom.ha.perv))]

## kill some artifacts
biom.dat[bos.biom30m<1 | live.Mgbiom.ha.ground<1, ground.gfact:=0]
biom.dat[bos.biom30m<1 | live.Mgbiom.ha.forest<1, forest.gfact:=0]
biom.dat[bos.biom30m<1 | live.Mgbiom.ha.perv<1, perv.gfact:=0]

summary(biom.dat[aoi>800, ground.gfact])
summary(biom.dat[aoi>800, forest.gfact])
summary(biom.dat[aoi>800, perv.gfact])

# View(biom.dat[is.finite(bos.biom30m) & is.na(forest.gfact),]) ## good news -- these don't exist, just 111 weird out-of-aoi pixels that have biomass readings

biom.dat[live.Mgbiom.ha.ground>60 & live.Mgbiom.ha.ground<150, .(mean(ground.gfact), length(ground.gfact))] ## about 3% growth in pixels that fall near where the FIA data fell
biom.dat[live.Mgbiom.ha.forest>60 & live.Mgbiom.ha.forest<150, .(mean(forest.gfact), length(forest.gfact))]
biom.dat[live.Mgbiom.ha.forest>60 & live.Mgbiom.ha.forest<150, .(mean(forest.gfact), length(forest.gfact))]
biom.dat[bos.biom30m>0, length(bos.biom30m)] ## 108k pix have biomass
summary(live.plot[,total.biom0.Mg.ha]) ## contrast, FIA plots don't sample low biomass density below 70 Mg/ha

#### OK: This approach has us extrapolating from our field data to estimate (high) growth rates in the low-biomass parts of the city
#### In contrast, the previous equation-based approach saw most of our real forest falling off the edge of the age/density curve and predicting 0 growth, while partial forest looked "young"

## calculate npp from these growth factors
## regression coeff*biom.density(kg/ha, 3 approaches)-->growth factors (kg/kg) per cell
## growth factors * cell biomass (kg) --> npp (kg biomass per cell) 
biom.dat[,npp.kg.hw.ground:=bos.biom30m*ground.gfact]
biom.dat[,npp.kg.hw.forest:=bos.biom30m*forest.gfact]
biom.dat[,npp.kg.hw.perv:=bos.biom30m*perv.gfact]

### totals for aoi
biom.dat[aoi>800 & is.finite(isa.frac) & is.finite(can.frac), sum(npp.kg.hw.ground, na.rm=T)]/2000 ## 12k MgC/yr ground basis
biom.dat[aoi>800 & is.finite(isa.frac) & is.finite(can.frac), sum(npp.kg.hw.forest, na.rm=T)]/2000 ## 7.8k MgC/yr forest basis
biom.dat[aoi>800 & is.finite(isa.frac) & is.finite(can.frac), sum(npp.kg.hw.perv, na.rm=T)]/2000  ## 7.8k MgC/yr perv basis
biom.dat[aoi>800 & is.finite(isa.frac) & is.finite(can.frac), ((sum(npp.kg.hw.ground, na.rm=T)/2000)/sum(aoi))*1E4] ## 0.98 MgC/ha/yr ground basis
biom.dat[aoi>800 & is.finite(isa.frac) & is.finite(can.frac), ((sum(npp.kg.hw.forest, na.rm=T)/2000)/sum(aoi))*1E4] ## 0.63 MgC/ha/yr forest basis
biom.dat[aoi>800 & is.finite(isa.frac) & is.finite(can.frac), ((sum(npp.kg.hw.perv, na.rm=T)/2000)/sum(aoi))*1E4]  ## 0.64 MgC/ha/yr perv basis

biom.dat[aoi>800, .(median(npp.kg.hw.ground, na.rm=T), ## median is ~50-60 kg-biomass/pix
                    median(npp.kg.hw.forest, na.rm=T),
                    median(npp.kg.hw.perv, na.rm=T))]

hist(biom.dat[aoi>800, ((npp.kg.hw.forest/2000)/aoi)*1E4]) ### sombitch, it too craps out at about 3 MgC/ha/yr
hist(biom.dat[aoi>800, ((npp.kg.hw.ground/2000)/aoi)*1E4]) ### 
hist(biom.dat[aoi>800, ((npp.kg.hw.perv/2000)/aoi)*1E4]) ### 
## all of these seem to recapitulate the productivity distributions seen IN THE EQUATION-BASED FIA THING WHATT???

write.csv(biom.dat, "processed/npp.FIA.empirV23.csv")


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
### FIA V1: Equation for wood volume~time (proxy for stand biomass density)
a=123.67
b=0.04
AGE=seq(0,500)
y = a*(1-exp(-b*AGE))^3
plot(AGE, y, ylab="stand MgC/ha")
z = diff(y) ## yearly growth increment, absolute kg
z.rel <- z/y[2:501] ## yearly growth increment, % of previous year biomass
plot(y[2:501], z,  xlab="biomass MgC/ha", ylab="absolute growth rate, MgC/ha/yr") ## zero growth above ~125 MgC/ha
points(y[2:501], z.rel, pch=14, cex=0.5, col="red") #hyperbolic
plot(AGE[2:501], z, xlab="Stand age, yr", ylab="absolute growth rate, MgC/ha/yr") ## this is the gain curve over site age, near zero >200 yrs
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






