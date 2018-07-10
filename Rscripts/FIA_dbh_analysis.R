library(raster)
library(data.table)
library(ggplot2)

### process and assessment of FIA individual dbh records provided by Moreale
### individual FIA tree data from sites near Boston
# dat <- read.csv("data/FIA/MA_Tree_Data.csv")
dat <- read.csv("data/FIA/MA_Tree_Data_ID.csv")
dat <- as.data.table(dat)
names(dat)[1] <- c("TreeID")
names(dat)[2] <- c("PlotID")
spec <- read.csv("data/FIA/REF_SPECIES.csv")
dat <- merge(x=dat, y=spec[,c("SPCD", "GENUS")], by.x="SPECIES_CD", by.y="SPCD", all.x=T, all.y=F)
dat$GENUS <- as.character(dat$GENUS)
dat$GENUS <- as.factor(dat$GENUS)
table(dat$GENUS)
dat$GENUS.num <- as.numeric(dat$GENUS)

### using eastern hardwood defaults, but a lot of the record are Pinus or Tsuga. Might be good to go ahead and use the correct generic equation
b0 <- -2.48
b1 <- 2.4835 ## these are eastern hardwood defaults
biom.pred <- function(x){exp(b0+(b1*log(x)))}
dat$biom0 <- biom.pred(dat$DIAM_T0)
dat$biom1 <- biom.pred(dat$DIAM_T1)
dat$biom.delt <- dat$biom1-dat$biom0
dat$biom.rel <- (dat$biom.delt/4.8)/dat$biom0 ## annualized relative growth increment
dat$dbh.delt <- dat$DIAM_T1-dat$DIAM_T0

#### what is the deal with these scatter plots, artifact bug hunt
plot(dat$DIAM_T0, log(dat$biom.delt))
plot(log(dat$DIAM_T0), log(dat$biom.delt)) ## linear but only in the ones that are not falling (non negative)
plot(log(dat$DIAM_T0), log(dat$biom.delt),
     col=dat$GENUS.num) ## a lot of one particular genus
plot(log(dat$DIAM_T0), log(dat$biom.delt),
     col=dat$SPECIES_CD) ## no clear plot ID in lower parts
plot(dat$DIAM_T0, dat$biom.delt, col=dat$GENUS) ## some are diving
dat[biom.delt<(-500),] ## only 3 lose more than 500kg
View(dat[biom.delt<exp(1) & biom.delt>0,]) ## what are the tiny growth records
plot(dat[log(biom.delt)<1.45 & log(DIAM_T0)<3, log(DIAM_T0)], dat[log(biom.delt)<1.45 & log(DIAM_T0)<3, log(biom.delt)])
View(dat[log(biom.delt)<1.45 & log(DIAM_T0)<3, ]) ## a straight line! wtf???
dat[,diam.delt:=DIAM_T1-DIAM_T0]

### here's what's happening: The low straight-line changes are where they are detecting a 0.1" increase in circum. (minimum increase on dbh tape) -- locks things into a single growth line
table(dat[log(biom.delt)<1.5 & log(DIAM_T0)<3, GENUS]) ## dominated by a handfull of genera
table(dat[log(biom.delt)<1.5 & log(DIAM_T0)<3, PlotID]) ## these low growth small ones are everywhere
hist(dat[log(biom.delt)<1.5 & log(DIAM_T0)<3, DIAM_T0]) ## fairly even dbh dist between 12 and 20cm
plts <- unique(dat$PlotID)
for(i in 1:length(plts)){ ### in general, a positive trend of greater delta.biom with greater dbh to start
  plot(dat[PlotID==plts[i], log(DIAM_T0)], 
       dat[PlotID==plts[i],log(biom.delt)], 
       col=dat[PlotID==plts[i],GENUS.num], main=paste("plot", plts[i]))
}
### hang on: a lot of these are on the same slope
hist(dat[dbh.delt>0, dbh.delt], breaks=50)
dat[dbh.delt>2,]
dat[,dbh.delt.incr:=dbh.delt/0.254] ## 0.245 is the dbh tape increment
dat[dbh.delt.incr>0 & dbh.delt.incr<3,]
plot(dat[,log(DIAM_T0)], dat[,log(biom.rel)], col=as.factor(dat[,dbh.delt.incr]))
View(dat[log(DIAM_T0)<2.8&log(biom.rel)<(-4.1),]) ## a huge proportion of these are 0 growth
dat[biom.rel<=0,] ### 2k records show 0 to negative growth
plot(dat[log(DIAM_T0)<2.8&log(biom.rel)<(-4.1), log(DIAM_T0)], 
     dat[log(DIAM_T0)<2.8&log(biom.rel)<(-4.1), log(biom.rel)],
     col=as.factor(dat[,dbh.delt.incr]))
r <- dat[log(DIAM_T0)<2.55&log(biom.rel)<(-4.1),]
View(r[order(r$dbh.delt.incr),])
dat[, dbh.units0:=DIAM_T0/0.254] ## convert dbh to tape increment units
dat[, dbh.units1:=DIAM_T1/0.254]
plot(dat$dbh.units0, dat$dbh.units1) ### these are def integers
plot(dat[dbh.units1-dbh.units0>0, dbh.units0], dat[dbh.units1-dbh.units0>0, dbh.units1])
abline(a=0,b=1) ## positive growth
plot(dat[dbh.units1-dbh.units0<0, dbh.units0], dat[dbh.units1-dbh.units0<0, dbh.units1])
abline(a=0,b=1) ## negative growth
plot(dat[dbh.units1-dbh.units0==0, dbh.units0], dat[dbh.units1-dbh.units0==0, dbh.units1])
abline(a=0,b=1) ## zero growth

dat[,dbh.unit.delt:=dbh.units1-dbh.units0]
dat[log(DIAM_T0)<2.55&log(biom.rel)<(-4.1), unique(dbh.unit.delt)] ## for the low lines it's all 0 or 1 dbh increment
t.col <- c("blue", "purple", "red")
plot(dat[log(DIAM_T0)<2.8&log(biom.rel)<(-4.1), log(DIAM_T0)], ### everything is 1 dbh increment growth
     dat[log(DIAM_T0)<2.8&log(biom.rel)<(-4.1), log(biom.rel)],
     col=t.col[dat[log(DIAM_T0)<2.8&log(biom.rel)<(-4.1), dbh.unit.delt]])

View(dat[log(DIAM_T0)<2.8&log(biom.rel)<(-4.1) & log(biom.rel)>(-4.2), ])
View(dat[log(DIAM_T0)<2.8&log(biom.rel)<(-4.2), ])
plot(dat[DIAM_T0<exp(2.8) & biom.rel<exp(-4.2), DIAM_T0],
     dat[DIAM_T0<exp(2.8) & biom.rel<exp(-4.2), biom.rel],
     col=t.col[abs(dat[DIAM_T0<exp(2.8) & biom.rel<exp(-4.2), dbh.unit.delt])])

t.col <- c("blue", "red", "black", "green", "orange", "yellow", rep("salmon",100))
plot(dat[, DIAM_T0],
     dat[, biom.rel],
     col=t.col[abs(dat[, dbh.unit.delt])])
plot(dat[biom.rel>0, DIAM_T0],
     dat[biom.rel>0, biom.rel],
     col=t.col[abs(dat[biom.rel>0, dbh.unit.delt])])
plot(dat[dbh.unit.delt==1, DIAM_T0],
     dat[dbh.unit.delt==1, biom.rel], ylim=c(-0.05, 0.05), col="green")
points(dat[dbh.unit.delt==0, DIAM_T0],
     dat[dbh.unit.delt==0, biom.rel], col="black") ## look like they're on top of each other
points(dat[dbh.unit.delt==2, DIAM_T0],
       dat[dbh.unit.delt==2, biom.rel], col="blue") ## look like they're on top of each other
points(dat[dbh.unit.delt==3, DIAM_T0],
       dat[dbh.unit.delt==3, biom.rel], col="purple") ## look like they're on top of each other
points(dat[dbh.unit.delt==-1, DIAM_T0],
       dat[dbh.unit.delt==-1, biom.rel], col="orange") ## look like they're on top of each other
points(dat[dbh.unit.delt==-2, DIAM_T0],
       dat[dbh.unit.delt==-2, biom.rel], col="gold") ## look like they're on top of each other
points(dat[dbh.unit.delt==-3, DIAM_T0],
       dat[dbh.unit.delt==-3, biom.rel], col="pink") ## look like they're on top of each other
## congratulations, you've proved that the dbh~biomass allometric is exponential
## basically every linear line you get is a new click on the line: Higher => jumped more dbh units
## but in the log transform you are cutting off all the values that are negative dbh units

plot(dat[dbh.unit.delt%in%seq(-4,4), DIAM_T0],
     dat[dbh.unit.delt%in%seq(-4,4), biom.rel], col=dat[dbh.unit.delt%in%seq(-4,4), dbh.unit.delt+5])

plot(dat[dbh.unit.delt%in%seq(-10,10), DIAM_T0],
     dat[dbh.unit.delt%in%seq(-10,10), biom.rel], col=dat[dbh.unit.delt%in%seq(-10,10), dbh.unit.delt+11])

dat[dbh.unit.delt%in%seq(-10,10), ] ## 8018 records of 9402 (85%)
### log transforming it gives you this one (byproduct is it elminates the negative growth records entirely)
plot(dat[dbh.unit.delt%in%seq(1,10), log(DIAM_T0)],
     dat[dbh.unit.delt%in%seq(1,10), log(biom.rel)], col=dat[dbh.unit.delt%in%seq(1,10), dbh.unit.delt])
plot(log(dat$DIAM_T0), log(dat$biom.rel)) ## yep, each separate line is a different dbh increment class
plot(dat[dbh.unit.delt%in%seq(1,50), log(biom0)], ## same feckin story
     dat[dbh.unit.delt%in%seq(1,50), log(biom.rel)], col=dat[dbh.unit.delt%in%seq(1,50), dbh.unit.delt])
## so counting up from the bottom, each line represents 1, 2, 3, etc dbh increment gain, across a range of starting dbh (and bigger dbh class gives *relatively* less biomass gain at a given dbh increment class, but absolutely more)

plot(dat$DIAM_T0, dat$dbh.delt.incr) ## how much absolute dbh growth do you get?
plot(log(dat$DIAM_T0), log(dat$biom.delt)) ## feck same thing


### now to describe the growth curves for these mfs
### Here's the logic: 
## 1) The 1m biomass map shows how much is alive at time X, and we are using it to predict NPP up to time Y (1 year)
## 2) The Andy forest calc relies on predicting forest biomass gain *by area* using growth rate~dbh *per stem* for a selection of individual trees
## 2b) The growth rates per stem were determined *from trees that lived* -- no dbh in the forest area was measured for dead wood, no growth mortality rate was estimated for the stand
## 2c) So the Andy growth rate (kg/ha/yr) is estimated as if all biomass survives through the growth period
## 3) Same logic applies to the street tree records: Growth rate *per stem* is only measurable in trees that survived
## 3b) Ian has some model projections on mortality rates, but on an areal basis we only can start with the biomass we have (had, 2006/7) and project forward as if it all survives
## 4) The FIA equation approach is meant apply to a *forest area* directly, and may presumably incorporate mortalities and recruitment along the way
## 5) Alternatively, we can recalculate the areal biomass growth rate as a sum of stem growth~dbh, and determine the growth~dbh relationship from individual stem records
## 5b) To be comparable to results based on the other data sets, the growth~dbh relationship should be estimated based on *the trees that lived* and not on a generalized resample of all trees (as the FIA measures everything, living or dead)

hist(dat$DIAM_T0) ## ok, nothing gigantic, there's no old-growth monsters out there
dat[,median(DIAM_T0, na.rm=T)] ##22cm, perfectly comparable to the street trees or to the Andy trees
summary(dat$DIAM_T0) ## 9 to 91, pretty tight central range
summary(dat$biom.rel) ## median is 1.76% growth, range -16% to 52%, tight central range of 0.5 to 3.2% growth

mod.fia.live <- lm(log(biom.rel)~log(DIAM_T0), data=dat[biom.delt>0]) ## basic log-log growth~dbh
summary(mod.fia.live) ## diam0 is significant, but miserable predictive power R2 0.04

### what is mean growth rate at different intervals?
dat[biom.delt>0 & DIAM_T0<20, .(length(DIAM_T0), ## 2700, 2.76%, 4.1% +-8.7%
                                             median(biom.rel, na.rm=T),
                                             mean(biom.rel, na.rm=T),
                                             sd(biom.rel, na.rm=T)*1.96)]
dat[biom.delt>0 & DIAM_T0>20 & DIAM_T0<30, .(length(DIAM_T0), ## 2400, 2.22%, 2.8% +-4.8%
                                             median(biom.rel, na.rm=T),
                                             mean(biom.rel, na.rm=T),
                                             sd(biom.rel, na.rm=T)*1.96)]
dat[biom.delt>0 & DIAM_T0>30 & DIAM_T0<40, .(length(DIAM_T0), ## 1300, 2.2%, 2.4% +-3.2%
                                             median(biom.rel, na.rm=T),
                                             mean(biom.rel, na.rm=T),
                                             sd(biom.rel, na.rm=T)*1.96)]
dat[biom.delt>0 & DIAM_T0>40 & DIAM_T0<50, .(length(DIAM_T0), ## 500, 2.18%, 2.4% +-2.6%
                                             median(biom.rel, na.rm=T),
                                             mean(biom.rel, na.rm=T),
                                             sd(biom.rel, na.rm=T)*1.96)]
dat[biom.delt>0 & DIAM_T0>50 , .(length(DIAM_T0), ## 300, 1.9%, 2.1% +-2.4%
                                             median(biom.rel, na.rm=T),
                                             mean(biom.rel, na.rm=T),
                                             sd(biom.rel, na.rm=T)*1.96)]

### so in general there is a lot of variability but a growth rate of about 2ish percent is fine for any diameter -- the model adds litle predictive power to this
plot(dat[biom.delt>0, DIAM_T0], dat[biom.delt>0, biom.rel]) ## just some of the high boys in the low end of dbh
ggplot(dat[biom.delt>0,], aes((DIAM_T0), (biom.rel))) + geom_bin2d(bins=80) ## most values are low growth across the diameter range
ggplot(dat[biom.delt>0,], aes(log(DIAM_T0), log(biom.rel))) + geom_bin2d(bins=80)

## Any obvious taxonomic split?
plot(log(dat$DIAM_T0), log(dat$biom.rel), col=dat$GENUS.num) ## no
table(dat[biom.delt>0, GENUS]) ## mostly acer pinus quercus tsuga betula
dat[biom.delt>0 & GENUS%in%c("Quercus", "Pinus", "Betula", "Tsuga", "Acer"),] #90%
## disturbs me we are doing biomass on Tsuga and Pinus using default hardwood equation
## what is stem denstiy in each plot?
fack <- dat[,length(DIAM_T0), by=PlotID]
hist(fack$V1) ## up to 60 measured stems/plot, weak max at 20-30, surely the higher stem counts are the heavier plots

###
plot(log(dat$DIAM_T0), log(dat$biom.rel), col=dat$GENUS) ## OK, this is the plot to compare with the Andy and Street records
d <- summary(lm(log(biom.rel)~log(DIAM_T0), data=dat[biom.delt>0,])) ### miserable, R2 0.04
b0 <- d$coefficients[1]
b1 <- d$coefficients[2]
## sanity check: What do we predict for growth.rel across dbh?
par(mar=c(4,4,3,1))
plot(dat$DIAM_T0, dat$biom.rel, 
     col="gray45", cex=0.4, 
     xlab="DBH, cm", ylab="Biomass growth rate (kg/kg)", main="FIA growth rates")
abline(h=0, col="blue")
points(dat$DIAM_T0, exp(b0)*exp(b1*log(dat$DIAM_T0)), col="red", cex=0.3, pch=16) ## yeah, weak ass exponential 
median(exp(b0)*exp(b1*log(dat[DIAM_T0<20, DIAM_T0]))) ## 2.7% for <20cm
median(exp(b0)*exp(b1*log(dat[DIAM_T0>40 & DIAM_T0<60, DIAM_T0]))) ## 1.7% for 40-60cm
median(exp(b0)*exp(b1*log(dat[DIAM_T0>60, DIAM_T0]))) ## 1.4% for >60cm  ## OK these are all basically reasonable

###untransformed
plot(dat[,DIAM_T0], dat[,biom.rel], col=dat[,GENUS])
hist(dat[DIAM_T0>60, biom.rel]) ## at high biomass it stabilizes somewhat below the grand mean for growth
hist(dat[DIAM_T0<20, biom.rel]) ## get some higher rates at low dbh
summary(dat[DIAM_T0<20, biom.rel]) ## mean growth in small trees is close to grand mean, is a nls exponential model really more informative than the mean?

### OK: We know that trees in the FIA plots generally grow more slowly than the stem records we have for street and Andy forests.
### Let us look at measured growth (in biomass that survives from t0 to t1) as a %growth of living biomass in each plot 

## How much %biomass did each plot gain per year as a function of standing biomass in T0, and apply that 
options(scipen=000)

## stem records that might indicate mortality
dat[biom.delt==0,] ## 1500 records record no biomass gain -- any way to tell if these are dead or asleep?
dat[biom.delt<0,] ## another 500 records show biomass loss
dim(dat[biom.delt<=0,])[1]/dim(dat)[1] ## 21% of stems do not see biomass gain

### excluding non-growing stems, growth in each plot
biom.plot <- dat[biom.delt>0,
                 .(sum(biom.delt, na.rm=T), sum(biom0, na.rm=T)), by=PlotID]
names(biom.plot)[2:3] <- c("biom.growth", "total.biom0.kg")
biom.plot[,biom.growth.ann:=biom.growth/4.8]
biom.plot[,biom.growth.ann.rel:=biom.growth.ann/total.biom0.kg]
hist(biom.plot$total.biom0.kg) ## peak 5-10k kg
hist(biom.plot$biom.growth.ann.rel) ## most plots below 5%, a few are wildly productive, 25%
biom.plot[,median(biom.growth.ann.rel, na.rm=T)] ## median 2.5% gain per plot (in line with basic productivity model)
(biom.plot[,median(total.biom0.kg, na.rm=T)/2000]/675)*1E4 ## median 51 MgC/ha/yr
plot(biom.plot$total.biom0.kg, biom.plot$biom.growth.ann.rel) ## exponential negative, the super productive ones were very low biomass to begin
plot(biom.plot$total.biom0.kg, biom.plot$biom.growth.ann) ## linear-ish, more biomass --> more growth (no bottoming out)
biom.plot[,biom0.kg.ha:=(total.biom0.kg/675)*1E4]
hist(biom.plot$biom0.kg.ha/1000) ## some of these plots are very high biomass, well above the equation cut off for positive growth
plot(biom.plot$biom0.kg.ha/1000, biom.plot$biom.growth.ann) 
plot(biom.plot$biom0.kg.ha/1000, biom.plot$biom.growth.ann.rel) ##  contrast the equation version that says anything >120 Mg/ha doesn't gain live biomass
abline(h=0, col="red")

summary(biom.plot$biom.growth.ann.rel) ## 0.4-25%, fairly tight 2-3%
summary(biom.plot$biom0.kg.ha/1000) ## 0.7 <- 384, 66-154 middle part -- some are nearly not forest!
hist(log(biom.plot$biom.growth.ann.rel)) ## this is more normal
hist(log(biom.plot$total.biom0.kg)) ## this is more normal but not that normal, skewed left (a handful of super low biomass plots)
plot(log(biom.plot$total.biom0.kg), log(biom.plot$biom.growth.ann.rel))
m <- lm(log(biom.growth.ann.rel)~log(total.biom0.kg), data=biom.plot) ### not great
lines(log(biom.plot$total.biom0.kg), predict(m), col="red", lty=2)
summary(m) ## really not amazing, R2=0.2, but fuck it, not a big range, can apply this model to basically any pixel, covers that range fine
# this can be our model of biomass (kg/m2) to npp (kg/m2/yr)


## fun with non-transformed non-linear exponential model fitting!!1
### exponential model of growth~biomass (simple y=exp(a+bx))
mod <- nls(biom.growth.ann.rel ~ exp(a + b * total.biom0.kg), data=biom.plot, start=list(a=0, b=0)) ### OK THIS is the real exponential non-linear model that can handle the negatives
summary(mod)
plot(biom.plot$total.biom0.kg, biom.plot$biom.growth.ann.rel, col="gray55", pch=14, cex=0.3)
points(biom.plot$total.biom0.kg, predict(mod), col="red", pch=16, cex=0.3)

### the transformed log(y)~a+bx model -- R2 only 0.08
mod.lm <- lm(biom.plot[biom.growth.ann.rel>0, log(biom.growth.ann.rel)]~biom.plot[biom.growth.ann.rel>0, total.biom0.kg])
summary(mod.lm)
plot(biom.plot[biom.growth.ann.rel>0, total.biom0.kg], biom.plot[biom.growth.ann.rel>0, log(biom.growth.ann.rel)])
points(biom.plot[biom.growth.ann.rel>0, total.biom0.kg], predict(mod.lm), col="red", pch=14, cex=0.3)

### a log-log model
mod.loglog <- lm(biom.plot[biom.growth.ann.rel>0, log(biom.growth.ann.rel)]~biom.plot[biom.growth.ann.rel>0, log(total.biom0.kg)])
summary(mod.loglog) ## R2 0.14
plot(biom.plot[biom.growth.ann.rel>0, log(total.biom0.kg)], biom.plot[biom.growth.ann.rel>0, log(biom.growth.ann.rel)])
points(biom.plot[biom.growth.ann.rel>0, log(total.biom0.kg)], predict(mod.loglog), col="red", pch=14, cex=0.3)

## log log in nonlinear approach
mod2 <- nls(biom.growth.ann.rel ~ exp(a + b * log(total.biom0.kg)), data=biom.plot, start=list(a=0, b=0)) ### OK THIS is the real exponential non-linear model that can handle the negatives
summary(mod2)
plot(biom.plot$total.biom0.kg, biom.plot$biom.growth.ann.rel, col="gray55", pch=14, cex=0.3)
points(biom.plot$total.biom0.kg, predict(mod), col="red", pch=16, cex=0.3)
points(biom.plot$total.biom0.kg, predict(mod2), col="blue", pch=16, cex=0.3)

x=seq(0, 20, by=0.1)
a=4
b=(-0.5)
plot(x, exp(a+b*x))
