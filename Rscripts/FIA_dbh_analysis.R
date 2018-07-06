library(raster)
library(data.table)
### individual FIA tree data from sites near Boston
dat <- read.csv("data/FIA/MA_Tree_Data.csv")
dat <- read.csv("data/FIA/MA_Tree_Data_ID.csv")
plot(log(dat$DIAM_T0), dat$Adj.Biomass, ylim=c(0,100))
spec <- read.csv("data/FIA/REF_SPECIES.csv")
dat <- merge(x=dat, y=spec[,c("SPCD", "GENUS")], by.x="SPECIES_CD", by.y="SPCD", all.x=T, all.y=F)
dat$GENUS <- as.character(dat$GENUS)
dat$GENUS <- as.factor(dat$GENUS)
table(dat$GENUS)
dat$GENUS.num <- as.numeric(dat$GENUS)

b0 <- -2.48
b1 <- 2.4835 ## these are eastern hardwood defaults
biom.pred <- function(x){exp(b0+(b1*log(x)))}
dat$biom0 <- biom.pred(dat$DIAM_T0)
dat$biom1 <- biom.pred(dat$DIAM_T1)
dat$biom.delt <- dat$biom1-dat$biom0
plot(dat$DIAM_T0, dat$biom.delt, ylim=c(0, 500), col=dat$GENUS) ## some are diving
dat$biom.rel <- (dat$biom.delt/4.8)/dat$biom0
summary(dat$biom.rel) ## median is 1.75% growth
plot(log(dat$DIAM_T0), log(dat$biom.rel), xlim=c(log(5), log(100)), col=dat$GENUS) ## some are diving
dat <- as.data.table(dat)
summary(lm(log(biom.rel)~log(DIAM_T0), data=dat)) ### can't work, negatives 
summary(lm(log(biom.rel)~log(DIAM_T0), data=dat[biom.rel>0,])) ### removing morts and trees that lose biomass 
plot(dat[biom.rel>0, log(DIAM_T0)], dat[biom.rel>0, log(biom.rel)], xlim=c(log(5), log(100)), col=dat$GENUS) ## some are diving

freq <- dat[,length(DIAM_T0), by=GENUS]
dat.big <- dat[GENUS%in%as.character(freq[V1>100, GENUS]) & DIAM_T0>10,]
dat.big[,range(DIAM_T0)]
hist(dat.big$DIAM_T0)
dat.big$GENUS.num <- as.numeric(dat.big$GENUS)
dat.big$GENUS <- as.factor(as.character(dat.big$GENUS))
plot(log(dat.big$DIAM_T0), log(dat.big$biom.rel), col=dat.big$GENUS, xlim=c(2.5, 4.5))
legend(x = 4, y = -1, legend = unique(as.character(dat.big$GENUS)), fill=unique(dat.big$GENUS.num))

names(dat)[2] <- c("TreeID")
names(dat)[3] <- c("PlotID")
dat


range(dat$X)
range(dat$ID
      )
## can we just use the plot data like we did the Andy data? 
## How much %biomass did each plot gain per year as a function of standing biomass in T0, and apply that 

