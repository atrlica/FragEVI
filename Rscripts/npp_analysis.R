library(raster)
library(data.table)
library(tidyverse)
library(lme4)
library(nlme)

#### the summary of the underlying growth models (old)
#####
### final stem growth and areal growth equations
live <- as.data.table(read.csv("processed/fia.live.stem.dbh.growth.csv"))
andy <- as.data.table(read.csv("processed/andy.bai.dbh.pseudo.csv"))
street <- as.data.table(read.csv("processed/boston/street.trees.dbh.csv"))

## this handles the non-linear fit and includes hard/softwood factor, but no potential random effects
mod.fia.final <- nls(growth.ann.rel~exp(a+b*log(DIAM_T0)+as.numeric(type=="S")*c), data=live, start=list(a=0, b=0, c=0)) ## tiny upward correction for softwoods
summary(mod.fia.final)
save(mod.fia.final, file="processed/mod.fia.final.sav")

### this approach handles the random effects but has the unfortunate property of being log-transformed in the response and the predictor
mod.andy.final<- lmer(log(biom.rel.ann)~log(dbh.start)+seg.Edge + 
                          (1+log(dbh.start)|Plot.ID) + 
                          (1+log(dbh.start)|incr.ID) + 
                          (1+log(dbh.start)|interval) +
                          (1+seg.Edge|Plot.ID) +
                          # (1+seg.Edge|incr.ID) + ## there are no individual stems that are in both an edge and interior
                          (1+seg.Edge|interval), 
                        data=andy[dbh.start>=5,],
                        REML=FALSE)
drop1(mod.andy.final, test="Chisq") ## sig in dbh and in edge class (E vs I)
save(mod.andy.final, file="processed/mod.andy.final.sav")

### non-linear street trees -- does not account for any potential random effects
mod.street.final <- nls(npp.ann.rel ~ exp(a + b * log(dbh.2006)), data=street[record.good==1,], start=list(a=0, b=0)) ### OK THIS is the real exponential non-linear model that can handle the negatives
summary(mod.street.final) ## RSE is pretty high 0.33
save(mod.street.final, file="processed/mod.street.final.sav")

# ## compare coefficients
# summary(mod.fia.final)$coefficients #lowest
# summary(mod.andy.final)$coefficients #middle
# summary(mod.street.final)$coefficients #highest

###
## resultant areal-basis growth regressions
live.plot <- as.data.table(read.csv("processed/fia.live.plot.groomed.csv"))
mod.live.plot.final <- nls(biom.growth.ann.hw ~ exp(a + b * log(total.biom0.MgC.ha)),
                                              data=live.plot[hw.frac>0.25,], start=list(a=0, b=0))
summary(mod.live.plot.final) ## all significant
save(mod.live.plot.final, file="processed/mod.live.plot.final.sav")

andy.plot <- read.csv("processed/andy.plot.groomed.csv")
mod.andy.plot.final <- lm(rel.gain.ps~biom.MgC.ha+seg.F, data=andy.plot) #R2 0.67, only biom*edge not significant
save(mod.andy.plot.final, file="processed/mod.andy.plot.final.sav")


### summary table of all the coefficients in play
## plot growth models
mandy.stem <- summary(mod.andy.final)
mfia.stem <- summary(mod.fia.final)
mstreet.stem <- summary(mod.street.final)
mandy.stem.p <- drop1(mod.andy.final, test="Chisq")
stem.context <- c("Rural.stem", "Urban.stem", "Street.stem")
b0.stem <- c(mfia.stem$coefficients[1,1], mandy.stem$coefficients[1,1], mstreet.stem$coefficients[1,1])
b0.stem.se <- c(mfia.stem$coefficients[1,2], mandy.stem$coefficients[1,2], mstreet.stem$coefficients[1,2])
b0.stem.p <- round(c(mfia.stem$coefficients[1,4], NA, mstreet.stem$coefficients[1,4]),3)

b1.stem <- c(mfia.stem$coefficients[2,1], mandy.stem$coefficients[2,1], mstreet.stem$coefficients[2,1])
b1.stem.se <- c(mfia.stem$coefficients[2,2], mandy.stem$coefficients[2,2], mstreet.stem$coefficients[2,2])
b1.stem.p <- round(c(mfia.stem$coefficients[2,4], mandy.stem.p[2,4], mstreet.stem$coefficients[2,4]),3)

b2.stem <- c(mfia.stem$coefficients[3,1], mandy.stem$coefficients[3,1], NA)
b2.stem.se <- c(mfia.stem$coefficients[3,2], mandy.stem$coefficients[3,2], NA)
b2.stem.p <- round(c(mfia.stem$coefficients[3,4], mandy.stem.p[3,4], NA), 3)

b0.stem.tot <- paste(round(b0.stem, 2), " (", round(b0.stem.se, 2), ")", sep="")
b1.stem.tot <- paste(round(b1.stem, 2), " (", round(b1.stem.se, 2), ")", sep="")
b2.stem.tot <- paste(round(b2.stem, 2), " (", round(b2.stem.se, 2), ")", sep="")

sigma <- round(c(mfia.stem$sigma, mandy.stem$sigma, mstreet.stem$sigma), 3)

b0 <- cbind(b0.stem.tot, b0.stem.p)
b1 <- cbind(b1.stem.tot, b1.stem.p)
b2 <- cbind(b2.stem.tot, b2.stem.p)

stem.table <- cbind(stem.context, b0,b1,b2,sigma)

## plot growth models
mfia.plot <- summary(mod.live.plot.final)
mandy.plot <- summary(mod.andy.plot.final)
mandy.plot$coefficients[,1:2] <- mandy.plot$coefficients[,1:2]*100

plot.context <- c("Rural.plot", "Urban.plot")

b0.plot <- c(mfia.plot$coefficients[1,1], mandy.plot$coefficients[1,1])
b0.plot.se <- c(mfia.plot$coefficients[1,2], mandy.plot$coefficients[1,2])
b0.plot.p <- round(c(mfia.plot$coefficients[1,4], mandy.plot$coefficients[1,4]),3)

b0.plot <- c(mfia.plot$coefficients[1,1], mandy.plot$coefficients[1,1])
b0.plot.se <- c(mfia.plot$coefficients[1,2], mandy.plot$coefficients[1,2])
b0.plot.p <- round(c(mfia.plot$coefficients[1,4], mandy.plot$coefficients[1,4]),3)

b1.plot <- c(mfia.plot$coefficients[2,1], mandy.plot$coefficients[2,1])
b1.plot.se <- c(mfia.plot$coefficients[2,2], mandy.plot$coefficients[2,2])
b1.plot.p <- round(c(mfia.plot$coefficients[2,4], mandy.plot$coefficients[2,4]),3)

b2.plot <- c(NA, mandy.plot$coefficients[3,1])
b2.plot.se <- c(NA, mandy.plot$coefficients[3,2])
b2.plot.p <- round(c(NA, mandy.plot$coefficients[3,4]),3)

b0.plot.tot <- paste(round(b0.plot, 2), " (", round(b0.plot.se, 2), ")", sep="")
b1.plot.tot <- paste(round(b1.plot, 2), " (", round(b1.plot.se, 2), ")", sep="")
b2.plot.tot <- paste(round(b2.plot, 2), " (", round(b2.plot.se, 2), ")", sep="")

sigma <- round(c(mfia.plot$sigma, mandy.plot$sigma), 3)

b0 <- cbind(b0.plot.tot, b0.plot.p)
b1 <- cbind(b1.plot.tot, b1.plot.p)
b2 <- cbind(b2.plot.tot, b2.plot.p)

plot.table <- cbind(plot.context, b0,b1,b2,sigma)

write.csv(rbind(stem.table, plot.table), "processed/stem.plot.model.summstats.csv")

# ### just for context -- are the stem growth equations really indicated to be different depending on context?
# 
# dbh <- c(street[record.good==1, dbh.2006],
#          andy[,dbh.start],
#          live[,DIAM_T0])
# growth <- c(street[record.good==1, npp.ann.rel],
#             andy[, biom.rel.ann],
#             live[, growth.ann.rel])
# dat <- cbind(dbh, growth, 
#              c(rep("street", nrow(street[record.good==1,])),
#                rep("urban", nrow(andy)),
#                rep("fia", nrow(live))))
# dat <- as.data.frame(dat)
# dat$dbh <- as.numeric(as.character(dat$dbh))
# dat$growth <- as.numeric(as.character(dat$growth))
# dat <- as.data.table(dat)
# plot(dat$dbh, dat$growth)
# plot(dat[,log(dbh)], dat[,log(growth)], col=as.numeric(dat[,V3]))
# a <- lm(log(growth)~log(dbh), data=dat[growth>0,]) 
# summary(a) ## R2 0.18
# b <- lm(log(growth)~log(dbh)*V3, data=dat[growth>0,])
# summary(b) ## R2 0.52 -- all context factors significant!!!
#####

###
### Canopy configuration across study area
#####
## canopy and biomass configurations in this beast -- 1m data
library(raster)
library(data.table)
library(tidyverse)
library(lme4)
library(nlme)

a <- Sys.time()
aoi <- raster("processed/boston/bos.aoi.tif") 
lulc <- raster("processed/boston/bos.lulc.lumped.tif") %>% crop(aoi) %>% as.vector()
can <- raster("processed/boston/bos.can.redux.tif") %>% crop(aoi) %>% as.vector()
biom <- raster("data/dataverse_files/bostonbiomass_1m.tif") %>% crop(aoi) %>% as.vector()
ed10 <- raster("processed/boston/bos.ed10m.redux.tif") %>% crop(aoi) %>% as.vector()
isa <- raster("processed/boston/bos.isa.RR2.tif") %>% crop(aoi) %>% as.vector()
# aoi <-  as.vector(aoi)
Sys.time()-a ## this chunk takes 7 min to load

# aoi.tot <- sum(aoi[!is.na(aoi)])/1E4 ## 12k ha
aoi.tot <- length(lulc[!is.na(lulc)])/1E4 ## 12.4k ha *** more restrictive boundary of AOI
can.area.tot <- sum(can[!is.na(lulc)], na.rm=T)/1E4 # with new biomass>0 map, down to 3.1k ha (was 3.9k ha)
can.area.tot/aoi.tot ## down to 25.3% canopy area (was 31.8% canopy area)
ed.area.tot <- sum(ed10[!is.na(lulc)], na.rm=T)/1E4 ## down to 2.7k ha (was 3.2k ha)
ed.area.tot/can.area.tot ## up to 84.7% edge canopy (was 80.9% edge canopy in total area)
biom.tot <- sum(biom[!is.na(lulc)], na.rm=T)/(2*1E6) ### 357 MgCx1000
isa.tot <- sum(isa[!is.na(lulc)], na.rm=T)/1E4 ## 7127 ha

### land cover by lulc
lulc.area.tot <- numeric()
lulc.can.tot <- numeric()
lulc.ed.tot <- numeric()
lulc.biom.tot <- numeric()
lulc.isa.tot <- numeric()
for(i in 1:6){
  tmp <- length(lulc[lulc==i & !is.na(lulc)])/1E4 ## ha
  lulc.area.tot <- c(lulc.area.tot, tmp)
  
  tmp <- sum(can[lulc==i & !is.na(lulc)], na.rm=T)/1E4 ## ha
  lulc.can.tot <- c(lulc.can.tot, tmp)
  
  tmp <- sum(ed10[lulc==i & !is.na(lulc)], na.rm=T)/1E4 ## ha
  lulc.ed.tot <- c(lulc.ed.tot, tmp)
  
  tmp <- sum(biom[lulc==i & !is.na(lulc)], na.rm=T)/(2*1E6) ## in MgC x 1000
  lulc.biom.tot <- c(lulc.biom.tot, tmp)
  
  tmp <- sum(isa[lulc==i & !is.na(lulc)], na.rm=T)/1E4 ## ha
  lulc.isa.tot <- c(lulc.isa.tot, tmp)
}
fat <- cbind(c(1:6), lulc.area.tot, lulc.can.tot, lulc.ed.tot, lulc.biom.tot, lulc.isa.tot)

write.csv(fat, "processed/boston/bos.lulc.1mcover.summary.csv")
#####


###
### Comparison of pixel median distributions in the different raw maps
#####
## median value maps to look at distributions
street.med <- raster("processed/results/street/bos.street.V6.npp.med.tif")
andy.med <- raster("processed/results/andy/bos.andy.forest.V5.npp.tif")
fia.for.med <- raster("processed/results/fia/npp.FIA.empirV6.forest.median.tif")
fia.grd.med <- raster("processed/results/fia/npp.FIA.empirV6.ground.median.tif")

biom <- raster("processed/boston/bos.biom30m.tif")
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- crop(biom, aoi)

dat <- data.frame(biom=getValues(biom),
                     aoi=getValues(aoi),
                     street.med=getValues(street.med),
                     andy.med=getValues(andy.med),
                     fia.for.med=getValues(fia.for.med),
                     fia.grd.med=getValues(fia.grd.med))
dat <- as.data.table(dat)
### get growth factors
dat[, street.fact:=street.med/biom]
dat[, andy.fact:=andy.med/biom]
dat[, fia.for.fact:=fia.for.med/biom]
dat[, fia.grd.fact:=fia.grd.med/biom]
dat[aoi>800 & !is.na(biom) & biom>10, .(range(street.fact, na.rm=T),
                                        range(andy.fact, na.rm=T),
                                        range(fia.for.fact, na.rm=T),
                                        range(fia.grd.fact, na.rm=T))]

### general realtionship between productivity and biomass
par(mfrow=c(1,1))
plot(dat[aoi>800 & !is.na(biom) & biom>10, biom], 
     dat[aoi>800 & !is.na(biom) & biom>10, street.med])
plot(dat[aoi>800 & !is.na(biom) & biom>10, biom], 
     dat[aoi>800 & !is.na(biom) & biom>10, andy.med])
plot(dat[aoi>800 & !is.na(biom) & biom>10, biom], 
     dat[aoi>800 & !is.na(biom) & biom>10, fia.for.med])
plot(dat[aoi>800 & !is.na(biom) & biom>10, biom], 
     dat[aoi>800 & !is.na(biom) & biom>10, fia.grd.med])

### same thing, do it with growth factors
plot(dat[aoi>800 & !is.na(biom) & biom>10, biom], 
     dat[aoi>800 & !is.na(biom) & biom>10, street.fact], ylim=c(0, 0.1))
plot(dat[aoi>800 & !is.na(biom) & biom>10, biom], 
     dat[aoi>800 & !is.na(biom) & biom>10, andy.fact])
plot(dat[aoi>800 & !is.na(biom) & biom>10, biom], 
     dat[aoi>800 & !is.na(biom) & biom>10, fia.for.fact])
plot(dat[aoi>800 & !is.na(biom) & biom>10, biom], 
     dat[aoi>800 & !is.na(biom) & biom>10, fia.grd.fact])

### compare the histograms
lims <- c(0, 0.08)
par(mfrow=c(2,2))
hist(dat[aoi>800 & !is.na(biom) & biom>10 & street.fact<0.08, street.fact], xlim=lims)
hist(dat[aoi>800 & !is.na(biom) & biom>10, andy.fact], xlim=lims)
hist(dat[aoi>800 & !is.na(biom) & biom>10, fia.for.fact], xlim=lims)
hist(dat[aoi>800 & !is.na(biom) & biom>10, fia.grd.fact], xlim=lims)
## ok... the distributions look about what you'd expect.
## streeet is spiked about 0.05 (higher than the others) and has a long positive tail
## andy ramps up to 0.04-0.06, fia.forest is evenly humped at 0.02, and fia.ground is skewed up against the limit at about 0.04

### look at actual pixel productivity
lims=c(0, 500)
par(mfrow=c(2,2))
hist(dat[aoi>800 & !is.na(biom) & biom>10 & street.med<500 & street.med>0, street.med], xlim=lims)
hist(dat[aoi>800 & !is.na(biom) & biom>10 & andy.med>0, andy.med], xlim=lims)
hist(dat[aoi>800 & !is.na(biom) & biom>10 & fia.for.med>0, fia.for.med], xlim=lims)
hist(dat[aoi>800 & !is.na(biom) & biom>10 & fia.grd.med>0, fia.grd.med], xlim=lims)
## yeah, the andy forest and street just has much longer positive tails -- some very high productivity pixels
#####


### PREPARATION OF A FINAL ANDY+STREET HYBRID MAP
## HYBRID MAP WITH ERROR DISTRIBUTION
#####
street.res <- fread("processed/results/street/streettrees.npp.simulator.v6.results.random.csv")
street.res <- street.res[,-1]
street.res[biom<10 & aoi>800, 7:1006:=0] ### set low biomass to 0 growth
street.res[can<0.001 & aoi>800, 7:1006:=0] ## set anything with <1 pix of canopy to 0
street.res[aoi>800 & is.na(biom), 7:1006:=0] ## set no record of biomass to 0
flags <- as.data.table(read.csv("processed/boston/biom_street/results/streettrees.sim.v6.diagnostics.csv")) ## simulator status per pixel
bad.sim <- flags[num.sims<40, pix.ID] ## 9511 these are all pixels that failed complete-enough simulation
names(street.res)[c(2,4,5,6)] <- c("bos.biom30m", "bos.can30m", "bos.isa30m", "bos.lulc30m.lumped")
# summary(street.res[aoi>800 & pix.ID %in% bad.sim, bos.biom30m]); hist(street.res[aoi>800 & pix.ID %in% bad.sim, bos.biom30m]) ## the bad pixels apparently are mostly in the >20000 range
# summary(street.res[aoi>800 & pix.ID %in% bad.sim, bos.can30m]); hist(street.res[aoi>800 & pix.ID %in% bad.sim, bos.can30m]) ## most are fully canopied
# table(street.res[aoi>800 & pix.ID %in% bad.sim, bos.lulc30m.lumped]) ## most are forest or HD res

### select andy results that are either forest or high biomass or are failed sims
andy.res <- (fread("processed/andy.forest.results.V5.csv"))
forest <- andy.res[lulc==1, pix.ID]
big <- andy.res[biom>=20000, pix.ID]
gimme <- c(forest, big, bad.sim)
gimme <- unique(gimme); length(gimme) ## 17052 pix we are going to swap out with andy results
pick <- andy.res[pix.ID%in%gimme,] 
pick <- pick[,c(8,15:1014)] ## pix.ID and then results
names(pick)[2:1001] <- paste0("npp.iter.", 1:1000, ".hybrid"); dim(pick) ### 17052 pix

## finalize grooming the street pixels 
dim(street.res) ## 354068
street.res <- street.res[bos.lulc30m.lumped!=1,]; dim(street.res) ## 126437
street.res <- street.res[bos.biom30m<20000,]; dim(street.res) ## 122607
street.res <- street.res[!(pix.ID%in%bad.sim),]; dim(street.res) ### 119739
street.res <- street.res[!(pix.ID%in%gimme),]; dim(street.res) ### 119739
street.res <- street.res[,c(1,7:1006)]
names(street.res)[2:1001] <- paste0("npp.iter.", 1:1000, ".hybrid")
hybrid <- rbind(pick, street.res); dim(hybrid) ### 136791 pix; now need to merge these back into the complete map

library(raster)
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- raster("processed/boston/bos.biom30m.tif")
can <- raster("processed/boston/bos.can.redux30m.tif")
isa <- raster("processed/boston/bos.isa30m.tif")
lulc <- raster("processed/boston/bos.lulc30m.lumped.tif")
biom <- crop(biom, aoi)
can <- crop(can, aoi)
isa <- crop(isa, aoi)
lulc <- crop(lulc, aoi)
biom.dat <- as.data.table(cbind(as.data.frame(aoi), 
                                as.data.frame(biom), 
                                as.data.frame(can), 
                                as.data.frame(isa),
                                as.data.frame(lulc)))
biom.dat[, pix.ID:=seq(1:dim(biom.dat)[1])]
biom.dat <- merge(biom.dat, hybrid, by="pix.ID", all.x=T) # here's our map
biom.dat <- biom.dat[order(pix.ID),]; dim(biom.dat) ## 354068 pixels x 1006 cols
fwrite(biom.dat, "processed/results/hybrid.results.V7.csv") ## V7 includes andy V5 and street V6, all with biomass>0 canopy map incorporated

### write a tiff of the median values
# median.na <- function(x){median(x, na.rm=T)}
# rr <- aoi
# pix.med <- apply(as.matrix(biom.dat[,7:1006]), MARGIN=1, FUN=median.na)
# rr <- setValues(rr, pix.med)
# writeRaster(rr, filename="processed/results/hybrid.V7.median.tif", format="GTiff", overwrite=T)
# plot(rr)
#####

###
### analysis of model results
#####
library(raster)
library(data.table)
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- raster("processed/boston/bos.biom30m.tif")
can <- raster("processed/boston/bos.can.redux30m.tif")
isa <- raster("processed/boston/bos.isa30m.tif")
lulc <- raster("processed/boston/bos.lulc30m.lumped.tif")
biom <- crop(biom, aoi)
can <-  crop(can, aoi)
isa <- crop(isa, aoi)
lulc <- crop(lulc, aoi)
biom.dat <- as.data.table(cbind(as.data.frame(aoi), 
                                as.data.frame(biom), 
                                as.data.frame(can), 
                                as.data.frame(isa),
                                as.data.frame(lulc)))
biom.dat[, pix.ID:=seq(1:dim(biom.dat)[1])]
names(biom.dat)[1:5] <- c("aoi", "biom", "can", "isa", "lulc")

## put pixel medians into biom.dat
biom.dat[,fia.ground.med:=getValues(raster("processed/results/fia/npp.FIA.empirV6.ground.median.tif"))]
biom.dat[,fia.forest.med:=getValues(raster("processed/results/fia/npp.FIA.empirV6.forest.median.tif"))]
biom.dat[,hybrid.med:=getValues(raster("processed/results/hybrid.V7.median.tif"))]
biom.dat[,fia.ground.med.MgC.ha:=((fia.ground.med/2000)/aoi)*1E4]
biom.dat[,fia.forest.med.MgC.ha:=((fia.forest.med/2000)/aoi)*1E4]
biom.dat[,hybrid.med.MgC.ha:=((hybrid.med/2000)/aoi)*1E4]
# hist(biom.dat[aoi>800 & !is.na(biom),hybrid.med.MgC.ha]) ## much longer positive tail
# hist(biom.dat[aoi>800 & !is.na(biom),fia.forest.med.MgC.ha])

## compile a table of pixel median ranges
l <- biom.dat[aoi>800 & !is.na(biom), round(quantile(hybrid.med.MgC.ha, probs=0.05, na.rm=T), 1), by=lulc]
l <- l[order(lulc),]
m <- biom.dat[aoi>800 & !is.na(biom), round(quantile(hybrid.med.MgC.ha, probs=0.5, na.rm=T), 1), by=lulc]
m <- m[order(lulc),]
h <- biom.dat[aoi>800 & !is.na(biom), round(quantile(hybrid.med.MgC.ha, probs=0.95, na.rm=T), 1), by=lulc]
h <- h[order(lulc),]

fin <- cbind(c("Forest", "Dev", "HDRes", "LDRes", "OVeg", "Water", "NA"),
      paste0(m[,V1], " (", l[,V1], "-", h[,V1], ")"))

l <- biom.dat[aoi>800 & !is.na(biom), round(quantile(fia.forest.med.MgC.ha, probs=0.05, na.rm=T), 1), by=lulc]
l <- l[order(lulc),]
m <- biom.dat[aoi>800 & !is.na(biom), round(quantile(fia.forest.med.MgC.ha, probs=0.5, na.rm=T), 1), by=lulc]
m <- m[order(lulc),]
h <- biom.dat[aoi>800 & !is.na(biom), round(quantile(fia.forest.med.MgC.ha, probs=0.95, na.rm=T), 1), by=lulc]
h <- h[order(lulc),]

fin <- cbind(fin,
             paste0(m[,V1], " (", l[,V1], "-", h[,V1], ")"))

fin <- rbind(fin, c("Total",
                    paste0(biom.dat[aoi>800 & !is.na(biom), round(quantile(hybrid.med.MgC.ha, probs=0.5, na.rm=T), 1)],
                           " (",
                           biom.dat[aoi>800 & !is.na(biom), round(quantile(hybrid.med.MgC.ha, probs=0.05, na.rm=T), 1)],
                           "-",
                           biom.dat[aoi>800 & !is.na(biom), round(quantile(hybrid.med.MgC.ha, probs=0.95, na.rm=T), 1)],
                           ")"),
                    paste0(biom.dat[aoi>800 & !is.na(biom), round(quantile(fia.forest.med.MgC.ha, probs=0.5, na.rm=T), 1)],
                           " (",
                           biom.dat[aoi>800 & !is.na(biom), round(quantile(fia.forest.med.MgC.ha, probs=0.05, na.rm=T), 1)],
                           "-",
                           biom.dat[aoi>800 & !is.na(biom), round(quantile(fia.forest.med.MgC.ha, probs=0.95, na.rm=T), 1)],
                           ")"))
             
)
write.csv(fin, "processed/results/pix.med.lulc.spread.csv")

### compile the whole map sum distributions one model at a time (not enough memory to load up every shits at once)
median.na <- function(x){median(x, na.rm=T)}
sum.na <- function(x){sum(x, na.rm=T)}
files <- c("processed/results/hybrid.results.V7.csv", "processed/results/fia/npp.FIA.empirV6.forest.csv", "processed/results/fia/npp.FIA.empirV6.ground.csv")
nms <- c("hybrid", "fia.forest", "fia.ground")
burp <- data.frame(c("Forest", "Dev", "HDRes", "LDRes", "Oveg", "Water", "Total"))
for(t in 1:3){ ## perform this routine for every map model
  tmp <- fread(files[t])
  names(tmp)[grep(names(tmp), pattern="*lulc*")] <- "lulc"
  names(tmp)[grep(names(tmp), pattern="*aoi*")] <- "aoi"
  names(tmp)[grep(names(tmp), pattern="*biom*")] <- "biom"
  p <- character()
  for(l in 1:6){ ## split data up by lulc
    print(paste("getting totals for lulc =", l))
    tmp.l <- tmp[lulc==l & aoi>800 & !is.na(biom),]
    l.med <- apply(as.matrix(tmp.l[, min(grep(names(tmp), pattern="iter")):max(grep(names(tmp), pattern="iter"))]), MARGIN=2, FUN=sum.na)
    j <- paste0(round(median(l.med, na.rm=T)/(2000*1000),1),
           " (",
           round(quantile(l.med, probs=0.05, na.rm=T)/(2000*1000), 1),
           "-", 
           round(quantile(l.med, probs=0.95, na.rm=T)/(2000*1000), 1),
           ")")
    p <- rbind(p, j)
  }
  print(paste("getting model totals for", nms[t]))
  mod.med <- apply(as.matrix(tmp[aoi>800 & !is.na(biom), min(grep(names(tmp), pattern="iter")):max(grep(names(tmp), pattern="iter"))]), MARGIN=2, FUN=sum.na) ## the per-realization map sum of the model (x1000 per model)
  p <- rbind(p, 
             paste0(round(median(mod.med, na.rm=T)/(2000*1000),1),
                    " (",
                    round(quantile(mod.med, probs=0.05, na.rm=T)/(2000*1000), 1),
                    "-", 
                    round(quantile(mod.med, probs=0.95, na.rm=T)/(2000*1000), 1),
                    ")"))
  burp <- cbind(burp, p)
}
burp <- as.data.frame(burp)
names(burp) <- c("LULC", "Hybrid", "FIA.forest", "FIA.ground")
write.csv(burp, paste0("processed/results/npp.model.sums.spread.csv"))
read.csv(paste0("processed/results/npp.model.sums.spread.csv"))

### Results for ANDY FOREST, error distributed
#####
andy <- as.data.table(read.csv("processed/andy.forest.results.V4.csv"))
# andy <- as.data.table(read.csv("processed/andy.forest.results.V2.csv")) ## this is the 7000 MgC version
# andy.retest <- as.data.table(read.csv("processed/andy.forest.results.V2-retest.csv")) ## this is the 7000 MgC version
sum.na <- function(x){sum(x, na.rm=T)}
andy.npp.tot <- apply(andy[aoi>800,16:115], MARGIN=2, FUN=sum.na) 
# andy.npp.tot.retest <- apply(andy.retest[aoi>800,16:115], MARGIN=2, FUN=sum.na) 
# andy.npp.tot.newmod <- apply(andy.newmod[aoi>800,16:115], MARGIN=2, FUN=sum.na) 
hist((andy.npp.tot/2000))
# hist((andy.npp.tot.retest/2000))## no variance, something is broken in the model
# hist((andy.npp.tot.newmod/2000))
median((andy.npp.tot/2000)) ## about 8.5k
# median((andy.npp.tot.retest/2000)) ## about 6.2k
# median((andy.npp.tot.newmod/2000)) ## about 12.4k

## fuck trying to do a raster of the andy median

### model realizations by lulc
andy.lulc <- data.frame(1:100)
for(l in 1:6){
  tmp <- apply(andy[aoi>800 & lulc==l,16:115], MARGIN=2, FUN=sum.na)
  andy.lulc <- cbind(andy.lulc, tmp)
}
names(andy.lulc) <- c("iter", paste("andy.lulc", 1:6, ".npp.tot", sep=""))
rownames(andy.lulc) <- NULL
write.csv(andy.lulc, "processed/andy.v4.lulc.results.csv")

### summary stats
andy.lulc <- as.data.table(read.csv("processed/andy.v4.lulc.results.csv"))
andy.lulc[,X:=NULL]

## quantiles by lulc 
rangey <- function(x){quantile(x, probs=c(0.05, 0.5, 0.95))}
tot <- round(rangey(andy.npp.tot)/2000/1000,1)
lulc.tot <- round(apply(andy.lulc[,2:7], FUN=rangey, MARGIN=2)/2000/1000, 1)
results.andy <- character()
for(i in 1:6){
  results.andy <- c(results.andy,
                          paste0(lulc.tot[2,i], " (", lulc.tot[1,i], "-", lulc.tot[3,i], ")"))
}
t <- cbind(t, 
           c(results.andy, paste0(tot[2], " (", tot[1], "-", tot[3], ")")))
colnames(t)[4] <- "results.andy"
#####


### random bullshit trying to summarize and compare results of different approaches
#####
# ### summary table by lulc
# tot.area <- npp.dat[aoi>800 & !is.na(lulc), sum(aoi)]
# fox <- npp.dat[aoi>800 & !is.na(lulc), 
#                .((sum(aoi)/tot.area), 
#                  (median(can.frac, na.rm=T)), 
#                  (sum(bos.biom30m, na.rm=T)/1000), 
#                  (sum(bos.biom30m, na.rm=T)/1000)/(sum(aoi)/1E4)), by=lulc]
# names(fox) <- c("lulc", "area.frac", "med.can.frac", "tot.biomass.Mg", "biomass.Mg.ha")
# npp.dat[aoi>800 & !is.na(lulc), ## whole area
#                .((sum(aoi)/tot.area), 
#                  (median(can.frac, na.rm=T)), 
#                  (sum(bos.biom30m, na.rm=T)/1000), 
#                  (sum(bos.biom30m, na.rm=T)/1000)/(sum(aoi)/1E4))]
# 


# ## productivity summaries
# hyb.tot <- npp.dat[aoi>800 & !is.na(lulc), sum(hyb.npp, na.rm=T)/(1000*2)]
# hen <- npp.dat[aoi>800 & !is.na(lulc),
#                .((sum(hyb.npp, na.rm=T)/2000),
#                  (sum(hyb.npp, na.rm=T)/2000)/(sum(aoi, na.rm=T)/1E4),
#                  (sum(hyb.npp, na.rm=T)/2000)/hyb.tot), by=lulc]
# names(hen) <- c("lulc", "hyb.npp.MgC", "hyb.npp.MgC.ha", "hyb.npp.frac.tot")
# house <- merge(fox, hen, by="lulc")
# 
# flop.tot <- npp.dat[aoi>800 & !is.na(lulc), sum(fia.npp.ann.forest, na.rm=T)/(1000*2)]
# flop <- npp.dat[aoi>800 & !is.na(lulc), 
#                 .((sum(fia.npp.ann.forest, na.rm=T)/2000),
#                   (sum(fia.npp.ann.forest, na.rm=T)/2000)/(sum(aoi, na.rm=T)/1E4),
#                   (sum(fia.npp.ann.forest, na.rm=T)/2000)/flop.tot), by=lulc]
# names(flop) <- c("lulc", "fia.npp.MgC", "fia.npp.MgC.ha", "fia.npp.frac.tot")
# house <- merge(house, flop, by="lulc")
# write.csv(house, "processed/boston/results/summary.hybrid.fia.npp.csv")

## productivity across all lulc
# npp.dat[aoi>800, (sum(fia.npp.ann.ground, na.rm=T)/(1000*2))/sum(aoi, na.rm=T)]*1E4 ### 1.11 MgC/ha
# npp.dat[aoi>800, (sum(fia.npp.ann.forest, na.rm=T)/(1000*2))/sum(aoi, na.rm=T)]*1E4 ### 0.37 MgC/ha
# npp.dat[aoi>800, (sum(fia.npp.ann.perv, na.rm=T)/(1000*2))/sum(aoi, na.rm=T)]*1E4 ### 0.43 MgC/ha
# npp.dat[aoi>800, (sum(fia.empir.npp.kg.hw.forest, na.rm=T)/(1000*2))/sum(aoi, na.rm=T)]*1E4 ### 0.7 MgC/ha
# npp.dat[aoi>800, (sum(andy.npp.tot, na.rm=T)/(1000*2))/sum(aoi, na.rm=T)]*1E4 ### 1.11 MgC/ha
# npp.dat[aoi>800, (sum(st.med.ann.npp.all, na.rm=T)/(1000*2))/sum(aoi, na.rm=T)]*1E4 ### 1.07 MgC/ha (excluding non-retrieves)
# npp.dat[aoi>800, (sum(hyb.npp, na.rm=T)/(1000*2))/sum(aoi, na.rm=T)]*1E4 ### 1.17 MgC/ha
# npp.dat[aoi>800, (sum(hyb.npp, na.rm=T)/(1000*2))] ### 14,393
# 


### only thing is to decide whether or not to swap the high-biomass non-forest sim results with the andy forest equivalents
## I thikn this is a yes: Anything >20000 kg/pix is a "forest" in the andy sense
# 
# #### BRING IN NPP estimates based on the FIA empirical records
# npp.dat <- as.data.table(read.csv("processed/npp.estimates.V1.csv"))
# npp.tab <- npp.dat[aoi>800 & !is.na(lulc), .(round(sum(fia.empir.npp.kg.hw.ground, na.rm=T)/(1000*2), digits=0),
#                               round(sum(fia.empir.npp.kg.hw.forest, na.rm=T)/(1000*2), digits=0),
#                               round(sum(fia.empir.npp.kg.hw.perv, na.rm=T)/(1000*2), digits=0),
#                               round(sum(hyb.npp, na.rm=T)/(1000*2), digits=0)), by=lulc]
# npp.tab$lulc.name <- c("Dev", "Other Veg", "Forest", "Water", "HD Resid", "LD Resid")
# names(npp.tab) <- c("lulc", "FIA (Ground)", "FIA (canopy)", "FIA (pervious)", "Street+Urban Forest", "Type")
# npp.tab <- cbind(npp.tab[,2:5, with=F], npp.tab$Type)
# npp.tot <- apply(npp.tab[,1:4, with=F], 2, sum)
# npp.tot <- matrix(nrow=1, ncol=5, c(npp.tot, "Total"))
# npp.tab <- as.matrix(npp.tab)
# npp.tab <- rbind(npp.tab, npp.tot)
# npp.tab <- as.data.frame(cbind(npp.tab[,5], npp.tab[,1:4]))
# names(npp.tab)[1] <- "Type"
# 
# hist(npp.dat[bos.biom30m>20000 & aoi>800, can.frac])
# hist(npp.dat[lulc==1 & aoi>800, can.frac]) ## looks a lot like forest distribution generally
# summary(npp.dat[bos.biom30m>20000 & aoi>800, can.frac]) ## anything over 20k kg looks a lot like forest
# summary(npp.dat[lulc==1 & aoi>800, can.frac])

#########
### summary table (TABLE 1)
npp.dat <- as.data.table(read.csv("processed/npp.estimates.V4.csv"))
npp.dat[,biom.MgC.ha.ground:=((bos.biom30m/2000)/aoi)*1E4] ## per-pixel ground biomass density
npp.dat[,biom.MgC.ha.can:=((bos.biom30m/2000)/(aoi*can.frac))*1E4] ## per-pixel canopy biomass density
### these NPP flux densities are all in units of MgC/ha-ground (physical reality), even though they may differ in how *biom density* was calcualted prior to NPP estimation
npp.dat[,fia.empir.npp.ground.MgC.ha:=((fia.empir.npp.kg.hw.ground/2000)/aoi)*1E4] ## actual per pixel MgC/ha, ground basis calc
npp.dat[,fia.empir.npp.can.MgC.ha:=((fia.empir.npp.kg.hw.forest/2000)/aoi)*1E4] ## actual per pixel MgC/ha, canopy basis calc
# npp.dat[,hyb.npp.MgC.ha.ground:=((hyb.npp.mod.lin/2000)/aoi)*1E4]
npp.dat[,hyb.npp.MgC.ha:=((hyb.npp.mod.lin/2000)/(aoi))*1E4]
# npp.dat[,andy.npp.MgC.ha.ground:=((andy.npp.tot.mod.lin/2000)/aoi)*1E4]
npp.dat[,andy.npp.MgC.ha:=((andy.npp.tot.mod.lin/2000)/(aoi))*1E4]

## basic by-LULC summary table
tab1 <- npp.dat[aoi>800 & !is.na(lulc),.(round((sum(aoi, na.rm=T)/1E4), 0), ## tot area
                         round(median(can.frac, na.rm=T),2)*100, ## % canopy
                         round((sum(bos.biom30m, na.rm=T)/2000)/1000,0), ## biomass total
                         round(median(biom.MgC.ha.ground, na.rm=T),0), ## biom MgC/ha-ground
                         round(median(biom.MgC.ha.can, na.rm=T),0), ## biom MgC/ha-canopy
                         round(sum(fia.empir.npp.kg.hw.ground/2000, na.rm=T),0), ## Tot NPP FIA, ground basis calc
                         round(sum(fia.empir.npp.kg.hw.forest/2000, na.rm=T),0), ## Tot NPP FIA, canopy basis calc
                         round(sum(andy.npp.tot.mod.lin/2000, na.rm=T),0), ## tot NPP andy forest 
                         round(sum(hyb.npp.mod.lin/2000, na.rm=T),0), ## tot NPP HYBRID
                         ###### per-pixel MgC/ha/yr rates
                         round(median(fia.empir.npp.ground.MgC.ha, na.rm=T),1), ## FIA NPP MgC/ha-ground -- canopy basis calc
                         round(median(fia.empir.npp.can.MgC.ha, na.rm=T),1), ### FIA NPP MgC/ha-can -- canopy basis calc
                         # round(median(andy.npp.MgC.ha.ground, na.rm=T),1), ## andy forest NPP MgC/ha-ground
                         round(median(andy.npp.MgC.ha, na.rm=T),1), ### andy forest NPP MgC/ha-can
                         # round(median(hyb.npp.MgC.ha, na.rm=T),1), ## HYBRID NPP MgC/ha-ground
                         round(median(hyb.npp.MgC.ha, na.rm=T),1)), by=lulc] ## HYBRID NPP MgC/ha-can

names(tab1)[2:14] <- c("area.ha", "Mcan", "Tbiom.kMgC", "Mbiom.MgC.ha.ground", "Mbiom.MgC.ha.can", 
                      "Tfia.empir.npp.ground.MgC", "Tfia.empir.npp.forest.MgC", "Tandy.npp.MgC", "Thyb.npp.mod.lin.MgC",
                      "Mfia.npp.ground.MgC.ha", "Mfia.npp.can.MgC.ha", "Mandy.npp.MgC.ha", "Mhyb.npp.MgC.ha")
# View(tab1)
# area.tot <- npp.dat[aoi>800, sum(aoi, na.rm=T)/1E4] ## total area in ha
# biom.tot <- npp.dat[aoi>800, (sum(bos.biom30m, na.rm=T)/2000)/1000] ## total biomass in 1000 MgC
# hyb.tot <- npp.dat[aoi>800, sum(hyb.npp.mod.lin, na.rm=T)/2000]
# fia.empir.tot <- npp.dat[aoi>800, sum(fia.empir.npp.kg.hw.forest, na.rm=T)/2000]
# npp.dat[aoi>800, median(can.frac, na.rm=T)] ## whole-city median canopy
# npp.dat[aoi>800, median(fia.empir.npp.MgC.ha.ground, na.rm=T)]
# npp.dat[aoi>800, median(fia.empir.npp.MgC.ha.can, na.rm=T)]
# npp.dat[aoi>800, median(andy.npp.MgC.ha.ground, na.rm=T)]
# npp.dat[aoi>800, median(andy.npp.MgC.ha.can, na.rm=T)]
# npp.dat[aoi>800, median(hyb.npp.MgC.ha.ground, na.rm=T)]
# npp.dat[aoi>800, median(hyb.npp.MgC.ha.can, na.rm=T)]




### alternate layout: LULC | total NPP: FIAground, FIAcan, Andy-can, Hybrid | Median (95%ile) MgC/ha-ground: FIAcan, Hybrid
fia.MgC.ha.low <- npp.dat[aoi>800 & !is.na(lulc), quantile(fia.empir.npp.can.MgC.ha, probs=c(0.05),  na.rm=T), by=lulc]
fia.MgC.ha.hi <- npp.dat[aoi>800 & !is.na(lulc), quantile(fia.empir.npp.can.MgC.ha, probs=c(0.95),  na.rm=T), by=lulc]
hyb.MgC.ha.low <- npp.dat[aoi>800 & !is.na(lulc), quantile(hyb.npp.MgC.ha, probs=c(0.05),  na.rm=T), by=lulc]
hyb.MgC.ha.hi <- npp.dat[aoi>800 & !is.na(lulc), quantile(hyb.npp.MgC.ha, probs=c(0.95),  na.rm=T), by=lulc]


tab1.f <- data.frame(cbind(as.character(tab1$lulc),
                           as.character(tab1$Tfia.empir.npp.ground.MgC),
                           as.character(tab1$Tfia.empir.npp.forest.MgC),
                           as.character(tab1$Tandy.npp.MgC),
                           as.character(tab1$Thyb.npp.mod.lin.MgC),
                           paste(tab1$Mfia.npp.can.MgC.ha, " (",
                                 round(fia.MgC.ha.low$V1, 1), "-",
                                 round(fia.MgC.ha.hi$V1, 1), ")", sep=""),
                           paste(tab1$Mhyb.npp.MgC.ha, " (",
                                 round(hyb.MgC.ha.low$V1, 1), "-",
                                 round(hyb.MgC.ha.hi$V1, 1), ")", sep="")))

tots <- c("Total", 
  round(npp.dat[aoi>800 & !is.na(lulc), sum(fia.empir.npp.kg.hw.ground/2000, na.rm=T)], 0),
  round(npp.dat[aoi>800 & !is.na(lulc), sum(fia.empir.npp.kg.hw.forest/2000, na.rm=T)], 0),
  round(npp.dat[aoi>800 & !is.na(lulc), sum(andy.npp.tot.mod.lin/2000, na.rm=T)], 0),
  round(npp.dat[aoi>800 & !is.na(lulc), sum(hyb.npp.mod.lin/2000, na.rm=T)], 0),
  paste(round(npp.dat[aoi>800 & !is.na(lulc), median(fia.empir.npp.can.MgC.ha, na.rm=T)], 1),
        " (", round(npp.dat[aoi>800 & !is.na(lulc), quantile(fia.empir.npp.can.MgC.ha, probs=0.05, na.rm=T)], 1),
        "-", round(npp.dat[aoi>800 & !is.na(lulc), quantile(fia.empir.npp.can.MgC.ha, probs=0.95, na.rm=T)], 1), ")", sep=""),
  paste(round(npp.dat[aoi>800 & !is.na(lulc), median(hyb.npp.MgC.ha, na.rm=T)], 1),
        " (", round(npp.dat[aoi>800 & !is.na(lulc), quantile(hyb.npp.MgC.ha, probs=0.05, na.rm=T)], 1),
        "-", round(npp.dat[aoi>800 & !is.na(lulc), quantile(hyb.npp.MgC.ha, probs=0.95, na.rm=T)], 1), ")", sep=""))

## fuck it export
write.csv(tab1.f, "processed/boston/results/Table1_totals_only.csv")

### new layout: LULC, AREA, CANOPY, TOTAL BIOM | FIA(ground), HYBRID(ground) | FIA(can), HYBRID(can) | Tot.FIA, Tot.HYBRID
tab1.f <- data.frame(cbind(tab1$lulc,
                           paste(tab1$area.ha, " (", round((tab1$area/area.tot)*100,0), "%)", sep=""),
                           paste(tab1$med.can, "%", sep=""),
                           paste(tab1$Tbiom.kMgC, " (", round((tab1$Tbiom.kMgC/biom.tot)*100,0), "%)", sep=""),
                           tab1$fia.npp.MgC.ha.ground, tab1$hyb.npp.MgC.ha.ground,
                           tab1$fia.npp.MgC.ha.can, tab1$hyb.npp.MgC.ha.can,
                           paste(tab1$fia.empir.npp.MgC, " (", round((tab1$fia.empir.npp.MgC/fia.empir.tot)*100,0), "%)", sep=""),
                           paste(tab1$hyb.npp.mod.lin.MgC, " (", round((tab1$hyb.npp.mod.lin.MgC/hyb.tot)*100,0), "%)", sep="")
))
names(tab1.f) <- c("LULC",
                   "area.ha(%)",
                   "med.can",
                   "Tot.biom.kMgC(%)",
                   "FIA.MgC.ha.gnd",
                   "HYB.MgC.ha.gnd",
                   "FIA.MgC.ha.can",
                   "HYB.MgC.ha.can",
                   "FIA.tot.MgC(%)",
                   "HYB.tot.MgC(%)")

write.csv(tab1.f, "processed/boston/results/lulc_NPP_summary_simp1.csv")
aa <- read.csv("processed/boston/results/lulc_NPP_summary_simp1.csv")
View(aa)




### new layout: LULC, AREA, CANOPY, TOTAL BIOM | FIA(ground), HYBRID(ground) | FIA(can), HYBRID(can) | Tot.FIA, Tot.HYBRID
tab1.f <- data.frame(cbind(tab1$lulc,
                     paste(tab1$area.ha, " (", round((tab1$area/area.tot)*100,0), "%)", sep=""),
                     paste(tab1$med.can, "%", sep=""),
                     paste(tab1$Tbiom.kMgC, " (", round((tab1$Tbiom.kMgC/biom.tot)*100,0), "%)", sep=""),
                     tab1$fia.npp.MgC.ha.ground, tab1$hyb.npp.MgC.ha.ground,
                     tab1$fia.npp.MgC.ha.can, tab1$hyb.npp.MgC.ha.can,
                     paste(tab1$fia.empir.npp.MgC, " (", round((tab1$fia.empir.npp.MgC/fia.empir.tot)*100,0), "%)", sep=""),
                     paste(tab1$hyb.npp.mod.lin.MgC, " (", round((tab1$hyb.npp.mod.lin.MgC/hyb.tot)*100,0), "%)", sep="")
                     ))
names(tab1.f) <- c("LULC",
                   "area.ha(%)",
                   "med.can",
                   "Tot.biom.kMgC(%)",
                   "FIA.MgC.ha.gnd",
                   "HYB.MgC.ha.gnd",
                   "FIA.MgC.ha.can",
                   "HYB.MgC.ha.can",
                   "FIA.tot.MgC(%)",
                   "HYB.tot.MgC(%)")

write.csv(tab1.f, "processed/boston/results/lulc_NPP_summary_simp1.csv")
aa <- read.csv("processed/boston/results/lulc_NPP_summary_simp1.csv")
View(aa)
## final layout: LULC, area, med.can, total biomass, biomass density | FIA NPP results | HYB NPP RESULTS
tab1.f <- data.frame(cbind(tab1$lulc, 
                           paste(tab1$area.ha, " (", 
                                 round((tab1$area.ha/area.tot)*100,0),
                                 "%)", sep="")))
tab1.f <- cbind(tab1.f, 
                paste(tab1$med.can, "%", sep=""),
                paste(tab1$Tbiom.kMgC, " (", 
                      round((tab1$Tbiom.kMgC/biom.tot)*100,0),
                      "%)", sep=""))
tab1.f <- cbind(tab1.f, tab1$biom.MgC.ha.ground, tab1$biom.MgC.ha.can)
tab1.f <- cbind(tab1.f,
                paste(tab1$fia.empir.npp.MgC, " (",
                      round((tab1$fia.empir.npp.MgC/fia.empir.tot)*100, 0),
                "%)", sep=""))
tab1.f <- cbind(tab1.f,
                round(tab1$fia.npp.MgC.ha.ground, 2),
                round(tab1$fia.npp.MgC.ha.can, 2))
tab1.f <- cbind(tab1.f,
                paste(tab1$hyb.npp.mod.cap.MgC, " (",
                      round((tab1$hyb.npp.mod.cap.MgC/hyb.tot)*100, 0),
                      "%)", sep=""))
tab1.f <- cbind(tab1.f, 
                round(tab1$hyb.npp.MgC.ha.ground, 2),
                round(tab1$hyb.npp.MgC.ha.can, 2))

names(tab1.f) <- c("lulc", "area.ha", "med.can", "T.biom.kMgC",
                   "biom.MgC.ha.ground", "biom.MgC.ha.can", 
                   "fia.tot.npp.MgC", "fia.npp.MgC.ha.ground", "fia.npp.MgC.ha.can",
                   "hyb.tot.npp.MgC", "hyb.npp.MgC.ha.ground", "hyb.npp.MgC.ha.can")


totals <- c("ALL",
            area.tot, 
            npp.dat[aoi>800, median(can.frac, na.rm=T)],
            biom.tot,
            npp.dat[aoi>800, median(biom.MgC.ha.ground, na.rm=T)],
            npp.dat[aoi>800, median(biom.MgC.ha.can, na.rm=T)],
            fia.empir.tot,
            npp.dat[aoi>800, median(fia.empir.npp.MgC.ha.ground, na.rm=T)],
            npp.dat[aoi>800, median(fia.empir.npp.MgC.ha.can, na.rm=T)],
            hyb.tot,
            npp.dat[aoi>800, median(hyb.npp.MgC.ha.ground, na.rm=T)],
            npp.dat[aoi>800, median(hyb.npp.MgC.ha.can, na.rm=T)])
tab1.f <- as.matrix(tab1.f)
tab1.f <- rbind(tab1.f, totals)
write.csv(tab1.f, "processed/boston/results/lulc_NPP_summary.csv")



#########


### question: what is the magnitude of effect of npp on ppCO2 in local column of atmosphere?
## DOY range 135-258
gseas <- 258-135 ## number of days in growing season
hrs <- 8 ## hours per day photosynthesis works during growing season
ghrs <- hrs*gseas ## total photosynthesizing hours to distribute annual npp across
pix.npp <- 60 ## for instance, a typical is 200 kg-biomass/pix, 60 is median for FIA.forest.empir
pix.CO2 <- (pix.npp/2)*44/12 ## kg-CO2/pix
air.dens <- 1.293 ## kg/m3
air.molwt <- 28.97 # g/mol mixed air
vol.col <- 1000*30*30 ## volume of air column to 1km over 30m pixel footprint
CO2.molwt <- 44.01 ## gCO2/mol
CO2.ppm <- 400
CO2.col <- vol.col*air.dens*1000*(1/air.molwt)*CO2.ppm/1E6 ## 16k moles CO2 in column
CO2.col/400 ## 1 ppm CO2 in the column is ~40 moles of CO2

CO2.drawdown.hr <- (pix.CO2/ghrs)*1000*(1/CO2.molwt) ## avg mols-CO2/hr drawdown during growing periods

(CO2.drawdown.hr*hrs)/(CO2.col/400) 
## so a 200 kg-biomass/yr pix can draw about 1.7 ppm CO2/d out of the local column
## ## an 800 kg-biomass/yr pix can draw about 6.7 pp CO2/d (this is a 4.4 MgC/ha/yr pixels)
## max NPP is ~ 8 MgC/ha/yr 
# (8*2000)*(900/1E4) ## 1440 kg-biomass/pix
## at 1440 kg-biomass/pix, local drawdown is 12.1 ppm CO2/d
## at FIA npp.forest median (60 kg-biomass/pix), draw about 0.5 ppm CO2/d

## what about for atmos. column over entire city in aggregate?
npp.hyb <- npp.dat[aoi>800, sum(hyb.npp, na.rm=T)]
area.tot <- npp.dat[aoi>800, sum(aoi, na.rm=T)]
npp.co2 <- (npp.hyb/2)*(44/12) ## 52.7k tCO2
vol.col <- 1000*area.tot ## volume of air column to 1km over city footprint
CO2.molwt <- 44.01 ## gCO2/mol
CO2.ppm <- 400
CO2.col <- vol.col*air.dens*1000*(1/air.molwt)*(CO2.ppm/1E6) ## 2195 Mmol CO2 over city


CO2.drawdown.hr <- (npp.co2/ghrs)*1000*(1/CO2.molwt)
(CO2.drawdown.hr*hrs)  # daily avg drawdown of CO2 from 1km box over city
((CO2.drawdown.hr*hrs)/(vol.col*air.dens*1000*(1/air.molwt)))*1E6  # daily avg drawdown of CO2 from 1km box over city

### do you understand this shit about exponentials?
x <- seq(1,100)
test.b0 <- c(-1, 0, 1)
test.b1 <- c(-1, -1, -1)
plot(x, exp(test.b0[1]+(test.b1[1]*log(x))), col="blue")
for(e in 2:3){
  points(x, exp(test.b0[e]+(test.b1[e]*log(x))), col=e)
}

test.b0 <- c(0,0,0)
test.b1 <- c(-1, -0.5, -1.5)
plot(x, exp(test.b0[1]+(test.b1[1]*log(x))), col="blue")
for(e in 2:3){
  points(x, exp(test.b0[e]+(test.b1[e]*log(x))), col=e)
}


### how do the dbh growth rates seen in van Noorn 2018 translate to biomass growth rates?
dbh.test <- 25
gr.test <- 0.024
b0.test <- -2.48
b1.test <- 2.4835
biom0 <- exp(b0.test+(b1.test*log(dbh.test)))
biom1 <- exp(b0.test+(b1.test*log(dbh.test*(1+gr.test))))
(biom1-biom0)/biom0 ### this is a 3.5% biomass gain @ median 1.4%, at max 2.4% is 6.1% growth

### comparison NPP vs FFCO2 (ACES)


