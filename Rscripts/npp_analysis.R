library(raster)
library(data.table)

#### all NPP will be dealt with in MgC/ha/yr, all biomass in kg-biomass/ha

## maybe look at how this sorts re. lulc
lulc <- raster("processed/boston/bos.lulc30m.lumped.tif")
# aoi <- raster("processed/boston/bos.aoi30m.tif")
# lulc <- crop(lulc, aoi)
# lulc <- mask(lulc, aoi)
# biom.dat[,lulc:=getValues(lulc)]
# hist(biom.dat[,lulc])
# 
### FIA results, packaged with some exploratory calcs
fia <- read.csv("processed/npp.FIA.v3.csv")
fia <- as.data.table(fia)
## the current npp.ann figures are in MgC/pix/yr; correct to MgC/ha/yr
fia[,npp.ann.forest.MgC.ha:=(npp.ann.forest/aoi)*1E4/(2*1000)]
fia[,npp.ann.ground.MgC.ha:=(npp.ann.ground/aoi)*1E4/(2*1000)]
fia[,npp.ann.perv.MgC.ha:=(npp.ann.perv/aoi)*1E4/(2*1000)]
names(fia) <- paste("fia", names(fia), sep=".")
# hist(fia[fia.aoi>800, fia.npp.ann.forest.MgC.ha]) ## up to ~2.0
# hist(fia[fia.aoi>800, fia.npp.ann.ground.MgC.ha]) ## up to ~2.0
# hist(fia[fia.aoi>800, fia.npp.ann.perv.MgC.ha]) ## up to ~2.0
### FIA results, empirically derived npp estiamtes
fia.empir <- read.csv("processed/npp.FIA.empirV22.csv")
fia.empir <- as.data.table(fia.empir)
names(fia.empir) <- paste("fia.empir", names(fia.empir), sep=".")
  
## andy forest results and exploration of beta range
andy.res <- read.csv("processed/andy.forest.results.csv")
# andy.betas <- load("processed/andy.forest.beta.samples") ## comes in as list beta.track
# andy.samples <- read.csv("processed/andy.forest.npp.edge.int.samples.csv")
andy.res <- as.data.table(andy.res)
names(andy.res) <- paste("andy", names(andy.res), sep=".")

## street trees
st <- read.csv("processed/streettrees.npp.simulator.v4.results.csv")
st <- as.data.table(st)
st[,med.ann.npp.MgC.ha:=med.ann.npp.all*1E-3*(1/2)/aoi*1E4]
names(st) <- paste("st", names(st), sep=".")

## pull together and clean up
npp.dat <- merge(fia, andy.res, by.x="fia.pix.ID", by.y="andy.pix.ID", all=T)
npp.dat <- merge(npp.dat, st, by.x="fia.pix.ID", by.y="st.pix.ID", all=T)
npp.dat <- merge(npp.dat, fia.empir, by.x="fia.pix.ID", by.y="fia.empir.pix.ID", all=T)
npp.dat[, fia.X:=NULL]
npp.dat[, andy.X:=NULL]
npp.dat[, st.X:=NULL]
npp.dat[, fia.empir.X:=NULL]
# sum(npp.dat$fia.aoi-npp.dat$andy.aoi, na.rm=T)
# sum(npp.dat$fia.aoi-npp.dat$st.aoi, na.rm=T)
npp.dat[,andy.aoi:=NULL]
npp.dat[,st.aoi:=NULL]
# sum(npp.dat$fia.can.frac-npp.dat$andy.can, na.rm=T)
npp.dat[,andy.can:=NULL]
# sum(npp.dat$fia.can.frac-npp.dat$st.bos.can30m, na.rm=T)
npp.dat[,st.bos.can30m:=NULL]
# sum(npp.dat$fia.isa.frac-npp.dat$andy.isa, na.rm=T)
# mean(npp.dat$fia.isa.frac-npp.dat$andy.isa, na.rm=T) ## on avg they are extremely similar but some do deviate a lot
# plot(npp.dat$fia.isa.frac, npp.dat$andy.isa) ## andy dat doesn't need ISA, but is using the un-re-registered one, FIA uses the re-registered
npp.dat[,andy.isa:=NULL]
# sum(npp.dat$fia.bos.biom30m-npp.dat$andy.biom, na.rm=T)
npp.dat[,andy.biom:=NULL]
# sum(npp.dat$fia.bos.biom30m-npp.dat$st.bos.biom30m, na.rm=T)
npp.dat[,st.bos.biom30m:=NULL]
# sum(npp.dat$fia.bos.biom30m-npp.dat$st.biom.kg, na.rm=T)
npp.dat[,st.biom.kg:=NULL]
# plot(npp.dat$fia.bos.biom30m, npp.dat$st.med.biom.all)
# abline(a=0, b=1)
npp.dat[,fia.empir.bos.biom30m:=NULL]
npp.dat[,fia.empir.aoi:=NULL]
npp.dat[,fia.empir.can.frac:=NULL]
npp.dat[,fia.empir.isa.frac:=NULL]
npp.dat <- cbind(npp.dat, getValues(lulc))
names(npp.dat)[c(1:5, dim(npp.dat)[2])] <- c("pix.ID", "bos.biom30m", "aoi", "can.frac", "isa.frac", "lulc")

## why do this? fuck it, tiny biomass is going to come out with tiny values of npp anyway
# ### set all npp estimates to 0 in tiny biomass cells
# npp.dat[bos.biom30m<10 & is.finite(aoi), c(16:23):=0, with=F] ## FIA npp and derivatives thereof
# npp.dat[bos.biom30m<10 & is.finite(aoi), c(30,31):=0, with=F] ## andy total npps
# npp.dat[andy.ed.biom<10 & is.finite(aoi), andy.npp.edge:=0] ## andy edge fraction npp
# npp.dat[andy.int.biom<10 & is.finite(aoi), andy.npp.int:=0] ## andy interior fraction npp
# npp.dat[bos.biom30m<10 & is.finite(aoi), st.med.ann.npp.all:=0] ## median street npp estimate
# npp.dat[bos.biom30m<10 & is.finite(aoi), st.med.ann.npp.MgC.ha:=0] ## median street npp estimate MgC


### package up series of rasters for visual exploratory
# r <- lulc
# for(f in grep(names(npp.dat), pattern="fia.")){
#   a <- setValues(r, npp.dat[[f]])
#   writeRaster(a, filename=paste("processed/boston/results/fia/", names(npp.dat)[f], ".tif", sep=""), overwrite=T, format="GTiff")
#   print(paste("exported", names(npp.dat)[f]))
# }
#### FIA: raw area gives an overall younger forest, canopy basis is in general older than raw ground, perv gives some older estimates, a lot of non-retreives in perv
#### more extremes in biomass density with raw ground; canopy has a more even spread, and perv has some hot spots in the v. built up parts of town
## canopy basis tends to raise density broadly, and pervious raises density only in high biomass high pervious places (residential nhoods)
## npp.ground is maxed out almost everywhere except the highest density places, npp.forest is maxed out only in places with nearly continuous canopy an dmiddling otherwise, 
## general reduction, some severe, in npp with either forest or perv, forest reduces more in the deep canopy parts while perv leaves it untouched

# r <- lulc
# for(f in grep(names(npp.dat), pattern="andy.")){
#   a <- setValues(r, npp.dat[[f]])
#   writeRaster(a, filename=paste("processed/boston/results/andy/", names(npp.dat)[f], ".tif", sep=""), overwrite=T, format="GTiff")
#   print(paste("exported", names(npp.dat)[f]))
# }
### andy plots: The only places with meaningful interior biomass are the forest fragments, similar canopy fraction
## npp follows some reasonably predictable association with biomass
## some really dipshit high values for npp MgC/ha
# hist(npp.dat$andy.npp.tot)
# hist(npp.dat$andy.npp.tot.MgC.ha, breaks=1000, xlim=c(0,10)) ## realistically everything peters out before 8 tC/ha/yr
# npp.dat[andy.npp.tot.MgC.ha>10,] ## the insane values are all partial pixels
# hist(npp.dat[aoi>800, andy.npp.tot.MgC.ha], breaks=30) ## voila
# hist(npp.dat[aoi>800, fia.npp.ann.forest.MgC.ha], breaks=30) ## contrast: stops at ~2.0

# ### the street tree results
# r <- lulc
# for(f in grep(names(npp.dat), pattern="st.")){
#   a <- setValues(r, npp.dat[[f]])
#   writeRaster(a, filename=paste("processed/boston/results/street/", names(npp.dat)[f], ".tif", sep=""), overwrite=T, format="GTiff")
#   print(paste("exported", names(npp.dat)[f]))
# }
## process failures are mostly in the forest fragments, but not exclusively
## tree number follows the biomass contours ok
## basically every pixel that successfully runs gets all 100 samples done
## median dbh for the median npp set is pretty even across the city -- close to median for all street trees
### median npp is higher in the forest fragments (but these are full of no-sim holes), but there's some anomalies of high npp like the fenway, river valleys, whole southern part is hot

### how many successful estimates do we have for each approach?
npp.dat[is.finite(st.med.ann.npp.MgC.ha) & aoi>800,length(st.med.ann.npp.MgC.ha)] # 103850 complete pixel street tree retrievals
npp.dat[bos.biom30m<10, st.med.ann.npp.all:=0] ## manually fix the negligible biomass retreivals did not simulate to begin with
npp.dat[bos.biom30m<10, st.med.ann.npp.MgC.ha:=0]
npp.dat[is.finite(st.med.ann.npp.MgC.ha) & aoi>800,length(st.med.ann.npp.MgC.ha)] # up to 132896 retrievals
npp.dat[is.finite(fia.npp.ann.forest.MgC.ha) & aoi>800, length(fia.npp.ann.forest.MgC.ha)] #135705 c.p. fia retrievals
npp.dat[is.finite(andy.npp.tot.MgC.ha) & aoi>800, length(andy.npp.tot.MgC.ha)] ## 135705 c.p. andy retreivals
npp.dat[is.finite(fia.empir.npp.kg.hw.forest) & aoi>800, length(fia.empir.live.kgbiom.ha.forest)] ## 135705
npp.dat[!is.na(bos.biom30m) & aoi>800, length(bos.biom30m)] ## 135705 biomass records in relatively complete pixels
## get more retreivals in the forests for the street tree sim than for FIA (can't retrieve a good age using the equation approach)


## comparison: how much npp do we see by these different methods
hist(npp.dat[aoi>800, st.med.ann.npp.MgC.ha]) ## 0-5 MgC/ha/yr
hist(npp.dat[aoi>800, fia.npp.ann.forest.MgC.ha]) ## 0-2.5 MgC/ha/yr
hist(npp.dat[aoi>800, andy.npp.tot.MgC.ha]) ## 0-8 MgC/ha/yr
hist(npp.dat[aoi>800, fia.empir.npp.kg.hw.forest*(1E4/aoi)/2000]) ## up to~5 MgC/ha/yr

sum(npp.dat[, st.med.ann.npp.MgC.ha*(aoi/1E4)], na.rm=T) ## 13.2k tC/yr street tree sim (missing a bunch of high-value pixels)
sum(npp.dat[, fia.npp.ann.forest.MgC.ha*(aoi/1E4)], na.rm=T) ## 4.6k tC/yr ## fia forest method
sum(npp.dat[, andy.npp.tot.MgC.ha*(aoi/1E4)], na.rm=T) ## 13.8k tC/yr ## andy trees, avg. dbh
sum(npp.dat[, andy.npp.tot.ps.MgC.ha*(aoi/1E4)], na.rm=T) ## 12.7k tC/yr ## andy trees with pseudoreps
sum(npp.dat[, fia.empir.npp.kg.hw.forest/2000], na.rm=T) ## 8.7k tC/yr for fia empirical(forest) method


## visualize the differences in npp retreivals between them (in actual MgC per pixel)
## note that FIA npp is in MgC, not in kg biomass
### FIA, forest (canopy) areal basis
npp.dat[,diff.fia.forest.andy:=fia.npp.ann.forest-andy.npp.tot] ## FIA(forest) vs. andy
npp.dat[,diff.fia.forest.st:=fia.npp.ann.forest-st.med.ann.npp.all] ## FIA(forest) vs. street trees
npp.dat[,diff.andy.st:=andy.npp.tot-st.med.ann.npp.all]### street trees vs. andy

## FIA, ground basis
npp.dat[,diff.fia.ground.andy:=fia.npp.ann.ground-andy.npp.tot] ## FIA(ground) vs. andy
npp.dat[,diff.fia.ground.st:=fia.npp.ann.ground-st.med.ann.npp.all] ## FIA(ground) vs. street trees

#### what do these look like 
range(npp.dat[,fia.npp.ann.forest],  na.rm=T) ## up to 394 kg-biom/pix
range(npp.dat[,fia.npp.ann.ground],  na.rm=T) ## up to 395 kg-biom/pix
range(npp.dat[,andy.npp.tot], na.rm=T)## up to 1443 kg-biom/pix
range(npp.dat[,andy.npp.tot.ps], na.rm=T)## up to 1405 kg-biom/pix
range(npp.dat[,(st.med.ann.npp.all)], na.rm=T) ## up to 936 kg-biom/pix

### differences between estimate sets (downside is always bigger for fia than upside)
range(npp.dat[,diff.fia.forest.andy], na.rm=T) 
range(npp.dat[,diff.fia.forest.st], na.rm=T)
range(npp.dat[,diff.andy.st], na.rm=T)
range(npp.dat[,diff.fia.ground.andy], na.rm=T) 
range(npp.dat[,diff.fia.ground.st], na.rm=T) 

## export the rasters
r <- lulc
for(f in grep(names(npp.dat), pattern="diff.")){
  a <- setValues(r, npp.dat[[f]])
  writeRaster(a, filename=paste("processed/boston/results/", names(npp.dat)[f], ".tif", sep=""), overwrite=T, format="GTiff")
  print(paste("exported", names(npp.dat)[f]))
}

## fia.GROUND vs. street: Street estimates more in forest patches, but FIA tends to  be higher in a lot of the residential where it looks like young forest (both probably missing a fair amount of npp in the forest fragments neither can retrieve reliably)
### but forest fragments are also where the street tree set is struggling and sampling larger trees (iffy retreival quality)
## fia.GROUND vs. andy: hard to say, differs all over. Forest edges come in higher in andy, but most low/moderate biomass pixels come in higher in FIA so is it low-biomass common pix that drive npp, or high-biomass rare pix that drive npp

## andy vs. street: street predicts a little bit more in residential consistently, a LOT more in forest interiors, and less on forest edges

## fia.FOREST vs. andy: andy is generally a tiny bit higher except for leafier parts where it is a lot higher; the only parts with FIA higher are non-forest but still moderate biomass (where it looks like a young forest)
## fia.FOREST vs. street: basically same as andy forest, higher all over esp. in forest except for a few little spots with moderate biomass

##### HYBRID urban productivity map
### let us then combine the andy and street maps into a single "let's get real" map
## first where are the difficult sim pixels re. LULC
npp.dat[st.sim.incomp==1, length(pix.ID), by=lulc]
## bulk are forest (2k), but significant numbers in dev and hdres
npp.dat[st.sim.incomp==1 & aoi>800, mean(bos.biom30m, na.rm=T), by=lulc]
## failures are all hella high biomass except lowveg and dev
## how do andy FOREST compare to street FOREST
plot(npp.dat[lulc==1 & aoi>800, andy.npp.tot], npp.dat[lulc==1 & aoi>800, st.med.ann.npp.all])
abline(a=0, b=1) ## street trees trend higher except above about 800 kg/yr npp
## big kink at high biomass where the samples get weighted up
plot(npp.dat[lulc==1 & aoi>800, bos.biom30m], npp.dat[lulc==1 & aoi>800, andy.npp.tot])
## can see the different slopes
plot(npp.dat[aoi>800, bos.biom30m], npp.dat[aoi>800, andy.npp.tot])
## same across all pixels
plot(npp.dat[lulc==1 & aoi>800, bos.biom30m], npp.dat[lulc==1 & aoi>800, st.med.ann.npp.all])
## classic kink

### let's make a hybrid
### version using growth in andy avg.dbh
npp.dat[,hyb.npp:=andy.npp.tot]
npp.dat[lulc!=1, hyb.npp:=st.med.ann.npp.all]
npp.dat[aoi>800 & is.finite(hyb.npp), length(hyb.npp)] #134643
plot(npp.dat[aoi>800, bos.biom30m], npp.dat[aoi>800, hyb.npp], 
     col=npp.dat[aoi>800, lulc], main="avg.dbh")
# ### still get that kink in the street trees

npp.dat[,hyb.npp.ps:=andy.npp.tot.ps]
npp.dat[lulc!=1, hyb.npp.ps:=st.med.ann.npp.all]
npp.dat[aoi>800 & is.finite(hyb.npp.ps), length(hyb.npp.ps)] #134643
plot(npp.dat[aoi>800, bos.biom30m], npp.dat[aoi>800, hyb.npp.ps], 
     col=npp.dat[aoi>800, lulc], main="with pseudoreps")

# ### how can we reduce the artifacts in the street tree sim results, especially over 20k kg pixels?
# plot(npp.dat[aoi>800, bos.biom30m], npp.dat[aoi>800, st.med.ann.npp.all], col=npp.dat[aoi>800, st.max.wts])
# plot(npp.dat$bos.biom30m, npp.dat$st.max.wts) ## right at >20000 the shit starts to get super weighted
# 
# breaks <- seq(0, 60000, by=2000)
# # specify interval/bin labels
# labels <- paste("<", as.character(seq(2000, 60000, by=2000)), sep="")
# # bucketing data points into bins
# bins <- cut(npp.dat[aoi>800, bos.biom30m], breaks, include.lowest = T, right=FALSE, labels=labels)
# # inspect bins
# summary(bins)
# duur <- npp.dat[aoi>800,]
# duur$biom.bin <- bins
# duur[,median(st.max.wts, na.rm=T), by=biom.bin] ### so really at >20000 the simulator really begins to struggle
# plot(duur$bos.biom30m, duur$st.max.wts)
# duur[st.max.wts>50,length(st.max.wts)] ## 20k pix show signs of struggling to fit
# f <- duur[st.max.wts>50,length(st.max.wts), by=biom.bin] #almost half are in the ultra-low biomass category
# g <- duur[, length(st.max.wts), by=biom.bin]
# h <- merge(f, g, by="biom.bin", all=T)
# h$prop <- h$V1.x/h$V1.y ## nothing over 20% high weighting until >20000
# duur[bos.biom30m>20000,] ## 8041 pix over 20k -- let's just swap the fuckers out
# 
npp.dat[aoi>800 & bos.biom30m>20000, length(bos.biom30m)]/npp.dat[aoi>800, length(bos.biom30m)] ## 6% of pixels are over 20k
### swap the >20000kg biomass pixels with the andy forest npp estimates irrespective of lulc class; 20k/pix is 111 MgC/ha, about what a mature forest would be
npp.dat[bos.biom30m>20000, hyb.npp:=andy.npp.tot]
npp.dat[bos.biom30m>20000, hyb.npp.ps:=andy.npp.tot.ps]

### export the collated results
# write.csv(npp.dat, "processed/npp.estimates.V1.csv")
write.csv(npp.dat, "processed/npp.estimates.V2.csv") ## with pseudorep based andy trees

# npp.dat <- read.csv("processed/npp.estimates.V1.csv")
npp.dat <- read.csv("processed/npp.estimates.V2.csv")
npp.dat <- as.data.table(npp.dat)

## contrast (avg.dbh)
plot(npp.dat[aoi>800, bos.biom30m], npp.dat[aoi>800, hyb.npp], col=as.numeric(npp.dat[aoi>800,lulc]))
points(npp.dat[aoi>800, bos.biom30m], npp.dat[aoi>800, fia.npp.ann.ground], col="goldenrod", pch=5)
points(npp.dat[aoi>800, bos.biom30m], npp.dat[aoi>800, fia.npp.ann.forest], col="lightgreen", pch=7)
points(npp.dat[aoi>800, bos.biom30m], npp.dat[aoi>800, fia.npp.ann.perv], col="gray55", pch=9)
points(npp.dat[aoi>800, bos.biom30m], npp.dat[aoi>800, fia.empir.npp.kg.hw.forest], col="lightblue", pch=11)
## all the fia estimates are forced into a low hump at the bottom end of the scale, below 500 kg-biomass/yr
## fia empirical estimates are a lower arc below the hybrid output

### make a proper plot
par(mfrow=c(1,1), mar=c(4,4,1,1))
col.lulc <- c("forestgreen", "grey55", "salmon", "gold2", "lawngreen", "cadetblue3")
plot(npp.dat[aoi>800, (bos.biom30m/aoi)*1E4/2000], npp.dat[aoi>800, (hyb.npp/aoi)*1E4/2000], 
     col=col.lulc[npp.dat[aoi>800, lulc]], pch=14, cex=0.1, 
     xlab="biomass MgC/ha", ylab="NPP MgC/ha/yr")
legend(fill=col.lulc, legend=c("Forest", "Developed", "HDResid", "LDResid", "Other Veg", "Water"), x=200, y=4)
abline(v=111, lwd=2, lty=2)
# points(npp.dat[aoi>800, (bos.biom30m/aoi)*1E4/2000], npp.dat[aoi>800, (st.med.ann.npp.all/aoi)*1E4/2000],
#        col="purple", cex=0.2)


## contrast overlay the fia results
par(mfrow=c(1,1), mar=c(4,4,1,1))
col.lulc <- c("forestgreen", "grey55", "salmon", "gold2", "lawngreen", "cadetblue3")
plot(npp.dat[aoi>800, (bos.biom30m/aoi)*1E4/2000], npp.dat[aoi>800, (hyb.npp/aoi)*1E4/2000], 
     col="grey65", cex=0.1, 
     xlab="biomass MgC/ha", ylab="NPP MgC/ha/yr")
points(npp.dat[aoi>800, (bos.biom30m/aoi)*1E4/2000], npp.dat[aoi>800, (fia.npp.ann.forest/aoi)*1E4/2000],
       pch=14, col="forestgreen", cex=0.3)
points(npp.dat[aoi>800, (bos.biom30m/aoi)*1E4/2000], npp.dat[aoi>800, (fia.npp.ann.ground/aoi)*1E4/2000],
       pch=14, col="red", cex=0.3)
legend(fill=c("red", "forestgreen"), legend=c("FIA.ground", "FIA.canopy"), x=200, y=4)
abline(v=123, lwd=2, lty=4)


## how do it compare?
options(scipen=999)
npp.dat[aoi>800, .((sum(hyb.npp, na.rm=T)/1000),
                   (sum(fia.npp.ann.ground, na.rm=T)/1000),
                   (sum(fia.npp.ann.forest, na.rm=T)/1000), 
                   (sum(fia.npp.ann.perv, na.rm=T)/1000),
                   (sum(fia.empir.npp.kg.hw.forest, na.rm=T)/1000)), by=lulc] ## Mg-biomass per LULC
npp.dat[aoi>800, .((sum(hyb.npp, na.rm=T)/1000), 
                   sum(andy.npp.tot, na.rm=T)/1000, 
                   sum(fia.npp.ann.ground, na.rm=T)/1000, 
                   sum(fia.npp.ann.forest, na.rm=T)/1000, 
                   sum(fia.npp.ann.perv, na.rm=T)/1000,
                   sum(fia.empir.npp.kg.hw.forest, na.rm=T)/1000),] ## Mg-biomass total aoi
## total is close between hybrid and ground, but much higher than forest or perv; fia.empir is middle of extreme for fia.equation
## hybrid vs. ground (otherwise beats)
## lower in dev hdres lowveg, but beats by a lot in forest (i.e. the most comparable parts)


### summary table by lulc
tot.area <- npp.dat[aoi>800 & !is.na(lulc), sum(aoi)]
fox <- npp.dat[aoi>800 & !is.na(lulc), 
               .((sum(aoi)/tot.area), 
                 (median(can.frac, na.rm=T)), 
                 (sum(bos.biom30m, na.rm=T)/1000), 
                 (sum(bos.biom30m, na.rm=T)/1000)/(sum(aoi)/1E4)), by=lulc]
names(fox) <- c("lulc", "area.frac", "med.can.frac", "tot.biomass.Mg", "biomass.Mg.ha")

## productivity summaries
hyb.tot <- npp.dat[aoi>800 & !is.na(lulc), sum(hyb.npp, na.rm=T)/(1000*2)]
hen <- npp.dat[aoi>800 & !is.na(lulc),
               .((sum(hyb.npp, na.rm=T)/2000),
                 (sum(hyb.npp, na.rm=T)/2000)/(sum(aoi, na.rm=T)/1E4),
                 (sum(hyb.npp, na.rm=T)/2000)/hyb.tot), by=lulc]
names(hen) <- c("lulc", "hyb.npp.MgC", "hyb.npp.MgC.ha", "hyb.npp.frac.tot")
house <- merge(fox, hen, by="lulc")

flop.tot <- npp.dat[aoi>800 & !is.na(lulc), sum(fia.npp.ann.forest, na.rm=T)/(1000*2)]
flop <- npp.dat[aoi>800 & !is.na(lulc), 
                .((sum(fia.npp.ann.forest, na.rm=T)/2000),
                  (sum(fia.npp.ann.forest, na.rm=T)/2000)/(sum(aoi, na.rm=T)/1E4),
                  (sum(fia.npp.ann.forest, na.rm=T)/2000)/flop.tot), by=lulc]
names(flop) <- c("lulc", "fia.npp.MgC", "fia.npp.MgC.ha", "fia.npp.frac.tot")
house <- merge(house, flop, by="lulc")
write.csv(house, "processed/boston/results/summary.hybrid.fia.npp.csv")

## productivity across all lulc
npp.dat[aoi>800, (sum(fia.npp.ann.ground, na.rm=T)/(1000*2))/sum(aoi, na.rm=T)]*1E4 ### 1.11 MgC/ha
npp.dat[aoi>800, (sum(fia.npp.ann.forest, na.rm=T)/(1000*2))/sum(aoi, na.rm=T)]*1E4 ### 0.37 MgC/ha
npp.dat[aoi>800, (sum(fia.npp.ann.perv, na.rm=T)/(1000*2))/sum(aoi, na.rm=T)]*1E4 ### 0.43 MgC/ha
npp.dat[aoi>800, (sum(fia.empir.npp.kg.hw.forest, na.rm=T)/(1000*2))/sum(aoi, na.rm=T)]*1E4 ### 0.7 MgC/ha
npp.dat[aoi>800, (sum(andy.npp.tot, na.rm=T)/(1000*2))/sum(aoi, na.rm=T)]*1E4 ### 1.11 MgC/ha
npp.dat[aoi>800, (sum(st.med.ann.npp.all, na.rm=T)/(1000*2))/sum(aoi, na.rm=T)]*1E4 ### 1.07 MgC/ha (excluding non-retrieves)
npp.dat[aoi>800, (sum(hyb.npp, na.rm=T)/(1000*2))/sum(aoi, na.rm=T)]*1E4 ### 1.17 MgC/ha
npp.dat[aoi>800, (sum(hyb.npp, na.rm=T)/(1000*2))] ### 14,393



### only thing is to decide whether or not to swap the high-biomass non-forest sim results with the andy forest equivalents
## I thiknk this is a yes: Anything >20000 kg/pix is a "forest" in the andy sense

#### BRING IN NPP estimates based on the FIA empirical records
npp.dat <- as.data.table(read.csv("processed/npp.estimates.V1.csv"))
npp.tab <- npp.dat[aoi>800 & !is.na(lulc), .(round(sum(fia.empir.npp.kg.hw.ground, na.rm=T)/(1000*2), digits=0),
                              round(sum(fia.empir.npp.kg.hw.forest, na.rm=T)/(1000*2), digits=0),
                              round(sum(fia.empir.npp.kg.hw.perv, na.rm=T)/(1000*2), digits=0),
                              round(sum(hyb.npp, na.rm=T)/(1000*2), digits=0)), by=lulc]
npp.tab$lulc.name <- c("Dev", "Other Veg", "Forest", "Water", "HD Resid", "LD Resid")
names(npp.tab) <- c("lulc", "FIA (Ground)", "FIA (canopy)", "FIA (pervious)", "Street+Urban Forest", "Type")
npp.tab <- cbind(npp.tab[,2:5, with=F], npp.tab$Type)
npp.tot <- apply(npp.tab[,1:4, with=F], 2, sum)
npp.tot <- matrix(nrow=1, ncol=5, c(npp.tot, "Total"))
npp.tab <- as.matrix(npp.tab)
npp.tab <- rbind(npp.tab, npp.tot)
npp.tab <- as.data.frame(cbind(npp.tab[,5], npp.tab[,1:4]))
names(npp.tab)[1] <- "Type"

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


### comparison NPP vs FFCO2 (ACES)


