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
  
## andy forest results and exploration of beta range
andy.res <- read.csv("processed/andy.forest.results.csv")
andy.betas <- load("processed/andy.forest.beta.samples") ## comes in as list beta.track
andy.samples <- read.csv("processed/andy.forest.npp.edge.int.samples.csv")
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
npp.dat[, fia.X:=NULL]
npp.dat[, andy.X:=NULL]
npp.dat[, st.X:=NULL]
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
r <- lulc
for(f in grep(names(npp.dat), pattern="fia.")){
  a <- setValues(r, npp.dat[[f]])
  writeRaster(a, filename=paste("processed/boston/results/fia/", names(npp.dat)[f], ".tif", sep=""), overwrite=T, format="GTiff")
  print(paste("exported", names(npp.dat)[f]))
}
#### FIA: raw area gives an overall younger forest, canopy basis is in general older than raw ground, perv gives some older estimates, a lot of non-retreives in perv
#### more extremes in biomass density with raw ground; canopy has a more even spread, and perv has some hot spots in the v. built up parts of town
## canopy basis tends to raise density broadly, and pervious raises density only in high biomass high pervious places (residential nhoods)
## npp.ground is maxed out almost everywhere except the highest density places, npp.forest is maxed out only in places with nearly continuous canopy an dmiddling otherwise, 
## general reduction, some severe, in npp with either forest or perv, forest reduces more in the deep canopy parts while perv leaves it untouched

r <- lulc
for(f in grep(names(npp.dat), pattern="andy.")){
  a <- setValues(r, npp.dat[[f]])
  writeRaster(a, filename=paste("processed/boston/results/andy/", names(npp.dat)[f], ".tif", sep=""), overwrite=T, format="GTiff")
  print(paste("exported", names(npp.dat)[f]))
}
### andy plots: The only places with meaningful interior biomass are the forest fragments, similar canopy fraction
## npp follows some reasonably predictable association with biomass
## some really dipshit high values for npp MgC/ha
# hist(npp.dat$andy.npp.tot)
# hist(npp.dat$andy.npp.tot.MgC.ha, breaks=1000, xlim=c(0,10)) ## realistically everything peters out before 8 tC/ha/yr
# npp.dat[andy.npp.tot.MgC.ha>10,] ## the insane values are all partial pixels
# hist(npp.dat[aoi>800, andy.npp.tot.MgC.ha], breaks=30) ## voila
# hist(npp.dat[aoi>800, fia.npp.ann.forest.MgC.ha], breaks=30) ## contrast: stops at ~2.0

### the street tree results
r <- lulc
for(f in grep(names(npp.dat), pattern="st.")){
  a <- setValues(r, npp.dat[[f]])
  writeRaster(a, filename=paste("processed/boston/results/street/", names(npp.dat)[f], ".tif", sep=""), overwrite=T, format="GTiff")
  print(paste("exported", names(npp.dat)[f]))
}
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
npp.dat[!is.na(bos.biom30m) & aoi>800, length(bos.biom30m)] ## 135705 biomass records in relatively complete pixels
## get more retreivals in the forests for the street tree sim than for FIA (can't retrieve a good age using the equation approach)


## comparison: how much npp do we see by these different methods
hist(npp.dat[aoi>800, st.med.ann.npp.MgC.ha]) ## 0-5 MgC/ha/yr
hist(npp.dat[aoi>800, fia.npp.ann.forest.MgC.ha]) ## 0-2.5 MgC/ha/yr
hist(npp.dat[aoi>800, andy.npp.tot.MgC.ha]) ## 0-8 MgC/ha/yr

sum(npp.dat[, st.med.ann.npp.MgC.ha*(aoi/1E4)], na.rm=T) ## 13.2k tC/yr street tree sim (missing a bunch of high-value pixels)
sum(npp.dat[, fia.npp.ann.forest.MgC.ha*(aoi/1E4)], na.rm=T) ## 4.6k tC/yr ## fia forest method
sum(npp.dat[, andy.npp.tot.MgC.ha*(aoi/1E4)], na.rm=T) ## 13.8k tC/yr ## andy trees

## visualize the differences in npp retreivals between them (in actual MgC per pixel)
## note that FIA npp is in MgC, not in kg biomass
### FIA, forest (canopy) areal basis
npp.dat[,diff.fia.forest.andy:=fia.npp.ann.forest-andy.npp.tot] ## FIA(forest) vs. andy
npp.dat[,diff.fia.forest.st:=fia.npp.ann.forest-st.med.ann.npp.all] ## FIA(forest) vs. street trees
npp.dat[,diff.andy.st:=andy.npp.tot-st.med.ann.npp.all]### street trees vs. andy

## FIA, ground basis
npp.dat[,diff.fia.ground.andy:=fia.npp.ann.ground-andy.npp.tot] ## FIA(ground) vs. andy
npp.dat[,diff.fia.ground.st:=fia.npp.ann.ground-st.med.ann.npp.all] ## FIA(ground) vs. street trees

#### what do these look like (figures quoted below are fucked -- reran FIA to npp in kg-biomass)
range(npp.dat[,fia.npp.ann.forest],  na.rm=T) ## up to 0.2 MgC/pix
range(npp.dat[,fia.npp.ann.ground],  na.rm=T) ## up to 0.2 MgC/pix
range(npp.dat[,(andy.npp.tot.MgC.ha*(aoi/1E4))], na.rm=T) ## up to 0.72
range(npp.dat[,(st.med.ann.npp.MgC.ha*(aoi/1E4))], na.rm=T) ## up to 0.47

range(npp.dat[,diff.fia.forest.andy], na.rm=T) ## -0.44 to 0.13
range(npp.dat[,diff.fia.forest.st], na.rm=T) ## -0.46 to 0.09
range(npp.dat[,diff.andy.st], na.rm=T) ## -0.16 to 0.38
range(npp.dat[,diff.fia.ground.andy], na.rm=T) ## -0.48 to 0.13
range(npp.dat[,diff.fia.forest.st], na.rm=T) ## -0.46 to 0.09

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
npp.dat[,hyb.npp:=andy.npp.tot]
npp.dat[lulc!=1, hyb.npp:=st.med.ann.npp.all]
npp.dat[aoi>800 & is.finite(hyb.npp), length(hyb.npp)] #134643
plot(npp.dat[aoi>800, bos.biom30m], npp.dat[aoi>800, hyb.npp], col=npp.dat[aoi>800, lulc])
# ### still get that kink in the street trees
# ## could just elminate the street sims above a certain biomass
# plot(npp.dat[aoi>800, bos.biom30m], npp.dat[aoi>800, st.med.ann.npp.all], col=npp.dat[aoi>800, st.max.wts])
# plot(npp.dat$bos.biom30m, npp.dat$st.max.wts) ## right at >20000 the shit starts to get super weighted

## contrast
points(npp.dat[aoi>800, bos.biom30m], npp.dat[aoi>800, fia.npp.ann.ground], col="goldenrod", pch=5)
points(npp.dat[aoi>800, bos.biom30m], npp.dat[aoi>800, fia.npp.ann.forest], col="lightgreen", pch=7)
points(npp.dat[aoi>800, bos.biom30m], npp.dat[aoi>800, fia.npp.ann.perv], col="gray55", pch=9)
## all the fia estimates are forced into a low hump at the bottom end of the scale, below 500 kg-biomass/yr

## how do it compare?
options(scipen=999)
npp.dat[aoi>800, sum(hyb.npp, na.rm=T)/1000, by=lulc] ## Mg-biomass
npp.dat[aoi>800, sum(fia.npp.ann.ground, na.rm=T)/1000, by=lulc] ## Mg-Biomass
npp.dat[aoi>800, .((sum(hyb.npp, na.rm=T)/1000),(sum(fia.npp.ann.ground, na.rm=T)/1000),(sum(fia.npp.ann.forest, na.rm=T)/1000), (sum(fia.npp.ann.perv, na.rm=T)/1000)), by=lulc] ## Mg-biomass
npp.dat[aoi>800, .((sum(hyb.npp, na.rm=T)/1000), sum(andy.npp.tot, na.rm=T)/1000, sum(fia.npp.ann.ground, na.rm=T)/1000, sum(fia.npp.ann.forest, na.rm=T)/1000, sum(fia.npp.ann.perv, na.rm=T)/1000),]
## total is close between hybrid and ground, but much higher than forest or perv
## hybrid vs. ground (otherwise beats)
## lower in dev hdres lowveg, but beats by a lot in forest (i.e. the most comparable parts)

### what is productivity rate in different LULC?
npp.dat[aoi>800 & !is.na(lulc), (sum(hyb.npp, na.rm=T)/(1000*2))/sum(aoi, na.rm=T), by=lulc]*1E4 ### GUAO, up to 3.5 tC/ha in forest, 1.4 tC/ha residential
npp.dat[aoi>800 & !is.na(lulc), median(can.frac, na.rm=T), by=lulc] ### near total coverage in forest, 60% LD, 37%HD, <2% dev, 14% lowveg
npp.dat[aoi>800 & !is.na(lulc), (sum(fia.npp.ann.forest, na.rm=T)/(1000*2))/sum(aoi, na.rm=T), by=lulc]*1E4 ### 

## productivity across all lulc
npp.dat[aoi>800, (sum(fia.npp.ann.ground, na.rm=T)/(1000*2))/sum(aoi, na.rm=T)]*1E4 ### 1.11 MgC/ha
npp.dat[aoi>800, (sum(fia.npp.ann.forest, na.rm=T)/(1000*2))/sum(aoi, na.rm=T)]*1E4 ### 0.37 MgC/ha
npp.dat[aoi>800, (sum(fia.npp.ann.perv, na.rm=T)/(1000*2))/sum(aoi, na.rm=T)]*1E4 ### 0.43 MgC/ha
npp.dat[aoi>800, (sum(andy.npp.tot, na.rm=T)/(1000*2))/sum(aoi, na.rm=T)]*1E4 ### 1.11 MgC/ha
npp.dat[aoi>800, (sum(st.med.ann.npp.all, na.rm=T)/(1000*2))/sum(aoi, na.rm=T)]*1E4 ### 1.07 MgC/ha (excluding non-retrieves)
npp.dat[aoi>800, (sum(hyb.npp, na.rm=T)/(1000*2))/sum(aoi, na.rm=T)]*1E4 ### 1.13 MgC/ha

### only thing is to decide whether or not to swap the high-biomass non-forest sim results with the andy forest equivalents





