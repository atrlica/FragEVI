library(raster)
library(data.table)

## the FIA use this equation for volume of wood in a stand by stand age
a=123.67
b=0.04
AGE=40
AGE=seq(0,150)
y = a*(1-exp(-b*AGE))^3
plot(AGE, y)
z = diff(y)
z.rel <- z/y[2:151]
length(y)
length(z) ## this is growth after 1 year, 2 years, etc
plot(y[2:151], z)
points(y[2:151], z.rel, pch=14, cex=0.5, col="red") #hyperbolic
abline(v=1)
abline(v=0)
plot(AGE[2:151], z) ## this is the gain curve over site age

### standing live wood C storage in local forest that resemble (?) what we'd have in Boston:
### i.e. N.red oak; red maple/oak; mixed upland hwoods; red maple uplands: Range is  94.7-105.1 MgC/ha
## read in summaries of C stock change in Q.alba/Q.rubra/Carya and Q.rubra forest
## both forest types include values for reforestation and afforestation conditions
ne.summ <- read.csv("docs/FIA_CTMANHRI.csv")
plot(ne.summ$Age.yrs, ne.summ$live.tree.c.inc)

mix.oak.ref <- read.csv("docs/FIA_QUALQURUCATO_Reforest.csv")
plot(mix.oak.ref$Age.yrs, mix.oak.ref$live.tree.c.inc)
mix.oak.aff <- read.csv("docs/FIA_QUALQURUCATO_Afforest.csv")
plot(mix.oak.aff$Age.yrs, mix.oak.aff$live.tree.c.inc)

red.oak.ref <- read.csv("docs/FIA_QURU_Reforest.csv")
plot(red.oak.ref$Age.yrs, red.oak.ref$live.tree.c.inc)
red.oak.aff <- read.csv("docs/FIA_QURU_Afforest.csv")
plot(red.oak.aff$Age.yrs, red.oak.aff$live.tree.c.inc)

fia.summ <- as.data.frame(cbind(ne.summ$Age.yrs, ne.summ$live.tree.c.inc, mix.oak.ref$live.tree.c.inc,
                                mix.oak.aff$live.tree.c.inc, red.oak.ref$live.tree.c.inc, red.oak.aff$live.tree.c.inc))
colnames(fia.summ) <- c("Age", "NE.total", "mix.oak.ref", "mix.oak.aff", "red.oak.ref", "red.oak.aff")
plot(fia.summ$Age, fia.summ$NE.total, pch=15, col="black")
points(fia.summ$Age, fia.summ$mix.oak.aff, pch=16, col="lightblue")
points(fia.summ$Age, fia.summ$mix.oak.ref, pch=16, col="blue")
points(fia.summ$Age, fia.summ$red.oak.aff, pch=17, col="pink")
points(fia.summ$Age, fia.summ$red.oak.ref, pch=17, col="red")
### the only thing that changes between reforestation and afforestation are values for forest floor and soil C

## what is relationship between standing live biomass-C and C increment
plot(ne.summ$mean.vol.m3, ne.summ$live.tree.c.inc)
plot(ne.summ$live.tree.tCha, ne.summ$live.tree.c.inc)
plot(ne.summ$Age.yrs, ne.summ$live.tree.tCha)
plot(ne.summ$Age.yrs, ne.summ$mean.vol.m3) ## basic sigmoid 0 to max at 100


#### OK with that out of the way: Let's try to figure out the tree distribution per 30m cell and uptake
biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m 
street <- read.csv("docs/ian/Boston_Street_Trees.csv") ## Street tree resurvey, biomass 2006/2014, species
# length(unique(street$Species)) ## 77 species!
street$ann.npp <- (street$biomass.2014-street$biomass.2006)/8
street <- as.data.table(street)
street[is.na(Species), Species:="Unknown unknown"]
street[Species=="Unknown", Species:="Unknown unknown"]
street[Species=="unknown", Species:="Unknown unknown"]
street[Species=="Purpleleaf plum", Species:="Prunus cerasifera"]
street[Species=="Bradford pear", Species:="Pyrus calleryana"]
street[Species=="Liriondendron tulipifera", Species:="Liriodendron tulipifera"]
street[Species=="Thornless hawthorn", Species:="Crataegus crus-galli"]
genspec <- strsplit(as.character(street$Species), " ")
gen <- unlist(lapply(genspec, "[[", 1))
street[,genus:=gen] ## 40 genera
a <- street[,length(ann.npp)/dim(street)[1], by=genus]
a <- a[order(a$V1, decreasing = T),]
sum(a$V1[1:12])
focus <- a$genus[1:12]
# focus

##########
### prototype process for using FIA aggregate data to figure npp from Raciti biomass
# 1) determine MgC/ha in 30m plot, normalize to %canopy (i.e. there is X many ha of forest there with Y much living carbon in place)
# ie. as far as FIA is concerned, that is a forest, in the area it occupies, that has Y much AGB
# 2) figure out how old FIA would say a plot like that is
# 3) Figure out the next-year tC/ha in the plot
# 4) apply this incrememnt to the area of "forest" (canopy) in the pixel

a=123.67
b=0.04
AGE=40
AGE=seq(0,120)
y = a*(1-exp(-b*AGE))^3
par(mfrow=c(2,1))
plot(AGE, y, ylab="tC/ha", xlab="Stand age, yr")
plot(AGE[1:120], diff(y))

## example: X=40 yrs,
x=40
a=123.67
b=0.04
y=a*(1-(exp(-b*x)))^3 
y ## 62.9 tC/ha
### what is equivalent age of stand if we know live biomass?
log(1-(y/a)^(1/3))/(-b) ## predicts 40 years
tC <- seq(0,120)
plot(tC, log(1-(tC/a)^(1/3))/(-b), xlab="live biomass, tC/ha", ylab="Site age, yr")

## what is age of a forest with 100% canopy and 30000 kg/pixel biomass
30000/1000/2*(10^4)/900 ### 166 MgC/ha --> out of bounds)
20000/1000/2*(10^4)/900  ### at 30m pixel values over about 20000 kg biomass, you're out of bounds for stand age -- too thick with biomass!

biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m kg-biomass to 30m pixel
aoi <- raster("processed/boston/bos.aoi30m.tif")
can <- raster("processed/boston/bos.can30m.tif")
isa <- raster("processed/boston/bos.isa.rereg30m.tif")
biom <- crop(biom, isa)
aoi <- crop(aoi, isa)
can <- crop(can, isa)
biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[,aoi:=as.vector(getValues(aoi))]
can.dat <- as.data.table(as.data.frame(can))
biom.dat[, can.frac:=can.dat$bos.can30m]
isa.dat <- as.data.table(as.data.frame(isa))
biom.dat[, isa.frac:=isa.dat$bos.isa.rereg30m]

### live MgC by area of GROUND in each pixel
biom.dat[, live.MgC.ha.ground:=(bos.biom30m/aoi)*(1/2)*(1/1000)*(10^4)] ## convert kg-biomass/pixel based on size of pixel (summed by 1m2 in "aoi"), kgC:kgbiomass, Mg:kg, m2:ha
biom.dat[aoi>800, range(live.MgC.ha.ground, na.rm=T)] ## up to 284 MgC/ha, as we'd expect from what we saw in Raciti
hist(biom.dat[aoi>800,live.MgC.ha.ground]) #v skewed, most below 50 MgC/ha
biom.MgC.ha.ground <- raster(biom)
biom.MgC.ha.ground <- setValues(biom.MgC.ha.ground, biom.dat[,live.MgC.ha.ground])
plot(biom.MgC.ha.ground)

### also possible to correct biomass figures for the amount of canopy or pervious cover present per cell
### i.e. we assume FIA is measuring trees in "forest" with essentially continuous canopy, so that differences in tC/ha as a function of age are purely a function of tree growth and not differences in tree %coverage
### live MgC/ha for "forest" fraction in each pixel
biom.dat[,live.MgC.ha.forest:=(bos.biom30m/(aoi*can.frac))*(1/2)*(1/1000)*(10^4)]
biom.dat[can.frac<=0.01, live.MgC.ha.forest:=0]
range(biom.dat[aoi>800,live.MgC.ha.forest],  na.rm=T) ## 0 - 284
hist(biom.dat[aoi>800,live.MgC.ha.forest]) ## correcting for canopy cover, more mid-rage values
biom.dat[live.MgC.ha.forest<100, length(live.MgC.ha.forest)]/dim(biom.dat[!is.na(can.frac),])[1] ## 84% of pixels are below 100 MgC/ha
biom.MgC.ha.forest <- raster(biom)
biom.MgC.ha.forest <- setValues(biom.MgC.ha.forest, biom.dat[,live.MgC.ha.forest])
plot(biom.MgC.ha.forest)

## correct for pervious cover
biom.dat[,live.MgC.ha.perv:=(bos.biom30m/(aoi*(1-isa.frac)))*(1/2)*(1/1000)*(10^4)]
biom.dat[isa.frac>0.99, live.MgC.ha.perv:=0]
range(biom.dat[aoi>800 & isa.frac<0.98,live.MgC.ha.perv],  na.rm=T) ## 0 - 6786
hist(biom.dat[aoi>800 & isa.frac<0.98,live.MgC.ha.perv]) ## a small number of very extreme values
biom.dat[live.MgC.ha.perv<100, length(live.MgC.ha.perv)]/dim(biom.dat[!is.na(can.frac),])[1] ## 75% of pixels are below 100 MgC/ha
biom.MgC.ha.perv <- raster(biom)
biom.MgC.ha.perv <- setValues(biom.MgC.ha.perv, biom.dat[,live.MgC.ha.perv])
plot(biom.MgC.ha.perv)

### get delta figures
biom.dat[, delta.C.perv:=live.MgC.ha.perv-(live.MgC.ha.ground)]
biom.dat[, delta.C.forest:=live.MgC.ha.forest-(live.MgC.ha.ground)]
plot(biom.dat[isa.frac<0.9, isa.frac], biom.dat[isa.frac<0.9, delta.C.perv])
plot(biom.dat[can.frac>0.07, can.frac], biom.dat[can.frac>0.07, delta.C.forest])


### figure out forest "age" for the cells (using coefficients for NE total)
## age based on ground area
biom.dat[,age.ground:=log(1-(live.MgC.ha.ground/a)^(1/3))/(-b)] ## some of these will produce infinities with too dense biomass
biom.dat[age.ground>120, age.ground:=120] ## fix the divergent ones to just "old, not growing"
biom.dat[!is.finite(age.ground), age.ground:=120] ## again, fix the ones that got fucked to "old, not growing"
biom.dat[is.na(aoi), age.ground:=NA] # cancel places out of bounds
biom.dat[bos.biom30m<=10, age.ground:=NA] ## cancel places with no biomass
biom.dat[is.na(bos.biom30m), age.ground:=NA]

ground.age <- raster(biom)
ground.age <- setValues(ground.age, biom.dat[,age.ground])
plot(ground.age)
hist(biom.dat[,age.ground]) ## got a lot of "old" ones, indicating high density of biomass

## age based on canopy area
biom.dat[,age.forest:=log(1-(live.MgC.ha.forest/a)^(1/3))/(-b)] ## some of these will produce infinities with too dense biomass
biom.dat[age.forest>120, age.forest:=120] ## fix the divergent ones to just "old, not growing"
biom.dat[!is.finite(age.forest), age.forest:=120] ## again, fix the ones that got fucked to "old, not growing"
biom.dat[is.na(aoi), age.forest:=NA] # cancel places out of bounds
biom.dat[bos.biom30m<=10, age.forest:=NA] ## cancel places with no biomass
biom.dat[is.na(bos.biom30m), age.forest:=NA]

forest.age <- raster(biom)
forest.age <- setValues(forest.age, biom.dat[,age.forest])
plot(forest.age)
hist(biom.dat[,age.forest]) ## many more old forest, peak has moved older

## age based on pervious area
biom.dat[,age.perv:=log(1-(live.MgC.ha.perv/a)^(1/3))/(-b)] ## some of these will produce infinities with too dense biomass
biom.dat[age.perv>120, age.perv:=120] ## fix the divergent ones to just "old, not growing"
biom.dat[!is.finite(age.perv), age.perv:=120] ## again, fix the ones that got fucked to "old, not growing"
biom.dat[is.na(aoi), age.perv:=NA] # cancel places out of bounds
biom.dat[bos.biom30m<=10, age.perv:=NA] ## cancel places with no biomass
biom.dat[is.na(bos.biom30m), age.perv:=NA]

perv.age <- raster(biom)
perv.age <- setValues(perv.age, biom.dat[,age.perv])
plot(perv.age)
hist(biom.dat[,age.perv]) ## old forests now dominant

biom.dat[,delta.age.perv:=age.perv-age.ground]
biom.dat[,delta.age.forest:=age.forest-age.ground]

## frequency distributions of different methods
par(mfrow=c(3,1), mar=c(4,4,2,1))
hist(biom.dat$age.ground,  main="Forest age, unadjusted", xlab="Age, yrs")
hist(biom.dat$age.forest,  main="Forest age, canopy area adj.", xlab="Age, yrs")
hist(biom.dat$age.perv,  main="Forest age, pervious area adj.", xlab="Age, yrs")

a=123.67
b=0.04

## figure out next annual increment possible for the "forest" average present in each cell, based on projected "age" and corrected for area
biom.dat[,npp.ann.ground:=((a*(1-(exp(-b*(age.ground+1))))^3)-(a*(1-(exp(-b*(age.ground))))^3))*(1E-4)*aoi] ## by ground area
biom.dat[,npp.ann.forest:=((a*(1-(exp(-b*(age.forest+1))))^3)-(a*(1-(exp(-b*(age.forest))))^3))*(1E-4)*(aoi*can.frac)] ## by canopy area
biom.dat[,npp.ann.perv:=((a*(1-(exp(-b*(age.perv+1))))^3)-(a*(1-(exp(-b*(age.perv))))^3))*(1E-4)*(aoi*(1-isa.frac))] ## by pervious area
biom.dat[bos.biom30m<10, npp.ann.ground:=0] ## manually set negligible biomass cells to 0
biom.dat[bos.biom30m<10, npp.ann.forest:=0]
biom.dat[bos.biom30m<10, npp.ann.forest:=0]
biom.dat[isa.frac>0.99, npp.ann.perv:=0]
biom.dat[can.frac<0.01, npp.ann.forest:=0]


### look at some plots 
plot(biom.dat$npp.ann.ground, biom.dat$npp.ann.forest)
plot(biom.dat$npp.ann.ground, biom.dat$npp.ann.perv) ## hard to say, up to 0.2 MgC/pix
par(mfrow=c(3,1))
hist(biom.dat$npp.ann.ground, main="NPP, raw area", xlab="MgC/pix/yr")
hist(biom.dat$npp.ann.forest, main="NPP, canopy area", xlab="MgC/pix/yr")
hist(biom.dat$npp.ann.perv, main="NPP, pervious area", xlab="MgC/pix/yr")

### aggregated stats
biom.dat[,sum(npp.ann.ground, na.rm=T)] #13.8k tC/yr by raw ground area
(biom.dat[,sum(npp.ann.ground, na.rm=T)]/biom.dat[,sum(aoi, na.rm=T)])*1E4 ### 1.11 MgC/ha/yr raw ground area
biom.dat[,sum(npp.ann.forest, na.rm=T)] #4.7k tC/yr by canopy area
(biom.dat[,sum(npp.ann.forest, na.rm=T)]/biom.dat[,sum(aoi, na.rm=T)])*1E4 ### 0.38 MgC/ha/yr canopy area area
biom.dat[,sum(npp.ann.perv, na.rm=T)] #5.5k tC/yr by raw ground area
(biom.dat[,sum(npp.ann.perv, na.rm=T)]/biom.dat[,sum(aoi, na.rm=T)])*1E4 ### 0.44 MgC/ha/yr pervious area
## contrast 10.3-8.9 = 1.4 NEP for boston region in Hardiman

### median age for each method
biom.dat[, median(age.ground, na.rm = T)] ##20.2
biom.dat[, median(age.forest, na.rm = T)] ##39.7
biom.dat[, median(age.perv, na.rm = T)] ##37.3


######################
### Different FIA factors to estimate 30m annual NPP (MgC/yr)

library(data.table)
library(raster)

## cleaned up code, handles different growth parameters for different forests (note "Afforestation" and "Reforestation" values are same viz. live biomass growth rates)
## initialize parameters for different forest types that ?? resemble tree species distributions in Boston
for.type <- c("NEdefault","Mixoak", "Redoak")
a <- c(123.67, 130.81, 123.09)
b <- c(0.04, 0.03, 0.04)

## test limits of the live biomass~stand age function
for(f in 1:length(for.type)){
  x=seq(0,120)
  liveC=a[f]*(1-(exp(-b[f]*x)))^3
  plot(x,liveC, main=for.type[f], ylab="live tC/ha", xlab="stand age")
  x <- seq(0, 120) ## inverse: model age vs. live biomass
  st.age=log(1-(x/a[f])^(1/3))/(-b[f]) ##
  plot(x, st.age, main=for.type[f], ylab="live tC/ha", xlab="stand age")
}
diff(st.age) ## lagged differences --> yearly increment in C gain
diff(liveC)
## conservatively, none of the models is particularly stable over 100 yrs stand age

biom.dat[, delta.npp.forest:=npp.ann.forest-npp.ann.ground]
biom.dat[, delta.npp.per:=npp.ann.perv-npp.ann.ground]

## package up some summary rasters
biom.dat.r <- biom.dat
for(g in 5:19){
  r <- raster(biom)
  r <- setValues(r, biom.dat.r[[g]])
  writeRaster(r, filename=paste("processed/boston/fia/fia", colnames(biom.dat)[g], "tif", sep="."),
              format="GTiff", overwrite=T)
}



