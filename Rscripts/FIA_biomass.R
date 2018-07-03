library(raster)
library(data.table)

## the FIA use this equation for volume of wood in a stand by stand age
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


par(mfrow=c(1,1))
### live MgC by area of GROUND in each pixel
biom.dat[, live.MgC.ha.ground:=(bos.biom30m/aoi)*(1/2)*(1E-03)*(1E4)] ## convert kg-biomass/pixel based on size of pixel (summed by 1m2 in "aoi"), kgC:kgbiomass, Mg:kg, m2:ha
biom.dat[aoi>800, range(live.MgC.ha.ground, na.rm=T)] ## up to 284 MgC/ha, as we'd expect from what we saw in Raciti
hist(biom.dat[aoi>800,live.MgC.ha.ground]) #v skewed, most below 50 MgC/ha
biom.MgC.ha.ground <- raster(biom)
biom.MgC.ha.ground <- setValues(biom.MgC.ha.ground, biom.dat[,live.MgC.ha.ground])
plot(biom.MgC.ha.ground)

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
summary(biom.dat$npp.ann.ground)
summary(biom.dat$npp.ann.forest)
summary(biom.dat$npp.ann.perv) ## a handful of extra NAs in perv

View(biom.dat[is.finite(npp.ann.ground) & !is.finite(npp.ann.perv),]) ## all partial pix with NA isa, fine
biom.dat[aoi>800 & is.na(npp.ann.ground),] #962 non retreivs, all missing biomass
biom.dat[aoi>800 & is.na(npp.ann.forest),] #962 non retreivs, all missing biomass
biom.dat[aoi>800 & is.na(npp.ann.perv),] #972 non retreivs, all missing biomass
View(biom.dat[aoi>800 & is.na(npp.ann.perv) & !is.na(npp.ann.ground),]) #972 non retreivs, all missing biomass

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
# ### Different FIA factors for different forest types to estimate 30m annual NPP (MgC/yr)
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


