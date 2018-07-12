library(raster)
library(data.table)

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



#######
####
#### FIA V2: empirical NPP~biomass function
### process and assessment of FIA individual dbh records provided by Moreale
### individual FIA tree data from sites near Boston
### V2.2: 1) uses species-specific biometrics; 2) models hardwood growth separately from trees in general; 3) uses nls to avoid dumping negative growth records
### V2.3 MIGHT want to be more careful about excluding sites that are only partially forested
spp.allo <- read.csv("data/FIA/spp_allometrics.csv") ## manually entered selected map from spp to b0+b1
live <- read.csv("data/FIA/MA_Tree_Data_ID_NOMORT.csv")
live <- as.data.table(live)
names(live)[1] <- c("TreeID")
names(live)[2] <- c("PlotID")
spec <- read.csv("data/FIA/REF_SPECIES.csv")
live <- merge(x=live, y=spec[,c("SPCD", "GENUS", "SPECIES")], by.x="SPECIES_CD", by.y="SPCD", all.x=T, all.y=F)
live$GENUS <- as.character(live$GENUS)
live$GENUS <- as.factor(live$GENUS)
live$GENUS.num <- as.numeric(live$GENUS)

### calculate species specific allometrics
live[,spp:=paste(substr(GENUS, 1,1), ".", SPECIES, sep="")]
live <- merge(x=live, y=spp.allo[,c("spp", "b0", "b1")], by="spp", all.x=T)
live[is.na(b0), b0:=(-2.48)]
live[is.na(b1), b1:=2.4835]
biom.pred2 <- function(b0, b1, x){exp(b0+(b1*log(x)))}
live[,biom0.spp:=biom.pred2(b0, b1, DIAM_T0)]
live[,biom1.spp:=biom.pred2(b0, b1, DIAM_T1)]
## class as hard or soft wood
live[,type:="H"]
live[spp%in%c("P.strobus", "P.resinosa", "T.canadensis", "A.balsamea"), type:="S"]
live[,type:=as.factor(type)]

## compare the models of areal biomass growth with raw data and using only hardwoods
live.plot <- live[,.(sum(biom1.spp-biom0.spp, na.rm=T),
                     sum(biom0.spp, na.rm=T),
                     length(DIAM_T0)), by=PlotID] ## we are missing one plot -- all dead?
names(live.plot)[2:4] <- c("biom.growth.spp", "total.biom0.spp.kg", 
                           "num.stems")
### growth rates by hard/soft wood
hwood <- live[type=="H", .(sum(biom1.spp-biom0.spp, na.rm=T), 
                           sum(biom0.spp, na.rm=T)), 
              by=PlotID]
names(hwood)[2:3] <- c("biom.growth.spp.hw", "total.biom0.spp.kg.hw")
swood <- live[type=="S", .(sum(biom1.spp-biom0.spp, na.rm=T), 
                           sum(biom0.spp, na.rm=T)),
              by=PlotID]
names(swood)[2:3] <- c("biom.growth.spp.sw", "total.biom0.spp.kg.sw")

## plot-level growth
### all hard+softwood
live.plot[,growth.ann.rel.spp:=(biom.growth.spp/4.8)/total.biom0.spp.kg]
summary(live.plot$biom.growth.spp) ## a few that decline!
### log-log dummy model (excludes negative growth)
mod.spp <- lm(log(growth.ann.rel.spp)~log(((total.biom0.spp.kg/675)*1E4)), data=live.plot)
summary(mod.spp) ## slightly sloppier, R2 0.19, coefficients about the same
plot(live.plot[,log(((total.biom0.spp.kg/675)*1E4))], live.plot[,log(growth.ann.rel.spp)], cex=0.5, col="forestgreen")
### a marginal shift , doesn't apear to be pivotal
live.plot[growth.ann.rel.spp<0,] ## two plots 255 47 show a decline in biomss! One of them is very low biomass (255)

### Isolate the hardwood growth~total.biom0 if so (few softwoods in urban forest)
live.plot <- merge(live.plot, hwood, by="PlotID")
live.plot <- merge(live.plot, swood, by="PlotID")
live.plot[,growth.ann.rel.hw:=(biom.growth.spp.hw/4.8)/total.biom0.spp.kg.hw]
live.plot[,growth.ann.rel.sw:=(biom.growth.spp.sw/4.8)/total.biom0.spp.kg.sw]
summary(live.plot$growth.ann.rel.hw) # 1.7 to 3% growth for the hardwoods, a slight negative
summary(live.plot$growth.ann.rel.sw) # 1.9 to 4% for the softwoods, some real losses
# ## make sure the distribution of sizes is the same
# par(mfrow=c(2,1))
# hist(live[type=="H", DIAM_T0])
# hist(live[type=="S", DIAM_T0])
# summary(live[type=="H", DIAM_T0]) # 10-84, middle 17-30cm
# summary(live[type=="S", DIAM_T0]) # 9-91, middle 18-37cm -- softwoods are slightly larger if anything

par(mfrow=c(1,1), mar=c(4,3,1,1))
plot(live.plot[,log(((total.biom0.spp.kg/675)*1E4))], live.plot[,log(growth.ann.rel.spp)], cex=0.5, col="gray55", ylim=c(-6, -2))
points(live.plot[,log(((total.biom0.spp.kg/675)*1E4))], live.plot[,log(growth.ann.rel.hw)], cex=0.5, col="red", pch=14)
points(live.plot[,log(((total.biom0.spp.kg/675)*1E4))], live.plot[,log(growth.ann.rel.sw)], cex=0.5, col="blue", pch=16)
## dummy log-log models
mod.hw <- lm(log(growth.ann.rel.hw)~log(((total.biom0.spp.kg/675)*1E4)), data=live.plot[growth.ann.rel.hw>0,])
summary(mod.hw) ## v. miserable, R2 0.06
mod.sw <- lm(log(growth.ann.rel.sw)~log(((total.biom0.spp.kg/675)*1E4)), data=live.plot[growth.ann.rel.sw>0,])
summary(mod.sw) ## v. miserable, R2 0.13
lines(live.plot[growth.ann.rel.hw>0, log(((total.biom0.spp.kg/675)*1E4))], predict(mod.hw), col="red")
lines(live.plot[growth.ann.rel.sw>0, log(((total.biom0.spp.kg/675)*1E4))], predict(mod.sw), col="blue")
## so pines growth faster at low density, but slow down at higher densities harder

## this is the relationship of HARDWOOD relative growth rate to total forest density
## it is significant but does not vary much across the range of densities (37-200 MgC/ha)
## it is also based on a fair amount of data that might not actually be wall-to-wall "forest"
# dim(live.plot[num.stems<20,])[1]/dim(live.plot)[1] ## 35% of our sample might be from areas without full tree coverage
# ##20 stems per 1/6 acre plot is 1 stem per 34 m2, or a minimum of 3000 kg-biomass/ha = 1.5MgC/ha
# hist(live.plot$num.stems)
## so this threshold will have to somehow be found because we intend to treat this as if it indicated a density~growth relationship in 100% canopy

#######
## final exponential model fit, hardwood growth~biomass
plot((live.plot$total.biom0.spp.kg/675)*1E4, live.plot$growth.ann.rel.spp, ylim=c(-0.03, 0.16)) ## bulk growth~biomass up to 400k kg/ha = 200 MgC
points((live.plot$total.biom0.spp.kg/675)*1E4, live.plot$growth.ann.rel.hw, cex=0.4, pch=14, col="red") ## just the hardwoods
summary((live.plot$total.biom0.spp.kg.hw*(1E4/675))/(live.plot$total.biom0.spp.kg*(1E4/675))) ## most plots are more than half hardwood
live.plot[,hw.frac:=total.biom0.spp.kg.hw/total.biom0.spp.kg]
live.plot[hw.frac<0.1,] ## the high ones are all low fraction of HW
live.plot[growth.ann.rel.hw>0.1, hw.frac] ## but there's on that is 40% HW
### fuck it, life is messy
mod.exp.all <- nls(growth.ann.rel.spp ~ exp(a + b * (total.biom0.spp.kg*(1E4/675))),
                   data=live.plot, start=list(a=0, b=0))
mod.exp.hw <- nls(growth.ann.rel.hw ~ exp(a + b * (total.biom0.spp.kg*(1E4/675))),
                  data=live.plot, start=list(a=0, b=0))
v <- summary(mod.exp.all) ## a slightly better model than just hardwoods, everything is sig
w <- summary(mod.exp.hw) ### the b coefficient is barely worth including
x=seq(0,400000)
lines(x, exp(w$coefficients[1]+w$coefficients[2]*x), cex=0.3, col="red")
lines(x, exp(v$coefficients[1]+v$coefficients[2]*x), cex=0.3, col="gray55")
legend(fill=c("red", "gray55"), x = 20000, y = 0.1, legend = c("Hardwoods", "All trees"))
### OK Here's our model to work through the pinche biomass data
### NOTE: this model predicts relative growth (kg/kg) as a function of BIOMASS DENSITY IN KG/HA
## function maxes out at about 3% growth

## reload the biomass data and reprocess
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
biom.dat[, live.kgbiom.ha.ground:=(bos.biom30m/aoi)*(1E4)] ## convert kg-biomass/pixel based on size of pixel (summed by 1m2 in "aoi"), 
biom.dat[,live.kgbiom.ha.forest:=(bos.biom30m/(aoi*can.frac))*(1E4)]
biom.dat[bos.biom30m==0, live.kgbiom.ha.forest:=0] ## have to manually fix this because of 0 canopy pix
biom.dat[can.frac<0.01, live.kgbiom.ha.forest:=0]
biom.dat[,live.kgbiom.ha.perv:=(bos.biom30m/(aoi*(1-isa.frac)))*(1E4)]
biom.dat[bos.biom30m==0, live.kgbiom.ha.perv:=0] ## have to manually fix this because of isa=1 pix
biom.dat[isa.frac>0.99, live.kgbiom.ha.perv:=0]
hist(biom.dat[aoi>800, live.kgbiom.ha.ground]) ## up to 600k kgbiom/ha = 300 MgC/ha
hist(biom.dat[aoi>800, live.kgbiom.ha.forest]) ## same, more medium sized
hist(biom.dat[aoi>800 & isa.frac<0.98, live.kgbiom.ha.perv]) ## up to 800k kgbiom/ha
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
biom.dat[,ground.gfact:=exp(w$coefficients[1]+w$coefficients[2]*live.kgbiom.ha.ground)]
biom.dat[,forest.gfact:=exp(w$coefficients[1]+w$coefficients[2]*live.kgbiom.ha.forest)]
biom.dat[,perv.gfact:=exp(w$coefficients[1]+w$coefficients[2]*live.kgbiom.ha.perv)]
hist(biom.dat$ground.gfact) ## distribution of growth factors with different density calcs
hist(biom.dat$forest.gfact)
hist(biom.dat$perv.gfact)
### the biomass growth rate regression is valid up to ~400k kgbiom/ha

## calculate npp from these growth factors
## regression coeff*biom.density(kg/ha, 3 approaches)-->growth factors (kg/kg) per cell
## growth factors * cell biomass (kg) --> npp (kg biomass per cell) 
biom.dat[,npp.kg.hw.ground:=bos.biom30m*ground.gfact]
biom.dat[,npp.kg.hw.forest:=bos.biom30m*forest.gfact]
biom.dat[,npp.kg.hw.perv:=bos.biom30m*perv.gfact]

### totals for aoi
biom.dat[aoi>800 & is.finite(isa.frac) & is.finite(can.frac), sum(npp.kg.hw.ground, na.rm=T)]/2000 ## 9.1k MgC/yr ground basis
biom.dat[aoi>800 & is.finite(isa.frac) & is.finite(can.frac), sum(npp.kg.hw.forest, na.rm=T)]/2000 ## 8.7k MgC/yr forest basis
biom.dat[aoi>800 & is.finite(isa.frac) & is.finite(can.frac), sum(npp.kg.hw.perv, na.rm=T)]/2000  ## 8.4k MgC/yr perv basis
write.csv(biom.dat, "processed/npp.FIA.empirV22.csv")

# rr <- biom
# rr <- setValues(rr, biom.dat$npp.kg.hw.forest)
# plot(rr)
# plot(biom)
biom.dat[aoi>800, .(median(npp.kg.hw.ground, na.rm=T), ## median is ~50-60 kg-biomass/pix
                    median(npp.kg.hw.forest, na.rm=T),
                    median(npp.kg.hw.perv, na.rm=T))]
