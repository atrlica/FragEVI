library(data.table)
library(reticulate)
library(raster)

hf.dbh <- read.csv("docs/ian/Harvard_Forest_biomass.csv") ## edge plto data for HF
length(unique(hf.dbh$Spp)) ## 18 species ID'd

andy.dbh <- read.csv("docs/ian/Reinmann_Hutyra_2016_DBH.csv") ##edge plot data for Andy's suburban sites
length(unique(andy.dbh$Species_Code)) ## 17 species ID'd

andy.bai <- read.csv("docs/ian/Reinmann_Hutyra_2016_BAI.csv")
dim(andy.bai) ## 195 individuals, dbh by year -- convert to biomass by year? Edge plot data from surburb sites

rope <- read.csv("docs/ian/Rope_Plots_Data_Entry.csv") ## Rope plots, probably not useful, just a description of 

### should we also see if Britt's biomass increment~biomass is worth trying in here?

street <- read.csv("docs/ian/Boston_Street_Trees.csv") ## Street tree resurvey, biomass 2006/2014, species
length(unique(street$Species)) ## 77 species!
street$ann.npp <- (street$biomass.2014-street$biomass.2006)/8
mean.na <- function(x){mean(x, na.rm=T)}
tapply(street$ann.npp, INDEX = street$Species, FUN = mean.na)
tapply(street$ann.npp, INDEX = street$Species, FUN=length)
# genspec <- strsplit(as.character(street$Species), " ")
# gen <- unlist(lapply(genspec, "[[", 1))
# spec <- unlist(lapply(genspec, "[[", 2))
# a <- regexpr("^[[:alpha:]]+", street$Species)
# b <- a+attr(a, "match.length")-1
# unique(gen)
# street$genus <- gen

### clean up taxonomies
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
## dominated by Acer, Gleditsia, Tilia, Zlekova, Ulmus, Fraxinus

plotme <- street[standing.dead.2014==0 & genus %in% focus,]
plotme <- plotme[ann.npp>0,]
plot(street$dbh.2014, street$ann.npp, ylim=c(0, 110))
plot(plotme$dbh.2014, plotme$ann.npp)
plot(plotme$biomass.2006, plotme$ann.npp)
plot(plotme$biomass.2014, plotme$ann.npp) ### this is biomass kg per tree, related to biomass gain per year
biom.top <- plotme[,mean(biomass.2014)+1.96*(sd(biomass.2014))]
biom.bot <- plotme[,mean(biomass.2014)-1.96*(sd(biomass.2014))]
plot(plotme[biomass.2014<biom.top, biomass.2014], plotme[biomass.2014<biom.top, ann.npp])
plot(plotme[biomass.2014<biom.top, dbh.2014], plotme[biomass.2014<biom.top, ann.npp]) ## we don't have dbh

lin.mod <- lm(formula = plotme[biomass.2014<biom.top, ann.npp]~plotme[biomass.2014<biom.top, biomass.2014])
summary(lin.mod)
pol.mod <- lm(formula = plotme[biomass.2014<biom.top, ann.npp]~poly(plotme[biomass.2014<biom.top, biomass.2014], 2))
summary(pol.mod)

### biomass in 2014 doesn't really predict annual growth much in the preceding years, R2 0.43 for poly2

for(g in 1:length(focus)){
  plot(plotme[genus==focus[g], biomass.2014], plotme[genus==focus[g], ann.npp],
       main=focus[g])
  print(summary(lm(plotme[genus==focus[g], ann.npp]~plotme[genus==focus[g], dbh.2014])))
  
}
## everything grows faster at higher dbh, some are sloppy (R2 0.3), some are tighter (0.6-0.8)
## combined top genera have 0.65*dbh2014 for ann.prod, R2 0.46
plot(plotme$biomass.2014, plotme$ann.npp)

## mixed hardwood allometric
test <- 1:100
b0 <- -2.48
b1 <- 2.48
plot(test, exp(b0+(b1*log(test))))

### process biomass to 30m aggregate
### reran modified bos1m_agg to get summed 30m biomass
### need to think about what the raster values mean when summed to 30m grid from 1m

biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m 
biom1m <- raster("data/dataverse_files/bostonbiomass_1m.tif") ## this is in MgC/ha
can1m <- raster("data/dataverse_files/bostoncanopy_1m.tif")
## conversion factors: the summed biomass/30m (900m2) cell 

(51240/1E04) ## a total of about 5 MgC/cell at max (convert MgC/ha to /m2 -- grid cell value is sum of the 1m2 MgC/ha values)
(51240/1E04)*2 ## a total of about 10 Mg.biomass per cell at max
51240*(1E-04)*(1E04/900) ## or about 57 MgC/ha max, in 900m2 cell footprint

## compare to 1m map of biomass (MgC/ha) ## max of 84 MgC/ha at 1m
84*(1E-04)*1000 ### in 1m grid, max of 8.4 kgC/m2
84*(1E-04)*1000*900 ### in 900m2 30m grid cell, 7560 kgC/cell = 7.6 MgC/cell
biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[, npp.ss:=5.104+(bos.biom30m/1E04)*2*1000*0.029] ### based on linear model of street tree growth~
biom.dat[bos.biom30m==0, npp.ss:=0] ## manually correct for places w no trees
npp.ss.r <- raster(biom)
npp.ss.r <- setValues(npp.ss.r, biom.dat$npp.ss) ## ok, just a linear transform, 0-300 kg biomass/cell/yr, equiv to 1.66 MgC/ha/yr

## check: can we get figures same as reported in Raciti et al?
plot(biom1m) ## Raciti lists this as MgC/ha, but it's really Mg-biomass/ha
biom1m.dat <- as.data.table(as.data.frame(biom1m))
names(biom1m.dat) <- c("biom.raw")
biom1m.dat[, MgC.ha.corr:=2*biom.raw]
biom1m.dat[, Mgbiom:=biom.raw/(10^4)]
biom1m.dat[, mean(Mgbiom, na.rm=T)]
mean.biom.raw <- biom1m.dat[, mean(biom.raw, na.rm=T)] ## 5.75 Mgbiom/ha
sum.biom.raw <- biom1m.dat[, sum(biom.raw, na.rm=T)]
area <- biom1m.dat[!is.na(biom.raw), length(biom.raw)] ## 124 km2
(340000-315000)*(4695000-4680000)/(10^6) ## size of raw raster approx: 3.75E8 = 375km2
mean.biom.raw*(area/10^4) ## 71412 Mgbiom/ha
sum.biom.raw/(10^4) ## 71.4k Mgbiom
(sum.biom.raw/(10^4))/(area/(10^4)) ## 5.75 Mg biom/ha

## biomass formula, basically just lidar canopy height linear transform h(m)-->bm(kg)
biomass=(2.1015*h)+0.8455
biomass=82 ## max biomass in this 1m raster
h=27 ### actual canopy height in the highest biomass pixels (brief visual search of a single area) is about 27m
(2.1015*h)+0.8455 ### a single pixel of 27m height gives biomass of 57kg (this is consistent with the biom1m raster value)
biomass=57
## at max biomass/pixel
biomass=82
(biomass-0.8455)/2.1015 ### this would imply a 39m canopy height

### test, single open-grown tree
## a tree of 14m diameter, max height ~19m, area = 153m, ~30 kg biomass per pixel average
153*30 ## about 4.5 Mg biomass
### actual extract
library(rgeos)
library(rgdal)
tree.foot <- readOGR(dsn="processed/boston/tree.footprint.shp", layer="tree.footprint")
crs(tree.foot)
crs(biom1m)
tree.foot <- spTransform(tree.foot, crs(biom1m))
### this is a open-grown tree, 14m diameter canopy roughly circular, max height is 19m
a <- extract(biom1m, tree.foot)
a <- unlist(a)
sum(a) ### 3941 kg biomass in this tree
## general equation for biomass (Ian street tree)
biomass=exp(-2.48+(2.4835(log(dbh))))
biomass=sum(a)
exp((log(biomass)+2.48)/2.4835) ## this would imply a 76cm dbh tree with 19m height
## at 19m high, this tree, according to the transform in Raciti
h=19
(2.1015*h)+0.8455 ## about 40 kg/pixel
40*153 ## roughly correct, somewhat too high (max pixel height)

### CONCLUSION: RACITI'S 1m BIOMASS DATA is correct transform of lidar height, yielding kg-biomass/1m.cell
sum.biom.raw <- biom1m.dat[, sum(biom.raw, na.rm=T)]
area <- biom1m.dat[!is.na(biom.raw), length(biom.raw)] ## 124 km2
(sum.biom.raw/10^3)/(area/(10^4)) ### average of 57.5 Mgbiomass/ha
(sum.biom.raw/10^3)/(area/(10^4))/2 ### average of 28.8 MgC/ha --> same as raciti
(sum.biom.raw/10^3)/2/(10^3) ### 357 GgC total storage --> same as raciti
### the scale on the corrected figures from Raciti is what you get when you take
### kgbiomass/m2 * 1E4m2/ha * 1kgC/2kgbiomass * 1MgC/1E3kgC (e.g. 82 kgbiomass/m2 --> 410 MgC/ha)

plot(can1m)
can1m.dat <- as.data.table(as.data.frame(can1m))
names(can1m.dat) <- c("can")
can.totcan <- can1m.dat[!is.na(can), sum(can)]
can.totar <- can1m.dat[!is.na(can), length(can)]
can.totcan/can.totar ## avg cover 25.8%, contrast 25.5% reported 
### summary: canopy data looks about like what was reported in Raciti



#### OK with that out of the way: Let's try to figure out the tree distribution per 30m cell and uptake
biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m 
plot(biom)
street <- read.csv("docs/ian/Boston_Street_Trees.csv") ## Street tree resurvey, biomass 2006/2014, species
length(unique(street$Species)) ## 77 species!
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
focus

### approach 0: get npp biomass increment as function of median tree biomass, FAI growth rates
## take median value of tree dbh, ignore genus (using inline biomass ~ common hardwood allometric)
dbh <- street[standing.dead.2014==0 & genus %in% focus, median(dbh.2014, na.rm=T)]
biomass <- street[standing.dead.2014==0 & genus %in% focus, median(biomass.2014, na.rm=T)]
hist(street[standing.dead.2014==0 & genus %in% focus, biomass.2014]) ## highly skewed towards <1000kg
biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[,num.trees.med:=bos.biom30m/biomass]
biom.dat[,range(num.trees.med, na.rm=T)] ## max 280 trees per cell i.e. a tree per every 3m2

## the FIA use this equation for volume of wood in a stand by stand age
a=123.67
b=0.04
AGE=40
y = a*(1-exp(-b*AGE))^3
y
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

biom.dat
biom.dat[,Mgbiom.ha:=((bos.biom30m/900)*1E4)/1E3]
hist(biom.dat[,Mgbiom.ha]) ## skewed, most below 100 Mgbiom/ha (~50 MgC/ha), max 569 = 284 MgC/ha, contrast: Max 1m = 469 MgC/ha
biom.dat[,MgC.ha:=Mgbiom.ha/2]

## what is relationship between standing live biomass-C and C increment
plot(ne.summ$mean.vol.m3, ne.summ$live.tree.c.inc)
plot(ne.summ$live.tree.tCha, ne.summ$live.tree.c.inc)
plot(ne.summ$Age.yrs, ne.summ$live.tree.tCha)
plot(ne.summ$Age.yrs, ne.summ$mean.vol.m3) ## basic sigmoid 0 to max at 100

### this is the controlling volume/live biomass C equation for all the figures produced in the FIA report
## coefficients below apply to NE total summary
## live C
a=123.67
b=0.04
x=seq(0,200) ## age of stand
y=a*(1-(exp(-b*x)))^3 
plot(x,y)



##########
### prototype process for using FIA aggregate data to figure npp from Raciti biomass
# 1) determine MgC/ha in 30m plot, normalize to %canopy (i.e. there is X many ha of forest there with Y much living carbon in place)
# ie. as far as FIA is concerned, that is a forest, in the area it occupies, that has Y much AGB
# 2) figure out how old FIA would say a plot like that is
# 3) Figure out the next-year tC/ha in the plot
# 4) apply this incrememnt to the area of "forest" (canopy) in the pixel

## example: X=40 yrs,
x=40
a=123.67
b=0.04
y=a*(1-(exp(-b*x)))^3 
y ## 62.9 tC/ha
### what is equivalent age of stand if we know live biomass?
log(1-(y/a)^(1/3))/(-b) ## predicts 40 years

## what is age of a forest with 100% canopy and 30000 kg/pixel biomass
y=30000/1000/2*(10^4)/900 ### 166 MgC/ha --> out of bounds)
y=20000/1000/2*(10^4)/900  ### at 30m pixel values over about 20000 kg biomass, you're out of bounds for stand age -- too thick with biomass!

biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m kg-biomass to 30m pixel
aoi <- raster("processed/boston/bos.aoi30m.tif")
plot(aoi)
plot(biom)
biom <- crop(biom, aoi) ## biomass was slightly buffered, need to clip to match canopy fraction raster
biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[,aoi:=as.vector(getValues(aoi))]
### live MgC by area in each pixel
biom.dat[, live.MgC.ha:=(bos.biom30m/aoi)*(1/2)*(1/1000)*(10^4)] ## convert kg-biomass/pixel based on size of pixel (summed by 1m2 in "aoi"), kgC:kgbiomass, Mg:kg, m2:ha
biom.dat[, range(live.MgC.ha, na.rm=T)] ## up to 284 MgC/ha, as we'd expect from what we saw in Raciti
biom.MgC.ha <- raster(biom)
biom.MgC.ha <- setValues(biom.MgC.ha, biom.dat[, live.MgC.ha])
plot(biom.MgC.ha)
hist(biom.dat[,live.MgC.ha]) #v skewed, most below 50 MgC/ha

## correct biomass figures for the amount of canopy present per cell
### i.e. we assume FIA is measuring trees in "forest" with essentially continuous canopy, so that differences in tC/ha as a function of age are purely a function of tree growth and not differences in tree %coverage
can <- raster("processed/boston/bos.can30m.tif")
can.dat <- as.data.table(as.data.frame(can))
biom.dat[, can.frac:=can.dat$bos.can30m]
### live MgC/ha for "forest" fraction in each pixel
biom.dat[,live.MgC.ha.forest:=live.MgC.ha/can.frac]
biom.dat[can.frac<=0.005, live.MgC.ha.forest:=0]
range(biom.dat[,live.MgC.ha.forest],  na.rm=T)
hist(biom.dat[,live.MgC.ha.forest]) ## correcting for canopy cover, most pixels are <100 MgC/ha (good -- the FIA only really handles stuff up to ~120 MgC/ha, about 100 years old)
biom.dat[live.MgC.ha.forest<100, length(live.MgC.ha.forest)]/dim(biom.dat[!is.na(can.frac),])[1] ## 84% of pixels are below 100 MgC/ha
biom.MgC.ha.forest <- raster(biom)
biom.MgC.ha.forest <- setValues(biom.MgC.ha.forest, biom.dat[,live.MgC.ha.forest])
plot(biom.MgC.ha.forest)

### figure out forest "age" for the cells (using coefficients for NE total)
a=123.67
b=0.04
biom.dat[,age:=log(1-(live.MgC.ha.forest/a)^(1/3))/(-b)] ## some of these will produce infinities with too dense biomass
biom.dat[age>100, age:=100] ## fix the divergent ones to just "old, not growing"
biom.dat[!is.finite(age), age:=100] ## again, fix the ones that got fucked to "old, not growing"
biom.dat[is.na(aoi), age:=NA] # cancel places out of bounds
biom.dat[bos.biom30m<=10, age:=NA] ## cancel places with no biomass
biom.dat[is.na(bos.biom30m), age:=NA]

forest.age <- raster(biom)
forest.age <- setValues(forest.age, biom.dat[,age])
plot(forest.age)
hist(biom.dat[,age]) ## got a lot of "old" ones, indicating high density of biomass

# ### test to see the range of this shit
# x=seq(1:500)
# plot(x, log(1-(x/a)^(1/3))/(-b)) ## dies after about 120

## figure out next annual increment possible for the "forest" average present in each cell, based on projected "age"
biom.dat[,npp.ann:=(a*(1-(exp(-b*(age+1))))^3)-(a*(1-(exp(-b*(age))))^3)]
## note these units are in MgC/ha/yr
npp.ann <- raster(biom)
npp.ann <- setValues(npp.ann, biom.dat[,npp.ann])
plot(npp.ann)

## now correct for the size of the pixel
biom.dat[, npp.ann.actual:=(npp.ann/10^4)*900*can.frac]
hist(biom.dat[,npp.ann.actual])
npp.actual <- raster(biom)
npp.actual <- setValues(npp.actual, biom.dat[,npp.ann.actual]) ## units are in MgC/yr for each pixel
plot(npp.actual)
biom.dat[,sum(npp.ann.actual, na.rm=T)] ## 4887 MgC/yr citywide
biom.dat[,sum(aoi, na.rm=T)/(10^4)] ## city area in ha
biom.dat[,sum(npp.ann.actual, na.rm=T)]/biom.dat[,sum(aoi, na.rm=T)/(10^4)]
## predit about 0.39 MgC/ha/yr NPP across the city: contrast 10.3-8.9 = 1.4 NEP for boston region in Hardiman
writeRaster(npp.actual, filename="processed/boston/bos.npp.actual.fia.NEdefault.tif", format="GTiff", overwrite=T)


######################
### FIA factors to estimate 30m annual NPP (MgC/yr)

library(data.table)
library(raster)

## cleaned up code, handles different growth parameters for different forests (note "Afforestation" and "Reforestation" values are same viz. live biomass growth rates)
## initialize parameters for different forest types that ?? resemble tree species distributions in Boston
for.type <- c("NEdefault","Mixoak", "Redoak")
a <- c(123.67, 130.81, 123.09)
b <- c(0.04, 0.03, 0.04)

# ## test limits of the live biomass~stand age function
# for(f in 1:length(for.type)){
#   x=seq(0,150)
#   y=a[f]*(1-(exp(-b[f]*x)))^3
#   plot(x,y, main=for.type[f], ylab="live tC/ha", xlab="stand age")
#   x <- seq(0, 140) ## inverse: model age vs. live biomass
#   y=log(1-(x/a[f])^(1/3))/(-b[f]) ## 
#   plot(x, y, main=for.type[f], ylab="live tC/ha", xlab="stand age")
# }
# ## conservatively, none of the models is particularly stable over 100 yrs stand age

## as before, procedure is:
# 1) scale up pixel biomass to MgC/ha 
# 2) adjust live C Mg/ha by %canopy of pixel (i.e. MgC ha *in forested area of pixel*)
# 3) figure out how old FIA would say a forest with that much C/ha is; set anything with age>100 to 100 (old, not growing)
# 4) Figure out the next-year tC/ha in the plot
# 5) apply this incrememnt to the area of "forest" (canopy) in the pixel, adjust by pixel size

biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m kg-biomass to 30m pixel
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- crop(biom, aoi) ## biomass was slightly buffered, need to clip to match canopy fraction raster
biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[,aoi:=as.vector(getValues(aoi))]
biom.dat[aoi<800, bos.biom30m:=NA] ### cancel values in cells that are not at least ~90% complete coverage in AOI

### live MgC by area in each pixel
biom.dat[, live.MgC.ha:=(bos.biom30m/aoi)*(1/2)*(1/1000)*(10^4)] ## convert kg-biomass/pixel based on size of pixel (summed by 1m2 in "aoi"), kgC:kgbiomass, Mg:kg, m2:ha
# biom.dat[, range(live.MgC.ha, na.rm=T)] ## up to 284 MgC/ha, as we'd expect from what we saw in Raciti
# biom.MgC.ha <- raster(biom)
# biom.MgC.ha <- setValues(biom.MgC.ha, biom.dat[, live.MgC.ha]) 
# plot(biom.MgC.ha, main="biomass Mg/ha, whole pixels") ## compare to Raciti, looks about right
# hist(biom.dat[,live.MgC.ha]) #v skewed, most below 50 MgC/ha -- so resembles either very young forest (mostly below 30 yrs) or fragmented forest with incomplete coverage

## correct biomass figures for the amount of canopy present per pixel
can <- raster("processed/boston/bos.can30m.tif")
biom.dat[, can.frac:=getValues(can)]
biom.dat[,live.MgC.ha.forest:=live.MgC.ha/can.frac] ### live MgC/ha for "forest" fraction in each pixel
biom.dat[can.frac<=0.01, live.MgC.ha.forest:=0] ## for pixels with very little canopy coverage, cancel biomass effect altogether (otherwise tiny biomass+tiny canopy could imply super old tiny forest patch)
# range(biom.dat[,live.MgC.ha.forest],  na.rm=T)
# hist(biom.dat[,live.MgC.ha.forest]) ## correcting for canopy cover, most pixels are <100 MgC/ha (good -- the FIA only really handles stuff up to ~120 MgC/ha, about 100 years old)
# biom.dat[live.MgC.ha.forest<120, length(live.MgC.ha.forest)]/dim(biom.dat[!is.na(can.frac),])[1] ## 90% of pixels are below 120 MgC/ha, where the thing starts to fall apart
# biom.MgC.ha.forest <- raster(biom)
# biom.MgC.ha.forest <- setValues(biom.MgC.ha.forest, biom.dat[,live.MgC.ha.forest])
# plot(biom.MgC.ha.forest, main="biomass MgC/ha, forest fraction of pixel")

### back-project forest "age" for the cells, calculate annual increment based on age, readjust to pixel dimensions and %canopy
for(f in 1:length(for.type)){
  biom.dat[,age:=log(1-(live.MgC.ha.forest/a[f])^(1/3))/(-b[f])] ## some of these will produce infinities with too dense biomass
  biom.dat[age>100, age:=100] ## fix the divergent ones to just "old, not growing"
  biom.dat[!is.finite(age), age:=100] ## again, fix the ones that got fucked to "old, not growing"
  biom.dat[is.na(aoi), age:=NA] # cancel places out of bounds
  biom.dat[bos.biom30m<=10, age:=NA] ## cancel places with very little biomass
  biom.dat[is.na(bos.biom30m), age:=NA]
  # forest.age <- raster(biom)
  # forest.age <- setValues(forest.age, biom.dat[,age])
  # plot(forest.age) ## places with a fair amount of dense biomass cover look like stable non-growing forests
  # hist(biom.dat[,age]) ## peak in the 20-40 range, but a lot of "old" forest pixels, indicating high density of biomass
  ## figure out next annual increment possible for the "forest" average present in each cell, based on projected "age"
  biom.dat[,npp.ann:=(a[f]*(1-(exp(-b[f]*(age+1))))^3)-(a[f]*(1-(exp(-b[f]*(age))))^3)] ## note these units are in MgC/ha/yr, treating the whole pixel as a forest
  # npp.ann <- raster(biom)
  # npp.ann <- setValues(npp.ann, biom.dat[,npp.ann]) ## reduced in the densely forested parts compared to the mod.developed; nothing in the very developed)
  # plot(npp.ann, main="Annual NPP tC/ha/yr, whole pixel basis")
  
  ## correct for the size and forested area of the pixel
  biom.dat[, npp.ann.actual:=(npp.ann/10^4)*aoi*can.frac] ## this is now MgC/yr for each pixel
  # hist(biom.dat[,npp.ann.actual]) ## mostly below 0.1 MgC/yr (i.e. 100 kgC per pixel, compare to >30000 kg biomass in some pixels)
  npp.actual <- raster(biom)
  npp.actual <- setValues(npp.actual, biom.dat[,npp.ann.actual]) ## units are in MgC/yr for each pixel
  plot(npp.actual, main=paste("Annual NPP MgC/yr,", for.type[f], ", forested area of each pixel"))
  # biom.dat[,sum(npp.ann.actual, na.rm=T)] ## total MgC/yr uptake in NPP, citywide
  # biom.dat[,sum(aoi, na.rm=T)/(10^4)] ## city area in ha
  # biom.dat[,sum(npp.ann.actual, na.rm=T)]/biom.dat[,sum(aoi, na.rm=T)/(10^4)] ## avg citywide NPP
  writeRaster(npp.actual, format="GTiff", overwrite=T,
              filename=paste("processed/boston/bos.npp.actual.fia", for.type[f], "tif", sep="."))
}
### immediate post-hoc notes: The dense forest parts can be
## 1) fairly low in NPP -- they exceed the threshold of biomass density, get classed as "very old", and have low annual increment even though canopy cover is high
## 2) high in NPP -- they have high canopy cover and moderate biomass so they have the most capacity for uptake ove the greatest area



##########
##########
### Approach 2: Per-cell growth rate derived from local street tree increment
### 2a: Assume distribution of trees that are all median dbh/biomass
library(data.table)
library(raster)

# figure out median dbh/biomass of a street tree
street <- read.csv("docs/ian/Boston_Street_Trees.csv") ## Street tree resurvey, biomass 2006/2014, species
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
# sum(a$V1[1:12]) ### top 12 is 94% of trees present
focus <- a$genus[1:12]

plot(street$dbh.2006, street$biomass.2006)
plot(street$dbh.2014, street$biomass.2014)
street[dbh.2014<=45, length(dbh.2014)]/street[!is.na(dbh.2014), length(dbh.2014)] #80% of trees below 45cm dbh
hist(street$dbh.2014) # steep drop north of 40cm
hist(street$dbh.2006) ## looks somehow higher than 2014
post.dbh <- street[dead.by.2014==0, median(dbh.2014, na.rm=T)] #28.4cm
street[dead.by.2014==0, median(biomass.2014, na.rm=T)] #172 kg/tree
### working biomass eq (generalized Eastern hardwood, Smith et al. 2018)
st.biom.med <- exp(-2.48+2.4835*log(post.dbh)) # the biomass of a median tree is 341kg

## determine the number of standard street trees present in each cell
biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m kg-biomass to 30m pixel
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- crop(biom, aoi) ## biomass was slightly buffered, need to clip to match canopy fraction raster
biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[,aoi:=as.vector(getValues(aoi))]
biom.dat[aoi<800, bos.biom30m:=NA] ### cancel values in cells that are not at least ~90% complete coverage in AOI

### live MgC by area in each pixel
biom.dat[, live.MgC.ha:=(bos.biom30m/aoi)*(1/2)*(1/1000)*(10^4)] ## convert kg-biomass/pixel based on size of pixel (summed by 1m2 in "aoi"), kgC:kgbiomass, Mg:kg, m2:ha
biom.dat[,num.streettree.med:=bos.biom30m/st.biom.med]
hist(biom.dat[,num.streettree.med]) ## right skewed, median 5.8 trees/pixel

######
### check calculations and data quality for street tree data
dim(street) # 4259 records
dim(street[is.finite(dbh.2014) & is.finite(dbh.2006),]) # 2681 records with both measurements
dim(street[is.finite(dbh.2014) & is.finite(dbh.2006) & dead.by.2014==0,]) # 2676 records with both live measurements
street[, delta.dbh:=dbh.2014-dbh.2006]
range(street$delta.dbh, na.rm=T) ## up to 23cm growth, some negative
dim(street[is.finite(dbh.2014) & is.finite(dbh.2006) & dead.by.2014==0 & delta.dbh>=0,]) # 2467 records with seemingly valid growth recorded
street[, record.good:=0]
street[is.finite(dbh.2014) & is.finite(dbh.2006) & dead.by.2014==0 & delta.dbh>=0, record.good:=1] # 2467 good growth records marked
street[Pruning==2, record.good:=0] ## down to 2161 trees that haven't been hacked to pieces
plot(street[record.good==1, dbh.2006], street[record.good==1, dbh.2014])
abline(a = 0, b=1, col="red") ## most obvious growth is taking place in the 0-80cm range
hist(street[record.good==1, dbh.2014]) #reasonably broad peak 20-40cm, skews to 120
dim(street[record.good==1 & dbh.2014<=45,])/2161 ## 77% of records below 45cm

### check the dbh calculations against published equation in Smith
b0 <- -2.48
b1 <- 2.4835 ## these are eastern hardwood defaults
biom.pred <- function(x){exp(b0+(b1*log(x)))}
biom.pred(20)
street[record.good==1, at.biom.2014:=biom.pred(dbh.2014)]
street[record.good==1, at.biom.2006:=biom.pred(dbh.2006)]
plot(street$dbh.2006, street$at.biom.2006) ## seems legit
plot(street$dbh.2006, street$biomass.2006) ## legit, but slower growth
plot(street$biomass.2006, street$at.biom.2006) ## a different equation was used for these (at = higher slope)
abline(a=0, b=1, col="red")
plot(street$biomass.2006*2, street$at.biom.2006) ## a different equation was used for these (at = higher slope)
abline(a=0, b=1, col="red") ### nm, Ian's figure are just kgC instead of kgbiomass

plot(street$dbh.2014, street$at.biom.2014) ## seems legit
plot(street$biomass.2014, street$at.biom.2014) ## a different equation was used for these (at = higher slope)
abline(a=0, b=1, col="red")
## check do you get the same median dbh and matching biomass.2014?
street[record.good==1, median(dbh.2014, na.rm=T)] ## 31.7cm
hist(street[record.good==1, dbh.2014])
hist(street[record.good==1, dbh.2006]) ## crept up the distribution a bit in 8 years
street[record.good==1, median(biomass.2014, na.rm=T)] ## median tree is 225.5 kg in Ian's data
street[record.good==1, median(at.biom.2014, na.rm=T)] ## median tree is 447.5 kg in my reckoning
biom.pred(31.7) ## 447.5kg is what the E. hardwoods equation would predict at median dbh
## what is that growth looking like
hist(street[record.good==1, delta.dbh]) ## like ~5-10 cm increase for most
hist(street[record.good==1, dbh.2006])
hist(street[record.good==1, at.biom.2006]) ## way more skewed than the dbh
biom.pred(60)
hist(street[record.good==1 & dbh.2006<=60, at.biom.2006]) ## still super skewed, most trees < 500kg
street[record.good==1 & dbh.2006<=60, length(at.biom.2006)]/street[record.good==1, length(at.biom.2006)] ### 95% of tree individuals are <60cm
street[record.good==1 & dbh.2006<=60, sum(at.biom.2006)]/street[record.good==1, sum(at.biom.2006)] ## but only 68% of total biomass is <60cm (i.e. largest 5% of trees is 32% of biomass)


#######
### OK, working from the equation Ian used in manuscript, determine relationship between biomass and npp
street[record.good==1, at.delta.biom:=at.biom.2014-at.biom.2006]
street[record.good==1, at.npp.ann:=at.delta.biom/8]

### looking at 2014 recollect data
street[record.good==1, median(at.biom.2014)] #446 kg
street[record.good==1, median(dbh.2014)] #31.7cm
street[record.good==1 & dbh.2014<=60, length(at.biom.2006)]/street[record.good==1, length(at.biom.2014)] ### 95% of tree individuals are <60cm
street[record.good==1 & dbh.2014<=60, sum(at.biom.2006)]/street[record.good==1, sum(at.biom.2014)] ## but only 68% of total biomass is <60cm (i.e. largest 5% of trees is 32% of biomass)

## ann.npp~dbh.2014
plot(street[record.good==1, dbh.2014], street[record.good==1, at.delta.biom]) ## whacktastic
plot(street[record.good==1, dbh.2014], street[record.good==1, at.npp.ann]) ## looks a little non-linear
summary(lm(at.npp.ann~dbh.2014, data=street[record.good==1,])) ## R2 = 0.47, might improve with an exponential form 
## filter for huge trees
summary(lm(at.npp.ann~dbh.2014, data=street[record.good==1 & at.biom.2014<=2700,])) ## R2=0.44, slope and intercept similar to full-data case
## can try a 2nd order poly here if curious

## ann.npp~biom.2014
plot(street[record.good==1, at.biom.2014], street[record.good==1, at.npp.ann]) ## linear, serious outliers/leverage
summary(lm(at.npp.ann~at.biom.2014, data=street[record.good==1,])) ## R2=0.40
street[record.good==1, mean(at.biom.2014)+1.96*sd(at.biom.2014)] ## upper 5% of biomass is ca. >2700
# filter for huge trees
summary(lm(at.npp.ann~at.biom.2014, data=street[record.good==1 & at.biom.2014<=2700,])) ## R2=0.43 (no big improvement), slope is x2 without including giants, but covers all trees <65cm (above this do we even trust the biomass equations?)

## same for 2006 measures
street[record.good==1, median(at.biom.2006)] #198.7kg
street[record.good==1, median(dbh.2006)] #22.9cm

# at.npp.ann~dbh.2006
plot(street[record.good==1, dbh.2006], street[record.good==1, at.delta.biom]) ## pretty whacktastic
summary(lm(at.npp.ann~dbh.2006, data=street[record.good==1,])) ## R2=0.25, at least 5kg/yr, b1=0.93 kg/cm but low predictive power, R2=0.25
## try with smaller trees
summary(lm(at.npp.ann~dbh.2006, data=street[record.good==1 & at.biom.2014<=2700,])) ## even worse, R2=0.18

# at.npp.ann~biom.2014
plot(street[record.good==1, at.biom.2006], street[record.good==1, at.npp.ann]) ## lots of huge outliers
summary(lm(at.npp.ann~at.biom.2006, data=street[record.good==1,])) ## R2=0.20
summary(lm(at.npp.ann~at.biom.2006, data=street[record.good==1 & at.biom.2006<2700,])) ## R2=0.15

## so.... 2006 dbh is an apparently bad predictor of future growth rate
## 2014 dbh is a better predictor of past growth rate
## is this because either the data from 2006 is sketchy or because growth is more variable than we'd like
## i.e. it's easy to tell that the big trees grew faster by the time you see them big, but you couldn't tell how fast they would grow earlier
## so which relationship(s) to use: 2014 back projection or 2006 forward (crappy)
## and does Raciti's tree map more resemble the trees measured in 2005/6 or 2014?
## his quickbird data was 2006/7, paired with LiDAR 2005
## fuck 2006 is probably our guy
x=0:100
y=5/x
plot(x, y)
street[record.good==1, dbh.inc:=delta.dbh/8]
street[record.good==1, dbh.inc.rel:=dbh.inc/dbh.2006]
plot(street[record.good==1, dbh.2006], street[record.good==1, dbh.inc]) # possibly hyperbolic?
plot(street[record.good==1, dbh.2006], street[record.good==1, dbh.inc.rel]) # def looks hyperbolic-ish

## 1/log of both looks good
plot(street[record.good==1 & dbh.2006<100, log(1/dbh.2006)], street[record.good==1 & dbh.2006<100, log(1/dbh.inc.rel)])
## hyperbolic, y=x/(a+Bx), not as good
plot(street[record.good==1 & dbh.2006<100, (1/dbh.2006)], street[record.good==1 & dbh.2006<100, (1/dbh.inc.rel)], ylim=c(0,500)) ##hyperbolic transform but looks fucked

## log-log looks ok (exponential function/transform, y=a*eBx, B<0)
plot(street[record.good==1 & dbh.2006<65, log(dbh.2006)], street[record.good==1 & dbh.2006<65, log(dbh.inc.rel)], main="Linear %change of dbh.2006")
mod.rel <- summary(lm(log(dbh.inc.rel)~log(dbh.2006), data=street[record.good==1 & dbh.inc.rel!=0,])) #r2=0.50!
plot(lm(log(dbh.inc.rel)~log(dbh.2006), data=street[record.good==1 & dbh.inc.rel!=0,]))

## ie. log(dbh relative annual increment) = a+b*log(dbh.2006)
a=mod.rel$coefficients[1]
b=mod.rel$coefficients[2]
## test: if you plug in a dbh.2006 do you get a sensible value out for dbh increment
exp(a+(b*log(30))) ## 30cm dbh in 2006 predicts a 2.1% dbh relative increment in 30cm trees
street[record.good==1 & dbh.2006>28 & dbh.2006<32, mean(dbh.inc.rel, na.rm=T)] ## mean is 2.6%
exp(a+b*log(10)) ## 8.7% dbh increment for 10cm tree in 2006
street[record.good==1 & dbh.2006>8 & dbh.2006<12, mean(dbh.inc.rel, na.rm=T)] ## mean is 11.4%
exp(a+b*(log(60))) ## 0.85 % dbh inrement for 60cm tree in 2006
street[record.good==1 & dbh.2006>58 & dbh.2006<62, mean(dbh.inc.rel, na.rm=T)] ## mean is 0.93%

## second test: will an annual dbh relative increment get you from dbh.2006 to dbh.2014
finish=street[record.good==1,dbh.2014][1] # 32.3
start=street[record.good==1,dbh.2006][1] # 15.24
grow=street[record.good==1,dbh.inc.rel][1] ## 0.1399
(grow*8*start)+start ## ok this predicts fine

## what if there's a compounding factor for growth in form T=T0*(1+grow)^t, where t=years between interval
t=8
finish=street[record.good==1,dbh.2014][5] # 13.9
start=street[record.good==1,dbh.2006][5] # 5.08
grow=((finish/start)^(1/t))-1 ## 0.134 compounding factor
start*((1+grow)^8) ## nails the finish perfectly

## now, does compounding factor (intrinsic rate of growth) track to dbh.2006 any better?
street[record.good==1, dbh.inc.comp:=((dbh.2014/dbh.2006)^(1/8))-1]
plot(street[record.good==1, dbh.2006], street[record.good==1, dbh.inc.comp]) ## this seems a cleaner exponential function
plot(street[record.good==1, log(dbh.2006)], street[record.good==1, log(dbh.inc.comp)], main="Compound of dbh.2006") ## this seems a cleaner exponential function
mod.comp <- summary(lm(log(dbh.inc.comp)~log(dbh.2006), data=street[record.good==1 & dbh.inc.comp!=0,])) ## r2=0.48
plot(lm(log(dbh.inc.comp)~log(dbh.2006), data=street[record.good==1 & dbh.inc.comp!=0,]))
## it's quite similar to the relative %change over time. %change works slightly better to predict future dbh.2014

## compare models of dbh growth~dbh.2006 to transformed biomass increment ### Doesn't perform as well as dbh relative increment~dbh.2006
plot(street[record.good==1 & dbh.2006<65, log(dbh.2006)], street[record.good==1 & dbh.2006<65, log(at.npp.ann)])
summary(lm(log(at.npp.ann)~log(dbh.2006), data=street[record.good==1 & at.npp.ann!=0])) #r2=0.27

### compare to models based on annualized biomass increment %change (i.e. a linear "short term" increase/yr as a fraction of starting biomass)
street[record.good==1, at.npp.ann.rel:=at.npp.ann/at.biom.2006]
hist(street$at.npp.ann.rel) ## most are below 1 factor
plot(street[record.good==1, at.biom.2006], street[record.good==1, at.npp.ann.rel]) ## looks like a pretty serious exponential function
## the tiny trees grew fastest as a fraction of biomass
hist(street[record.good==1, at.npp.ann]) ## most trees <50kg/yr
hist(street[record.good==1, at.npp.ann.rel]) ## most trees <100% biomass gain/yr
### log-log transform
plot(street[record.good==1, log(at.biom.2006)], street[record.good==1, log(at.npp.ann.rel)]) ## looks a lot like the dbh.increment~dbh.2006
mod.biom.rel <- summary(lm(log(at.npp.ann.rel)~log(dbh.2006), data=street[record.good==1 & ann.npp!=0])) #r2=0.54! so relative %change biomass per year is the most predictive so far

#test: does this predict future biomass with any sensitivity
## ie. log(dbh relative annual increment) = a+b*log(dbh.2006)
a=mod.biom.rel$coefficients[1]
b=mod.biom.rel$coefficients[2]
## test: if you plug in a dbh.2006 do you get a sensible value out for ann.npp %change per year
exp(a+(b*log(30))) ## 30cm dbh in 2006 predicts a 6.3% annualized increase in biomass
street[record.good==1 & dbh.2006>28 & dbh.2006<32, mean(at.biom.2006, na.rm=T)] ## mean biomass is 404kg
street[record.good==1 & dbh.2006>28 & dbh.2006<32, mean(at.npp.ann.rel, na.rm=T)] ## mean annualized %change is 8.1%
exp(a+(b*log(10))) ## 10cm dbh in 2006 predicts a 36.2% annualized increase in biomass
street[record.good==1 & dbh.2006>8 & dbh.2006<12, mean(at.biom.2006, na.rm=T)] ## mean biomass is 26.8kg
street[record.good==1 & dbh.2006>8 & dbh.2006<12, mean(at.npp.ann.rel, na.rm=T)] ## mean annualized %change is 56.7%
exp(a+(b*log(60))) ## 10cm dbh in 2006 predicts a 2.1% annualized increase in biomass
street[record.good==1 & dbh.2006>58 & dbh.2006<62, mean(at.biom.2006, na.rm=T)] ## mean biomass is 2190kg
street[record.good==1 & dbh.2006>58 & dbh.2006<62, mean(at.npp.ann.rel, na.rm=T)] ## mean annualized %change is 2.5%
### note there is some real tiny shit still in here -- like ~70 records below 5cm dbh

## second test: can you use ann.npp.rel to get back to biomass.2014
finish=street[record.good==1,at.biom.2014][1] # 468.8kg
start=street[record.good==1,at.biom.2006][1] # 72.6kg = 15.24cm, 392.3 kg gain
grow=street[record.good==1,at.npp.ann.rel][1] ## 68% biomass annualized growth
(grow*8*start)+start ## ok this delivers you back to the 2014 biomass
## best model so far: log-log at.npp.ann.rel~dbh.2006 (annualized %biomass change vs 2006 biomass)

## ok, so approach 2a: 
# 1) assume median tree dbh 
# 2) determine number of trees per cell to give biomass total, E hardwood default
# 3) figure annualized median % biomass increment based on biomass increment in median tree (or random/empirical distribution around median tree)
street[record.good==1 & dbh.2006>10 & dbh.2006<60, median(at.npp.ann.rel, na.rm=T)] ## median 6.4% annualized
hist(street[record.good==1 & dbh.2006>10 & dbh.2006<60, at.npp.ann.rel]) ## skewed, most below 20%, caps at 50
street[record.good==1, median(at.npp.ann.rel, na.rm=T)] ## among all trees, median is 10.2% annualized
hist(street[record.good==1, at.npp.ann.rel]) ## insanely skewed
## literal median tree
## first cut the tiniest trees
street[dbh.2006<5, record.good:=0]
street[record.good==1, median(dbh.2006, na.rm=T)] ## 23.1cm
street[record.good==1 & dbh.2006==median(dbh.2006, na.rm=T), median(at.npp.ann.rel, na.rm=T)] ## 28.7% annualized!
hist(street[record.good==1, at.npp.ann.rel]) ## in fucking sanely skewed
hist(street[record.good==1 & dbh.2006<60 & dbh.2006>20, at.npp.ann.rel])
street[record.good==1 & dbh.2006<60 & dbh.2006>20, median(at.npp.ann.rel, na.rm=T)]
length(street[record.good==1 & dbh.2006<60 & dbh.2006>20, at.npp.ann.rel]) #1219
length(street[record.good==1 & dbh.2006<60 & dbh.2006>0, at.npp.ann.rel]) #2045
length(street[record.good==1 & dbh.2006>20, at.npp.ann.rel]) #1324
length(street[record.good==1 & dbh.2006>60, at.npp.ann.rel]) # only 105 trees are above 60cm
length(street[record.good==1 & dbh.2006<20, at.npp.ann.rel]) # 826 are below 20cm
hist(street[record.good==1 & dbh.2006>20 & dbh.2006<60, at.npp.ann.rel]) #1324
street[record.good==1 & dbh.2006>20 & dbh.2006<60, median(at.npp.ann.rel, na.rm=T)] ## 5.96%
### OHHHH when we include all the tiny ones (and there are many) the median ann.npp goes up
street[record.good==1 & dbh.2006==median(dbh.2006, na.rm=T), median(at.npp.ann.rel, na.rm=T)] ## 28.7% annualized!
hist(street[record.good==1 & dbh.2006>21 & dbh.2006<25, at.npp.ann.rel]) ### skewed AF!


### OK analysis aside, let's do Approach 2a: standard trees
### Smith's dbh-->biomass equation
b0 <- -2.48
b1 <- 2.4835 ## these are eastern hardwood defaults
biom.pred <- function(x){exp(b0+(b1*log(x)))}

## so to get 1 year npp per tree now, you take the biomass.2006, multiply by at.npp.ann.rel*1 (1 year)
std.npp.rel <- median(street[record.good==1 & dbh.2006>21 & dbh.2006<25, at.npp.ann.rel], na.rm=T) ## 12.87% -- I can live with this
## first thing, figure out how many standard trees are in each cell
std.dbh.2006 <- street[record.good==1, median(dbh.2006, na.rm=T)] ## 23.1cm
## number of standard street trees per cell
biom.dat[bos.biom30m>50, num.std.trees:=bos.biom30m/biom.pred(std.dbh.2006)] ## at least some tree should be there
range(biom.dat[,num.std.trees], na.rm=T) ## up to 250 trees
#900/250 = 3.6m2/tree, or about a 1m radius around each (is this realistic?)
num.trees <- raster(biom)
num.trees <- setValues(num.trees, biom.dat$num.std.trees)
plot(num.trees, main="Number of street trees")
std.npp.rel*biom.pred(std.dbh.2006) ## each standard tree is growing by 26kgBiomass/yr
250*std.npp.rel*biom.pred(std.dbh.2006) # at peak density can get 6.6 Mgbiomass/cell/yr
75*std.npp.rel*biom.pred(std.dbh.2006) ## even at moderate you're seeing 2 Mgbiomass/cell/yr
biom.dat[num.std.trees>0, npp.std.streettree:=biom.pred(std.dbh.2006)*std.npp.rel*num.std.trees]
biom.dat[aoi>0 & is.na(npp.std.streettree), npp.std.streettree:=0]
biom.dat[, npp.std.streettree:=npp.std.streettree/1000/2] # convert kgBiomass to MgC 
npp.std.streettree <- raster(biom)
npp.std.streettree <- setValues(npp.std.streettree, biom.dat[, npp.std.streettree])
plot(npp.std.streettree, main="NPP, standard streettree")

### compare to FIA based
fia.ne <- raster("processed/boston/bos.npp.actual.fia.NEdefault.tif")
plot(fia.ne, main="NPP, FIA growth")
fia.ne.dat <- as.data.table(as.data.frame(fia.ne))
fia.ne.dat[,sum(bos.npp.actual.fia.NEdefault, na.rm=T)] ## 4794 MgC/yr
biom.dat[,sum(npp.std.streettree, na.rm=T)] ## 45,705 MgC/yr -- god damn

## approach 2b
# 1) distribute a floating number of trees with dbh/biomass drawn from empirical distribution, toss out combos that are not close enough to cell total biomass
# 2) apply npp.ann.rel~dbh.2006 all-genus relationship to determine biomass inrement in individual trees, sum per cell*tree population combo
target <- biom.dat[,median(bos.biom30m, na.rm=T)]
target <- 35000
biom.pred(std.dbh.2006)
a <- rnorm(n = 400, mean = target/biom.pred(std.dbh.2006)) ## e.g. a very dense pixel with 35k kg would have 171 stems at median dbh
hist(a)
dump <- numeric()
cage.dbh <- list()
x <- 0
start <- proc.time()
g <- 0
while(x<100){
  ## grab a random number of randomly selected trees
  n <- round(runif(1, min=0, max=80),0)
  grasp <- sample(street[record.good==1, at.biom.2006], size = n)
  test <- sum(grasp)
  if(test<(1.10*target) & test>(0.90*target)){ ## keep if they sum to close to the biomass in the cell
    dump <- c(dump, test)
    cage.dbh[[x+1]] <- grasp
    x <- x+1
    print(x)
  }
  g <- g+1
}
proc.time()-start
g
showme <- unlist(lapply(cage.dbh, length)) # number of trees in each test
hist(showme) ## ok this is a distribution of 0-14 trees from the empirical set that gives you biomass within +/- 10% of the Raciti biomass
median(showme)

17*biom.dat[!is.na(bos.biom30m), length(bos.biom30m)]/60/60 ### could take 651 hours
sqrt(10/pi)
## main problems: if we don't restrict the search space, for small biomass it takes forever for it to stumble upon a low-n sample
## figure how to narrow the search envelope re. number of trees based on biomass
## problem: the median number for correct samples seems to track biomass pretty cleanly, up into unreasonably high density
## basically we are independently sampling the entire distribution for every cell, so a lot of median trees are getting put down in the workable samples
## and basically the median n value that works is the number of median trees that can be crammed in to equal the biomass
## it's like a long-way round way of basically saying "put the std.median tree here, only also allow for some tails of weirder ones

## need some sort of independent constraint on the number of trees per pixel

## Andy's work implies that the maximum density they got near the edge was 3.6m2/ha of BA; 
## on a median street tree (~24cm dbh), this is 80 trees/pixel
## IN A FULLY CANOPIED PIXEL
## can constrain the search for biomass combinations that do not exceed
## critical stem density **in the canopied part of the pixel**


## approach 2b
# 1) distribute a floating number of trees with dbh/biomass drawn from empirical distribution, toss out combos that are not close enough to cell total biomass
# 2) apply npp.ann.rel~dbh.2006 all-genus relationship to determine biomass inrement in individual trees, sum per cell*tree population combo
# target <- 35000 ## example: a very dense patch with near perfect canopy
can <- raster("processed/boston/bos.can30m.tif")
can.dat <- as.data.table(as.data.frame(can))
biom.dat <- cbind(biom.dat, can.dat)

## standard street tree
target <- biom.dat[,median(bos.biom30m, na.rm=T)]
can.target <- can.dat[,median(bos.can30m, na.rm=T)]

## forest
target <- biom.dat[bos.biom30m<45000 & bos.biom30m>40000, median(bos.biom30m, na.rm=T)]
can.target <- biom.dat[bos.biom30m<45000 & bos.biom30m>40000,median(bos.can30m, na.rm=T)]

## low biomass
target <- biom.dat[bos.biom30m<400 & bos.biom30m>200, median(bos.biom30m, na.rm=T)]
can.target <- biom.dat[bos.biom30m<400 & bos.biom30m>200,median(bos.can30m, na.rm=T)]
## this one takes up to 5x longer, ~5s
## test flexible
hi.test <- 12000
lo.test <- 10000

target <- biom.dat[bos.biom30m<hi.test & bos.biom30m>lo.test, median(bos.biom30m, na.rm=T)]
can.target <- biom.dat[bos.biom30m<hi.test & bos.biom30m>lo.test,median(bos.can30m, na.rm=T)]

## set up a sample of targets to test
# set.seed(10)
biom.dat[bos.biom30m>0, .N] ## 107k pixels have nonzero biomass
biom.clean <- biom.dat[bos.biom30m>0 & bos.biom30m<12000,] ## restrict to less dense
biom.clean[,ind:=1:dim(biom.clean)[1]]
ind.samp <- sample(biom.clean[, ind], size = 100)
target <- biom.clean[ind.samp, bos.biom30m]
can.target <- biom.clean[ind.samp, bos.can30m]

### model for annual npp~dbh.2006
mod.biom.rel <- summary(lm(log(at.npp.ann.rel)~log(dbh.2006), data=street[record.good==1 & ann.npp!=0])) #r2=0.54! so relative %change biomass per year is the most predictive so far


cage.num.trees <- list()
cage.ann.npp <- list()
runtime <- numeric()
targs <- seq(1, dim(street[record.good==1])[1])


## initial idea: this gets hung badly on very large biomasses -- it's waiting for the
### rare event where it selects a lot of big trees
for(t in 1:100){
  start <- Sys.time()
  ann.npp <- numeric()
  num.trees <- numeric()
  x <- 0
  while(x<100){
    ## grab a random number of randomly selected trees
    ## this is probably not restricted enough, wasting a ton of time on samples that are going to be much too large to re biomass
    # n <- round(runif(1, min=1, max=(1800*can.target[t])),0) ## 1800 is max n per fully forested pixel at 5cm dbh
    n <- round(runif(1, min=1, max=60),0) ## gonna arbitrarily restrict this to hopefully produce higher # of workable matches faster
    dime <- sample(targs, size = n) ## grab n rows of street tree data
    grasp <- street[record.good==1, .(dbh.2006, at.biom.2006)][dime]
    if((sum(((grasp$dbh.2006/2)^2)*pi)/1E4)<(can.target[t]*3.6)){ ## if the BA density is low enough
      test <- sum(grasp$at.biom.2006)
      if(test<(1.10*target[t]) & test>(0.90*target[t])){ ## if biomass sum to close to the biomass in the cell
        ## apply all the npp annual relative factors to dbhs for this selection
        # hist(exp((mod.biom.rel$coefficients[2]*log(grasp$dbh.2006))+mod.biom.rel$coefficients[1]))
        ann.npp <- c(ann.npp, sum(biom.pred(grasp$dbh.2006)*exp((mod.biom.rel$coefficients[2]*log(grasp$dbh.2006))+mod.biom.rel$coefficients[1])))
        num.trees <- c(num.trees, length(grasp$dbh.2006))
        x <- x+1
        print(x)
      }
    }
  }
  cage.ann.npp[[t]] <- ann.npp
  cage.num.trees[[t]] <- num.trees
  runtime <- c(runtime, Sys.time()-start)
}

### optimize: restrict sample to data that is not single-tree bigger than the biomass present
### take a whole (large enough) sample of the combed data
### then sequentially add the biomasses and quit once you get to the threshold or exceed it (hit till you make it or bust)

## set up a sample of targets to test
# set.seed(10)
biom.dat[bos.biom30m>0 & !is.na(bos.can30m), .N] ## 107k pixels have nonzero biomass
biom.clean <- biom.dat[bos.biom30m>0 & bos.biom30m<12000,] ## restrict to less dense
biom.clean[,ind:=1:dim(biom.clean)[1]]
ind.samp <- sample(biom.clean[, ind], size = 100)
target <- biom.clean[ind.samp, bos.biom30m]
can.target <- biom.clean[ind.samp, bos.can30m]

### model for annual npp~dbh.2006
mod.biom.rel <- summary(lm(log(at.npp.ann.rel)~log(dbh.2006), data=street[record.good==1 & ann.npp!=0])) #r2=0.54! so relative %change biomass per year is the most predictive so far

cage.num.trees <- list()
cage.ann.npp <- list()
runtime <- numeric()
biom.track <- numeric()
targs <- seq(1, dim(street[record.good==1])[1])
clean <- street[record.good==1,]
setkey(clean, at.biom.2006)

## this loop works a hell of a lot faster than above, doesn't waste time testing samples that are too many trees and dynamically restricts sampling for low-mass pixels
## but it is incapable of handling large sizes -- it runs out of trees to pack in
for(t in 1:100){
  start <- Sys.time()
  ann.npp <- numeric()
  num.trees <- numeric()
  x <- 0
  q <- 0
  while(x<100 & q<4000){
    ## grab a random number of randomly selected trees
    ## this is probably not restricted enough, wasting a ton of time on samples that are going to be much too large to re biomass
    # n <- round(runif(1, min=1, max=(1800*can.target[t])),0) ## 1800 is max n per fully forested pixel at 5cm dbh
    # n <- round(runif(1, min=1, max=60),0) ## gonna arbitrarily restrict this to hopefully produce higher # of workable matches faster
    # dime <- sample(targs, size = n) ## grab n rows of street tree data
    # grasp <- street[record.good==1 & at.biom.2006<target[t], .(dbh.2006, at.biom.2006)][dime]
    grasp <- clean[at.biom.2006<target[t], .(dbh.2006, at.biom.2006)]
    n <- min(c(150, dim(grasp)[1])) ## if you get a really tiny biomass it might try to sample too many rows
    grasp <- grasp[sample(dim(grasp)[1], size=n),]
    w=grasp[1, at.biom.2006] ## keep cummulative tally of biomass
    d=1 # keep track of the number of trees
    while(w<(0.9*target[t])){ ## keep adding trees until you just get over the target biomass
      w=w+grasp[d+1, at.biom.2006]
      d=d+1
    }
    ### if you've gone too far or if you've packed them in too tight, ditch this sample
    if(grasp[1:d, sum((((dbh.2006/2)^2)*pi)/1E4)]<(can.target[t]*3.6) & w<(1.10*target[t])){ ## if the BA density is low enough & didn't overshoot biomass too much
        ann.npp <- c(ann.npp, sum(grasp[1:d, at.biom.2006]*exp((mod.biom.rel$coefficients[2]*log(grasp[1:d, dbh.2006]))+mod.biom.rel$coefficients[1])))
        num.trees <- c(num.trees, d)
        x <- x+1
        print(x)
    }
    q=q+1 ## if you can't find a combination that works, try q times to get develop a sample and if you can't fuck this pixel
  }
  if(q<4000){
    cage.ann.npp[[t]] <- ann.npp
    cage.num.trees[[t]] <- num.trees
    runtime <- c(runtime, Sys.time()-start)
    biom.track <- c(biom.track, target[t])
    print(paste("finished pixel", t))
  } else{
    cage.ann.npp[[t]] <- 9999 ## could not find a solution
    cage.num.trees[[t]] <- 9999
    runtime <- c(runtime, Sys.time()-start)
    biom.track <- c(biom.track, target[t])
    print(paste("finished pixel", t, "error"))
  }

}

hist(cage.ann.npp[[1]])
hist(cage.num.trees[[1]])
target[1]
hist(biom.track)



### version 3, second gear for very dense places
biom.dat[bos.biom30m>20000,length(bos.biom30m)]/biom.dat[,length(bos.biom30m)] ## only 7% of pixels are >10000, 2% above 20000, so say for that top 7% we can help it out

## e.g. take the top 20% of dbh and triple it, then draw from that larger sample

### set up a test dataset
clean <- street[record.good==1 & dbh.2006>=5,]
clean[,quantile(dbh.2006, probs= 0.80)] # 40.64
hist(clean[,dbh.2006])
clean <- rbind(clean, clean[dbh.2006>40,], clean[dbh.2006>40,], clean[dbh.2006>40,]) ## add the top 20% in another 3 times
hist(clean[,dbh.2006]) ## bumped up the moderately large more
clean[, quantile(dbh.2006, 0.8)] ## now upper 80% starts at 52
setkey(clean, at.biom.2006)

## pare down the biomass data to only v. large
bigtest <- biom.dat[bos.biom30m>10000,]
bigtest.grande <- bigtest[bos.biom30m>30000,]
biom.clean <- bigtest.grande
biom.clean <- bigtest
biom.clean[,ind:=1:dim(biom.clean)[1]]
ind.samp <- sample(biom.clean[, ind], size = 100)
target <- biom.clean[ind.samp, bos.biom30m]
can.target <- biom.clean[ind.samp, bos.can30m]

### model for annual npp~dbh.2006
mod.biom.rel <- summary(lm(log(at.npp.ann.rel)~log(dbh.2006), data=street[record.good==1 & ann.npp!=0])) #r2=0.54! so relative %change biomass per year is the most predictive so far

cage.num.trees <- list()
cage.ann.npp <- list()
cage.dbh <- list()
runtime <- numeric()
biom.track <- numeric()

for(t in 1:100){
  start <- Sys.time()
  ann.npp <- numeric()
  num.trees <- numeric()
  x <- 0
  q <- 0
  while(x<100 & q<4000){
    ## grab a random number of randomly selected trees
    ## this is probably not restricted enough, wasting a ton of time on samples that are going to be much too large to re biomass
    # n <- round(runif(1, min=1, max=(1800*can.target[t])),0) ## 1800 is max n per fully forested pixel at 5cm dbh
    # n <- round(runif(1, min=1, max=60),0) ## gonna arbitrarily restrict this to hopefully produce higher # of workable matches faster
    # dime <- sample(targs, size = n) ## grab n rows of street tree data
    # grasp <- street[record.good==1 & at.biom.2006<target[t], .(dbh.2006, at.biom.2006)][dime]
    grasp <- clean[at.biom.2006<target[t], .(dbh.2006, at.biom.2006)]
    n <- min(c(150, dim(grasp)[1])) 
    grasp <- grasp[sample(dim(grasp)[1], size=n),]
    w=grasp[1, at.biom.2006] ## keep cummulative tally of biomass
    d=1 # keep track of the number of trees
    while(w<(0.9*target[t])){ ## keep adding trees until you just get over the target biomass
      w=w+grasp[d+1, at.biom.2006]
      d=d+1
    }
    ### if you've gone too far or if you've packed them in too tight, ditch this sample
    if(grasp[1:d, sum((((dbh.2006/2)^2)*pi)/1E4)]<(can.target[t]*3.6) & w<(1.10*target[t])){ ## if the BA density is low enough & didn't overshoot biomass too much
      ann.npp <- c(ann.npp, sum(grasp[1:d, at.biom.2006]*exp((mod.biom.rel$coefficients[2]*log(grasp[1:d, dbh.2006]))+mod.biom.rel$coefficients[1])))
      num.trees <- c(num.trees, d)
      x <- x+1
      print(x)
    }
    q=q+1 ## if you can't find a combination that works, try q times to get develop a sample and if you can't fuck this pixel
  }
  if(q<4000){
    cage.ann.npp[[t]] <- ann.npp
    cage.num.trees[[t]] <- num.trees
    cage.dbh[[t]] <- grasp[1:d, dbh.2006]
    runtime <- c(runtime, Sys.time()-start)
    biom.track <- c(biom.track, target[t])
    print(paste("finished pixel", t))
  } else{
    cage.ann.npp[[t]] <- 9999 ## could not find a solution
    cage.num.trees[[t]] <- 9999
    cage.dbh[[t]] <- 9999
    runtime <- c(runtime, Sys.time()-start)
    biom.track <- c(biom.track, target[t])
    print(paste("finished pixel", t, "error"))
  }
  
}

hist(cage.ann.npp[[85]])
hist(cage.num.trees[[85]])
hist(cage.dbh[[85]])
hist(biom.track)
see <- data.frame(biom.track)
see$med.dbh <- lapply(cage.dbh, median)
see$med.num.trees <- lapply(cage.num.trees, median)
see$med.ann.npp <- lapply(cage.ann.npp, median)

View(see)
plot(see$biom.track[see$med.dbh<1000], see$med.dbh[see$med.dbh<1000])
plot(see$biom.track[see$med.dbh<1000], see$med.num.trees[see$med.dbh<1000])
plot(see$biom.track[see$med.dbh<1000], see$med.ann.npp[see$med.dbh<1000])
### let's just say we'll cut it off at 20000, and leave the few above this unknown
biom.dat[bos.biom30m>20000, length(bos.biom30m)]/biom.dat[!is.na(bos.biom30m), length(bos.biom30m)] ## 5.8% pixels above 20000 kg

################
#########
### final script for estimating tree distribution and running on the cluster
library(data.table)
library(raster)

### biomass equation biom~dbh
b0 <- -2.48
b1 <- 2.4835 ## these are eastern hardwood defaults
biom.pred <- function(x){exp(b0+(b1*log(x)))}

## get the street tree record prepped
street <- read.csv("docs/ian/Boston_Street_Trees.csv") ## Street tree resurvey, biomass 2006/2014, species
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
# a <- street[,length(ann.npp)/dim(street)[1], by=genus]
# a <- a[order(a$V1, decreasing = T),]
# # sum(a$V1[1:12]) ### top 12 is 94% of trees present
# focus <- a$genus[1:12]
street[, delta.dbh:=dbh.2014-dbh.2006]
street[, record.good:=0]
street[is.finite(dbh.2014) & is.finite(dbh.2006) & dead.by.2014==0 & delta.dbh>=0, record.good:=1] 
street[record.good==1, at.biom.2014:=biom.pred(dbh.2014)]
street[record.good==1, at.biom.2006:=biom.pred(dbh.2006)]
street[record.good==1, at.delta.biom:=at.biom.2014-at.biom.2006]
street[record.good==1, at.npp.ann:=at.delta.biom/8]
street[record.good==1, at.npp.ann.rel:=at.npp.ann/at.biom.2006]
street[record.good==1,range(dbh.2006, na.rm=T)]
street[dbh.2006<5, record.good:=0] ## filter for tiny trees, full record is now 2455

### biomass increment model based on dbh (exponential function)
mod.biom.rel <- summary(lm(log(at.npp.ann.rel)~log(dbh.2006), data=street[record.good==1 & ann.npp!=0])) #r2=0.54! so relative %change biomass per year is the most predictive so far
plot(street[record.good==1, log(at.biom.2006)], street[record.good==1, log(at.npp.ann.rel)])

## read in biomass and canopy data
biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m kg-biomass to 30m pixel
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- crop(biom, aoi) ## biomass was slightly buffered, need to clip to match canopy fraction raster
biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[,aoi:=as.vector(getValues(aoi))]
biom.dat[aoi<800, bos.biom30m:=NA] ### cancel values in cells that are not at least ~90% complete coverage in AOI
can <- raster("processed/boston/bos.can30m.tif")
can.dat <- as.data.table(as.data.frame(can))
biom.dat <- cbind(biom.dat, can.dat)
biom.dat[,index:=1:dim(biom.dat)[1]]

## prep street tree data and biomass data for processing (small biomass first)
clean <- street[record.good==1 & dbh.2006>=5,] # get a good street tree set ready
setkey(clean, at.biom.2006)

### first part: process the low biomass with distributions drawn from the unmodified street set
runme <- biom.dat[!is.na(bos.biom30m) & bos.biom30m<20000 & bos.biom30m>0,] #99k

## set up containers
cage.num.trees <- list()
cage.ann.npp <- list()
cage.dbh <- list()
index.track <- integer()
biom.track <- numeric()

## loop each row
for(t in 1:dim(runme)[1]){
  ann.npp <- numeric()
  num.trees <- numeric()
  x <- 0
  q <- 0
  while(x<100 & q<2000){ ## select 100 workable samples, or quit after 4000 attempts
    ## grab a random number of randomly selected trees
    grasp <- clean[at.biom.2006<runme[t, bos.biom30m], .(dbh.2006, at.biom.2006)]
    n <- min(c(80, dim(grasp)[1])) ## if you get a really tiny biomass it might try to sample too many rows
    grasp <- grasp[sample(dim(grasp)[1], size=n),]
    w=grasp[1, at.biom.2006] ## keep cummulative tally of biomass
    d=1 # keep track of the number of trees
    while(w<(0.9*runme[t, bos.biom30m])){ ## keep adding trees until you just get over the target biomass
      w=w+grasp[d+1, at.biom.2006]
      d=d+1
    }
    ### if you've gone too far or if you've packed them in too tight, ditch this sample
    if(grasp[1:d, sum((((dbh.2006/2)^2)*pi)/1E4)]<(runme[t, bos.can30m]*3.6) & w<(1.10*target[t])){ ## if the BA density is low enough & didn't overshoot biomass too much
      ann.npp <- c(ann.npp, sum(grasp[1:d, at.biom.2006]*exp((mod.biom.rel$coefficients[2]*log(grasp[1:d, dbh.2006]))+mod.biom.rel$coefficients[1])))
      num.trees <- c(num.trees, d)
      x <- x+1
    }
    q=q+1 ## if you can't find a combination that works, try q times to get develop a sample and if you can't fuck this pixel
  }
  if(q<2000){
    cage.ann.npp[[t]] <- ann.npp
    cage.num.trees[[t]] <- num.trees
    cage.dbh[[t]] <- grasp[1:d, dbh.2006]
    biom.track <- c(biom.track, runme[t, bos.biom30m])
    index.track <- c(index.track, runme[t, index])
    print(paste("finished pixel", t))
  } else{
    cage.ann.npp[[t]] <- 9999 ## could not find a solution
    cage.num.trees[[t]] <- 9999
    cage.dbh[[t]] <- 9999
    biom.track <- c(biom.track, runme[t, bos.biom30m])
    index.track <- c(index.track, runme[t, index])
    print(paste("pixel", t, "error"))
  }
}

# hist(unlist(cage.ann.npp[[4]]))
# hist(unlist(cage.dbh[[4]]))
# hist(unlist(cage.num.trees[[4]]))
# biom.track
## appears working
save(cage.ann.npp, file=paste("processed/boston/ann.npp.street.small"))
save(cage.num.trees, file=paste("processed/boston/num.trees.street.small"))
save(cage.dbh, file=paste("processed/boston/dbh.street.small"))
save(biom.track, file=paste("processed/boston/biom.track.street.small"))
save(cage.ann.npp, file=paste("processed/boston/index.track.street.small"))


####
### now specificially tackle the high-biomass areas
### modify the street tree record to amplify high-dbh trees
clean <- street[record.good==1 & dbh.2006>=5,]
bar <- clean[,quantile(dbh.2006, probs= 0.80)] ## what part of the high end to amplify
bar2 <- clean[,quantile(dbh.2006, probs= 0.90)]
clean <- rbind(clean, clean[dbh.2006>bar,], clean[dbh.2006>bar,], clean[dbh.2006>bar,], clean[dbh.2006>bar,], clean[dbh.2006>bar,]) ## add the top 20% in another 3 times
setkey(clean, at.biom.2006)

### first part: process the low biomass with distributions drawn from the unmodified street set
runme <- biom.dat[!is.na(bos.biom30m) & bos.biom30m>=20000 & bos.biom30m<30000,] #6k

## set up containers
cage.num.trees <- list()
cage.ann.npp <- list()
cage.dbh <- list()
index.track <- integer()
biom.track <- numeric()

## loop each row
for(t in 1:dim(runme)[1]){
  ann.npp <- numeric()
  num.trees <- numeric()
  x <- 0
  q <- 0
  while(x<100 & q<1000){ ## select 100 workable samples, or quit after 4000 attempts
    ## grab a random number of randomly selected trees
    grasp <- clean[, .(dbh.2006, at.biom.2006)]
    n <- min(c(80, dim(grasp)[1])) ## if you get a really tiny biomass it might try to sample too many rows
    grasp <- grasp[sample(dim(grasp)[1], size=n),]
    w=grasp[1, at.biom.2006] ## keep cummulative tally of biomass
    d=1 # keep track of the number of trees
    while(w<(0.9*runme[t, bos.biom30m])){ ## keep adding trees until you just get over the target biomass
      w=w+grasp[d+1, at.biom.2006]
      d=d+1
    }
    ### if you've gone too far or if you've packed them in too tight, ditch this sample
    if(grasp[1:d, sum((((dbh.2006/2)^2)*pi)/1E4)]<(runme[t, bos.can30m]*3.6) & w<(1.10*target[t])){ ## if the BA density is low enough & didn't overshoot biomass too much
      ann.npp <- c(ann.npp, sum(grasp[1:d, at.biom.2006]*exp((mod.biom.rel$coefficients[2]*log(grasp[1:d, dbh.2006]))+mod.biom.rel$coefficients[1])))
      num.trees <- c(num.trees, d)
      x <- x+1
    }
    q=q+1 ## if you can't find a combination that works, try q times to get develop a sample and if you can't fuck this pixel
  }
  if(q<1000){
    cage.ann.npp[[t]] <- ann.npp
    cage.num.trees[[t]] <- num.trees
    cage.dbh[[t]] <- grasp[1:d, dbh.2006]
    biom.track <- c(biom.track, runme[t, bos.biom30m])
    index.track <- c(index.track, runme[t, index])
    print(paste("finished pixel", t))
  } else{
    cage.ann.npp[[t]] <- 9999 ## could not find a solution
    cage.num.trees[[t]] <- 9999
    cage.dbh[[t]] <- 9999
    biom.track <- c(biom.track, runme[t, bos.biom30m])
    index.track <- c(index.track, runme[t, index])
    print(paste("pixel", t, "error"))
  }
}

hist(unlist(cage.ann.npp[[3]]))
hist(unlist(cage.dbh[[3]]))
hist(unlist(cage.num.trees[[3]]))
biom.track
## appears working
save(cage.ann.npp, file=paste("processed/boston/ann.npp.street.big"))
save(cage.num.trees, file=paste("processed/boston/num.trees.street.big"))
save(cage.dbh, file=paste("processed/boston/dbh.street.big"))
save(biom.track, file=paste("processed/boston/biom.track.street.big"))
save(cage.ann.npp, file=paste("processed/boston/index.track.street.big"))



################
### reconstitute the data from the pixel guesser run on the cluster
#### import results objects pulled from parallel processing on the cluster (chunks of 10k pixels)
### reconstitute the NPP estimates
### small biomass pixels first -- have consistent pix id scheme
obj.dump <- list.files("processed/boston/biom_street/")
npp.dump <- obj.dump[grep(obj.dump, pattern = "ann.npp.street.small*")]
npp.dump.chunks <- sub('.*small.', '', npp.dump)

## to get basal area
ba <- function(x){(((((x/2)^2)*pi)/900))} ## this is in m2/ha (correct for 900m2 pixel footprint)

### for each pixel, get median estimated npp, median # trees
container <- data.frame()
for(c in 1:length(npp.dump)){
  load(paste("processed/boston/biom_street/", npp.dump[c], sep=""))
  tmp.npp <- sapply(cage.ann.npp, FUN=median)
  tmp.npp[tmp.npp==9999] <- NA ## filter the NA flags
  
  ## grab the corresponding other dump files
  load(paste("processed/boston/biom_street/num.trees.street.small.", npp.dump.chunks[c], sep=""))
  tmp.num <- sapply(cage.num.trees, FUN=median)
  tmp.num[tmp.num==9999] <- NA ## filter the NA flags
  load(paste("processed/boston/biom_street/index.track.street.small.", npp.dump.chunks[c], sep="")) ## comes in as "index.track" object
  tmp.index <- index.track
  load(paste("processed/boston/biom_street/biom.track.street.small.", npp.dump.chunks[c], sep="")) ## comes in as "biom.track" object
  tmp.biom <- biom.track
  
  ### figure median dbh and median basal area for each pixel
  load(paste("processed/boston/biom_street/dbh.street.small.", npp.dump.chunks[c], sep="")) ## comes in as "index.track" object
  ba.track <- rep(9999, length(cage.dbh))
  dbh.track <- rep(9999, length(cage.dbh))
  print(paste("patience! Doing BA calculations on chunk", npp.dump.chunks[c]))
  for(b in 1:length(cage.dbh)){
    if(is.na(tmp.npp[b])){  ## if a valid NPP was retrieved only
      ba.track[b] <- NA
      dbh.track[b] <- NA
      } else{
      ba.track[b] <- median(sapply(sapply(cage.dbh[[b]], FUN=ba), FUN=sum)) ## median of the summed dbh-->ba values for each sample
      dbh.track[b] <- median(unlist(cage.dbh[[b]])) ## grand median dbh of all trees selected for all samples
      }
    if(b%%1000==0){print(b)}
  }
  
  ## bind to container
  g <- cbind(index.track, biom.track, tmp.npp, tmp.num, ba.track, dbh.track)
  container <- rbind(container, g)
  hist(tmp.npp, main=paste(npp.dump.chunks[c]))
  hist(tmp.num, main=paste(npp.dump.chunks[c]))
  hist(dbh.track, main=paste(npp.dump.chunks[c]))
  hist(ba.track, main=paste(npp.dump.chunks[c]))
}

## collect, export, make maps
names(container) <- c("id.proc1", "biom.kg", "ann.npp.street.sim", "tree.num.street.sim", 
                      "ba.street.sim", "med.dbh.street.sim")

# write.csv(container, "processed/boston/bos.street.trees.small.npp.simulatorv1.results.csv")
dim(container) # 98974 -- small only -- compare to 108974 when big trees included
sum(container$tree.num.street.sim, na.rm=T) #1.1M trees

### tedious reconstruction of the raster comenses here
### this is how the biomass raster was processed prior to the street tree simulator start
biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m kg-biomass to 30m pixel
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- crop(biom, aoi) ## biomass was slightly buffered, need to clip to match canopy fraction raster
biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[,aoi:=as.vector(getValues(aoi))]
biom.dat[aoi<800, bos.biom30m:=NA] ### cancel values in cells that are not at least ~90% complete coverage in AOI
can <- raster("processed/boston/bos.can30m.tif")
can.dat <- as.data.table(as.data.frame(can))
biom.dat <- cbind(biom.dat, can.dat)
biom.dat[,id.proc1:=1:dim(biom.dat)[1]] # 1:354068
biom.dat.proc <- biom.dat ## save a copy here

map <- merge(x=biom.dat.proc, y=container, by="id.proc1", all.x=T, all.y=T)
map[!is.na(bos.biom30m),] ## seems alright
plot(map[,bos.biom30m], map[,biom.kg]) ## lines up, but doesn't capture anything larger than 20k kg
range(map[,bos.biom30m], na.rm=T) ## original data still goes up to 50k kg


### big tree processing
### repeat processing for the ~6k "big biomass" pixels
obj.dump <- list.files("processed/boston/biom_street/")
npp.dump <- obj.dump[grep(obj.dump, pattern = "ann.npp.street.big*")]

container.big <- data.frame()
### for each pixel, get median estimated npp, median # trees
for(c in 1){
  load(paste("processed/boston/biom_street/", npp.dump[c], sep=""))
  tmp.npp <- sapply(cage.ann.npp, FUN=median)
  tmp.npp[tmp.npp==9999] <- NA ## filter the NA flags
  
  ## grab the corresponding other dump files
  load(paste("processed/boston/biom_street/num.trees.street.big"))
  tmp.num <- sapply(cage.num.trees, FUN=median)
  tmp.num[tmp.num==9999] <- NA ## filter the NA flags
  load(paste("processed/boston/biom_street/index.track.street.big")) ## comes in as "index.track" object
  tmp.index <- index.track
  load(paste("processed/boston/biom_street/biom.track.street.big")) ## comes in as "biom.track" object
  tmp.biom <- biom.track
  
  ### figure median dbh and median basal area for each pixel
  load(paste("processed/boston/biom_street/dbh.street.big")) ## comes in as "index.track" object
  ba.track <- rep(9999, length(cage.dbh))
  dbh.track <- rep(9999, length(cage.dbh))
  print(paste("patience! Doing BA calculations on chunk big"))
  for(b in 1:length(cage.dbh)){
    if(is.na(tmp.npp[b])){  ## if a valid NPP was retrieved only
      ba.track[b] <- NA
      dbh.track[b] <- NA
    } else{
      ba.track[b] <- median(sapply(sapply(cage.dbh[[b]], FUN=ba), FUN=sum)) ## median of the summed dbh-->ba values for each sample
      dbh.track[b] <- median(unlist(cage.dbh[[b]])) ## grand median dbh of all trees selected for all samples
    }
    if(b%%1000==0){print(b)}
  }
  
  ## bind to container
  g <- cbind(index.track, biom.track, tmp.npp, tmp.num, ba.track, dbh.track)
  container.big <- rbind(container.big, g)
  hist(tmp.npp, main=paste("big"))
  hist(tmp.num, main=paste("big"))
  hist(dbh.track, main=paste("big"))
  hist(ba.track, main=paste("big"))
}

## collect, export, make maps
names(container.big) <- c("id.proc1", "biom.kg", "ann.npp.street.sim", "tree.num.street.sim", 
                      "ba.street.sim", "med.dbh.street.sim")

### the runme file is altered for the big trees
runme <- biom.dat[!is.na(bos.biom30m) & !is.na(bos.can30m) & bos.biom30m>=20000 & bos.biom30m<30000,] #6k, ditch the small fraction that are too large to ever solve
dim(runme) # 6131
runme[,range(id.proc1)] # id.proc1 range 11385:350224
plot(runme[,id.proc1])
check <- as.vector(runme$id.proc1)
sum(check %in% container$id.proc1) ## none of these pixels appears in the list of small biom. pixels
range(container$id.proc1) # 1109:352514
dim(container) # 98974 small trees

range(container.big$id.proc1) # 251972:277668
plot(container.big$id.proc1) ## something is fucked about how index was assigned to the big processing
# runme.x gives proper index range in processing script
## saved the wrong object as index track, need to rerun

