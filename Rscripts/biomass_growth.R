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
biomass <- exp(-2.48+2.4835*log(post.dbh)) # the biomass of a median tree is 341kg

## determine the number of standard street trees present in each cell
biom <- raster("processed/boston/bos.biom30m.tif") ## this is summed 1m kg-biomass to 30m pixel
aoi <- raster("processed/boston/bos.aoi30m.tif")
biom <- crop(biom, aoi) ## biomass was slightly buffered, need to clip to match canopy fraction raster
biom.dat <- as.data.table(as.data.frame(biom))
biom.dat[,aoi:=as.vector(getValues(aoi))]
biom.dat[aoi<800, bos.biom30m:=NA] ### cancel values in cells that are not at least ~90% complete coverage in AOI

### live MgC by area in each pixel
biom.dat[, live.MgC.ha:=(bos.biom30m/aoi)*(1/2)*(1/1000)*(10^4)] ## convert kg-biomass/pixel based on size of pixel (summed by 1m2 in "aoi"), kgC:kgbiomass, Mg:kg, m2:ha

