library(data.table)
library(reticulate)

hf.dbh <- read.csv("docs/ian/Harvard_Forest_biomass.csv")
length(unique(hf.dbh$Spp)) ## 18 species ID'd

andy.dbh <- read.csv("docs/ian/Reinmann_Hutyra_2016_DBH.csv")
length(unique(andy.dbh$Species_Code)) ## 17 species ID'd

andy.bai <- read.csv("docs/ian/Reinmann_Hutyra_2016_BAI.csv")
dim(andy.bai) ## 195 individuals, dbh by year -- convert to biomass by year?

rope <- read.csv("docs/ian/Rope_Plots_Data_Entry.csv") ## probably not useful, just a description of 
street <- read.csv("docs/ian/Boston_Street_Trees.csv") ## biomass 2006/2014, species
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
summary(lm(plotme$ann.npp~plotme$dbh.2014))
for(g in 1:length(focus)){
  # plot(plotme[genus==focus[g], dbh.2014], plotme[genus==focus[g], ann.npp], 
  #      main=focus[g])
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



