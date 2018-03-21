clasSector <- rep(c("S01","S02","S03","S04","S05","S06","S07"),times=7)
Year <- rep(c("1950","1960","1970","1980","1990","2000","2010"),each=7)
Value <- runif(49, 10, 100)
df <- data.frame(Sector,Year,Value)

gg <- ggplot(df, aes(x=as.numeric(as.character(Year)), y=Value))
gg <- gg + geom_area(aes(colour=Sector, fill=Sector))
gg

evi.bin2 <- as.data.frame(evi.bin)


d <- reshape(evi.bin2, varying=c("m.forest", "m.lowveg", "m.lowres", "m.medres", "m.hdres", "m.dev"),
        idvar="bin", direction="long")
class(d$time)
d$time <- as.factor(d$time)

gg <- ggplot(d, aes(x=m.isa, y=m))
gg <- gg + geom_area(aes(colour=time, fill=time))

dat[, bin:=findInterval(isa, beta.range, all.inside=T)]
evi.bin <- dat[, .(m.isa=median(isa, na.rm=T), 
                             m.evi=median(evi, na.rm=T),
                             m.forest=median(forest, na.rm=T),
                             m.dev=median(dev, na.rm=T),
                             m.hdres=median(hdres, na.rm=T),
                             m.medres=median(medres, na.rm=T),
                             m.lowres=median(lowres, na.rm=T),
                             m.lowveg=median(lowveg, na.rm=T),
                   m.water=median(water, na.rm=T),
                   m.other=median(other, na.rm=T),
                             bin.count=.N), by=bin]
evi.bin <- evi.bin[order(m.isa),]
tot <- evi.bin[,sum(bin.count)]
evi.bin[,bin.frac:=bin.count/tot, by=bin]
evi.bin <- evi.bin[!is.na(m.isa),]


d <- reshape(evi.bin, varying=c("m.forest", "m.lowveg", "m.lowres", "m.medres", "m.hdres", "m.dev", "m.water", "m.other"),
             idvar="bin", direction="long")
class(d$time)
d$time <- as.factor(d$time)

gg <- ggplot(d, aes(x=m.isa, y=m))
gg <- gg + geom_area(aes(colour=time, fill=time))

