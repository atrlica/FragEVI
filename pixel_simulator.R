urb1 <- c(rep(0.1, 10))
urb2 <- c(rep(0.1,5), rep(0.2,3), rep(0.5, 6))
veg <- rep(0.6, 10)

box <- list()
for(m in 0:100){ ## 101 bins from 0-100%, m is fraction urban
  a <- vector()
  for(j in 1:100){ ### number of pixels per bin
    pix <- vector()
    ufrac <- m/100*900
    vfrac <- 900-ufrac
    pix <- c(pix, sample(urb1, size = ufrac, replace=T))
    pix <- c(pix, sample(veg, size=vfrac, replace=T))
    a[j] <- mean(pix) # simulated landsat pixels EVI value
  }
  box[[m+1]] <- a
}

d <- unlist(lapply(box, FUN = median))
b <- seq(0, 100, by=1)
par(mfrow=c(1,2))
plot(b, d, col="gray30", main="No urban vegetation",
     ylab="Simulated EVI", xlab="% urban", pch=15)

### urban with some vegetation
box <- list()
for(m in 0:100){ ## 101 bins from 0-100%, m is fraction urban
  a <- vector()
  for(j in 1:100){ ### number of pixels per bin
    pix <- vector()
    ufrac <- m/100*900
    vfrac <- 900-ufrac
    pix <- c(pix, sample(urb2, size = ufrac, replace=T))
    pix <- c(pix, sample(veg, size=vfrac, replace=T))
    a[j] <- mean(pix) # simulated landsat pixels EVI value
  }
  box[[m+1]] <- a
}

d <- unlist(lapply(box, FUN = median))
b <- seq(0, 100, by=1)
plot(b, d, col="gray60", main="Random urban vegetation",
     ylab="Simulated EVI", xlab="% urban", pch=15,
     ylim=c(0.08, 0.62))
### so if urban contains patches of vegetation, the decline in EVI with urbanization is slower than you'd expect

## what if urban changes in veg fraction over its urbanization intensity?
## 3 urban classes: 0.1 (super developed), 0.35 (street tree), 0.52 (park/undev)
box <- list()
for(m in 0:100){ ## 101 bins from 0-100%, m is fraction urban
  a <- vector()
  urb3 <- rep(0.1, 10*(m+1)) ## start off with some fully barren
  urb3 <- c(urb3, rep(0.35, round((m+1)*5))) ## add in an increasing number
  if(m>30 & m<80){urb3 <- c(urb3, rep(0.5, 100-m))}
  urb3 <- sample(urb3, size=10)
  for(j in 1:1000){ ### number of pixels per bin
    pix <- vector()
    ufrac <- m/100*900
    vfrac <- 900-ufrac
    pix <- c(pix, sample(urb3, size = ufrac, replace=T))
    pix <- c(pix, sample(veg, size=vfrac, replace=T))
    a[j] <- mean(pix) # simulated landsat pixels EVI value
  }
  box[[m+1]] <- a
}

d <- unlist(lapply(box, FUN = median))
b <- seq(0, 100, by=1)
plot(b, d, col="green", main="Gradient urban vegetation",
     ylab="Simulated EVI", xlab="% urban", pch=15,
     ylim=c(0.08, 0.62))
