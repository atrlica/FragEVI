urb1 <- c(rep(0.1, 10))
urb2 <- c(rep(0.1,5), rep(0.2,3), rep(0.5, 6))
veg <- rep(0.6, 10)

box <- list() 
box[[1]] <- c(0:100)
for(m in 0:100){
  pix <- vector()
  for(p in 1:500){
    urbfrac <- sample(urb1, m, replace=T)
    vegfrac <- sample(veg, 100-m, replace=T)
    pix[p] <- mean(c(urbfrac, vegfrac))
  }
  box[[m+2]] <- pix
  print(m)
}

res1 <- unlist(lapply(box, median))
plot(0:100, res[2:102], main="zero-veg urban", xlab="urban fraction", ylab="simulated EVI")


box <- list() 
for(m in 0:100){
  pix <- vector()
  for(p in 1:500){
    urbfrac <- sample(urb2, m, replace=T)
    vegfrac <- sample(veg, 100-m, replace=T)
    pix[p] <- mean(c(urbfrac, vegfrac))
  }
  box[[m+2]] <- pix
  print(m)
}

res2 <- unlist(lapply(box, median))
plot(0:100, res, main="urban w veg", xlab="urban fraction", ylab="simulated EVI")
