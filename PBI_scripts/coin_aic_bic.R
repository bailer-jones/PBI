##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Demonstration of AIC and BIC on coin modelling comparison problem

phi <- 1/2
n   <- c(5,10,25,50,100)

pdf("coin_aic.pdf", 4, 4)
par(mfrow=c(1,1), mgp=c(2.0,0.8,0), mar=c(3.5,3.5,1,1), oma=0.1*c(1,1,1,1))
plot(c(0,1), c(0,1), type="n", xlim=c(0,1), ylim=c(-11,3), xlab="f=r/n", 
     ylab=expression(paste("AIC(",M[2],") - AIC(",M[1],")")))
for(i in 1:length(n)) {
  r  <- 0:n[i]
  f  <- r/n[i]
  aicM1 <- -2*dbinom(x=r, size=n[i], prob=1/2,    log=TRUE)
  aicM2 <- -2*dbinom(x=r, size=n[i], prob=r/n[i], log=TRUE) + 2
  points(f, aicM2-aicM1, cex=0.5, pch=13+i)
  lines(f,  aicM2-aicM1, lty=6-i)
}
dev.off()

pdf("coin_bic.pdf", 4, 4)
par(mfrow=c(1,1), mgp=c(2.0,0.8,0), mar=c(3.5,3.5,1,1), oma=0.1*c(1,1,1,1))
plot(c(0,1), c(0,1), type="n", xlim=c(0,1), ylim=c(-9,5), xlab="f=r/n", 
     ylab=expression(paste("BIC(",M[2],") - BIC(",M[1],")")))
for(i in 1:length(n)) {
  r  <- 0:n[i]
  f  <- r/n[i]
  bicM1 <- -2*dbinom(x=r, size=n[i], prob=1/2,    log=TRUE) 
  bicM2 <- -2*dbinom(x=r, size=n[i], prob=r/n[i], log=TRUE) + log(n[i])
  points(f, bicM2-bicM1, cex=0.5, pch=13+i)
  lines(f,  bicM2-bicM1, lty=6-i)
}
dev.off()
