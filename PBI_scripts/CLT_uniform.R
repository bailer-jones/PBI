##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Apply central limit theorem to draws from uniform distribution

pdf("CLT_uniform.pdf", 12, 4)
par(mfrow=c(1,3), mgp=c(2,0.8,0), mar=c(3.5,3.5,1,0), oma=0.1*c(1,1,5,5), 
    cex=1.2)
set.seed(200)
for(N in c(1,2,5)) {
  Nsamp <- 10000
  z <- numeric(Nsamp)
  for (i in 1:Nsamp){
    z[i] <- mean(runif(N)*10)
  }
  binwidth <- 0.2
  hist(z, xlim=c(0,10), ylim=c(0,640), breaks=seq(0,10,binwidth), 
       ylab="frequency", main=paste("N =",N))
  cat("mean, sd = ", mean(z), sd(z), "\n")
  # overplot Gaussian with mean and sd from CLT
  gsd <- (10/sqrt(12))*1/sqrt(N)
  x <- seq(0,10,0.01)
  y <- dnorm(x, mean=5, sd=gsd)*Nsamp*binwidth
  lines(x, y, lwd=2)
}
dev.off()
