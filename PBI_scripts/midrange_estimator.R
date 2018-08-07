##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Show that midrange estimator of the uniform distribution has a
##### standard error which drops as ~1/N. Cf. mean which drops as 1/sqrt(N).

library(gplots) # for plotCI
midrange <- function(x) {0.5*sum(range(x))}
set.seed(150)

nmax  <- 1e2
nsamp <- 1e3
sigma <- 2/sqrt(12) # standard deviation of U(-1,+1)
nVec  <- seq(from=2, to=nmax, by=2)
est <- matrix(NA, nrow=nmax, ncol=4) # mean(mu), mean(mr), sd(mean), sd(mr)
for(n in nVec) { # samp is a n*nsamp matrix
  samp <- Vectorize(runif, "n")(n=rep(n, nsamp), min=-1, max=1)
  mu <- apply(samp, 2, mean)     # vector size nsamp, mean estimator
  mr <- apply(samp, 2, midrange) # vector size nsamp, midrange estimator
  est[n,1] <- mean(mu)
  est[n,2] <- sd(mu)
  est[n,3] <- mean(mr)
  est[n,4] <- sd(mr)
}

pdf("midrange_estimator.pdf", 8, 4)
par(mfrow=c(1,2), mgp=c(2.0,0.8,0), mar=c(3.5,3.5,1,1), oma=0.1*c(1,1,1,1))
plot(nVec, est[nVec,1], type="n", ylim=c(-0.45, 0.45), yaxs="i", 
     xlab="N", ylab="mean")
lines(nVec,  sigma*1/sqrt(nVec), col="grey60", lw=2.5)
lines(nVec, -sigma*1/sqrt(nVec), col="grey60", lw=2.5)
plotCI(nVec, est[nVec,1], uiw=est[nVec,2], gap=0, cex=0.5, add=TRUE)
plot(nVec, est[nVec,3], type="n", ylim=c(-0.45, 0.45), yaxs="i", 
     xlab="N", ylab="midrange")
lines(nVec,  sigma*1/sqrt(nVec), col="grey60", lw=2)
lines(nVec, -sigma*1/sqrt(nVec), col="grey60", lw=2)
lines(nVec,  sigma*sqrt(6/((nVec+2)*(nVec+1))), col="grey60", lw=6)
lines(nVec, -sigma*sqrt(6/((nVec+2)*(nVec+1))), col="grey60", lw=6)
plotCI(nVec, est[nVec,3], uiw=est[nVec,4], gap=0, cex=0.5, add=TRUE)
dev.off()
