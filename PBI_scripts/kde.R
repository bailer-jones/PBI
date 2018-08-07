##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Investigate kernel density estimation in 1D
##### Modified from the example given by Venables & Ripley (2002)

library(MASS) # for truehist and geyser data
attach(geyser)
pdf("kde.pdf", 8,4)
par(mfrow=c(1,2), mar=c(3.5,3.0,0.5,0.5), oma=c(0.5,0.5,0.5,0.5), 
    mgp=c(2.2,0.8,0), cex=1.0)
truehist(duration, h=0.25, xlim=c(0.5, 6), ymax=1.1, ylab="density", 
         col="white", lwd=0.5)
lines(density(duration, kernel="rectangular", bw=0.1, n=2^10), lwd=1.5)
lines(density(duration, kernel="rectangular", bw=0.2, n=2^10), lwd=2, 
      col="grey60")
truehist(duration, h=0.25, xlim=c(0.5, 6), ymax=1.1, ylab="density", 
         col="white", lwd=0.5)
lines(density(duration, kernel="gaussian",    bw=0.1, n=2^10), lwd=1.5)
lines(density(duration, kernel="gaussian",    bw=0.2, n=2^10), lwd=2, 
      col="grey60")
dev.off()
detach(geyser)
