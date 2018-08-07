##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Investigate histograms
##### Modified from the example given by Venables & Ripley (2002)

library(MASS) # for truehist and geyser data
attach(geyser)
# Look at sensitivity to bin width
pdf("histograms1.pdf", 12, 8)
par(mfrow=c(2,3), mar=c(3.5,3.5,1.5,0.5), oma=c(0.5,0.5,0.5,0.5), 
    mgp=c(2.2,0.8,0), cex=1.0)
for(h in c(2, 1, 1/2, 1/4, 1/8, 1/16)) {
 truehist(duration, h=h, x0=0.0, xlim=c(0,6), ymax=0.8, ylab="density", 
          col="white")
}
dev.off()
# Look at sensitivity to placing of bin centres
pdf("histograms2.pdf", 12, 8)
par(mfrow=c(2,3), mar=c(3.5,3.5,1.5,0.5), oma=c(0.5,0.5,0.5,0.5), 
    mgp=c(2.2,0.8,0), cex=1.0)
binwidth <- 0.5
for(x0 in binwidth*(0:5)/5) {
  truehist(duration, h=binwidth, x0=x0, xlim=c(0,6), ymax=0.8,
           ylab="density", col="white")
}
dev.off()
detach(geyser)
