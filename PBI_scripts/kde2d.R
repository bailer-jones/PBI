##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Investigate kernel density estimation in 2D
##### Modified from the example given by Venables & Ripley (2002)

library(fields) # for image.plot
library(MASS)
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(9, "Greys"), space="rgb", 
                              interpolate="linear", bias=2)
mycols <- mypalette(64)
geyser2 <- data.frame(as.data.frame(geyser)[-1, ], 
                      pduration=geyser$duration[-299])
attach(geyser2)
pdf("kde2d.pdf", 12, 4)
par(mfrow=c(1,3), mar=c(3.5,3.5,0.5,1), oma=0.5*c(1,1,1,1), 
    mgp=c(2.2,0.8,0), cex=1.0) 
plot(pduration, waiting, xlim=c(0.5, 6), ylim=c(35, 115), xaxs="i", 
     yaxs="i", xlab="previous duration", ylab="waiting")
f1 <- kde2d(pduration, waiting, n=50, h=c(0.75, 10), 
            lims=c(0.5, 6, 35, 115))
image.plot(f1, zlim=c(0, 0.075), xlim=c(0.5, 6), ylim=c(35, 115), xaxs="i",
           yaxs="i", xlab="previous duration", ylab="waiting", col=mycols)
persp(f1, phi=30, theta=20, d=5, xlim=c(0.5, 6), ylim=c(35, 115), xaxs="i",
      yaxs="i", xlab="previous duration", ylab="waiting", zlab="density")
dev.off()
detach(geyser2)
