##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Plot a Gaussian noise image

library(fields) # for image.plot
library(RColorBrewer) # for colorRampPalette
mypalette <- colorRampPalette(brewer.pal(9, "Greys"), space="rgb", 
                              interpolate="linear", bias=1)
mycols <- mypalette(64)
pdf("noisy_image.pdf", 5, 4)
par(mfrow=c(1,1), mgp=c(2.0,0.8,0), mar=c(1,1,1,1), oma=0.1*c(1,1,1,1))
set.seed(100)
x <- 1:1e2
y <- 1:1e2
z <- matrix(data=rnorm(length(x)*length(y))+10, nrow=length(y), 
            ncol=length(x))
image.plot(z=z, x=x, y=y, xaxt="n", xlab="", yaxt="n", ylab="", nlevel=1024, 
           zlim=c(-4.1,4.1)+10, col=mycols)
dev.off()
