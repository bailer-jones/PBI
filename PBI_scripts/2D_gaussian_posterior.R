##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Plot 2D posterior over Gaussian mu and sigma
##### for uniform prior on mu and Jeffreys prior on sigma

library(fields) # for image.plot
library(RColorBrewer) # for colorRampPalette
mypalette <- colorRampPalette(brewer.pal(9, "Greys"), space="rgb", 
                              interpolate="linear", bias=2.5)
mycols <- mypalette(64)

# Define function to return the unnormalized posterior
post <- function(mu, sigma, xbar, Vx, N) {
  (1/sigma^(N+1))*exp( (-N/(2*sigma^2)) * ((xbar-mu)^2 + Vx) )
}

# Define data and calculate posterior density on a dense grid
xbar  <- 0
Vx    <- 2^2
N     <- 10
mu    <- seq(from=-3,   to=3, length.out=1e3)
sigma <- seq(from=0.01, to=5, length.out=1e3)
postDen <- matrix(data=NA, nrow=length(mu), ncol=length(sigma))
for(i in 1:length(mu)) {
  for(j in 1:length(sigma)) {
    postDen[i,j] <- post(mu=mu[i], sigma=sigma[j], xbar=xbar, Vx=Vx, N=N)
  } 
}
postDen <- postDen/max(postDen) # scale so maximum is one

pdf("2D_gaussian_posterior.pdf", 5, 4)
par(mfrow=c(1,1), mar=c(3.5,3.5,0.5,1), oma=c(0.1,0.1,0.5,0.1), 
    mgp=c(2.2,0.8,0), cex=1.0) 
image.plot(z=postDen, x=mu, y=sigma, nlevel=1024, xlab=expression(mu), 
           ylab=expression(sigma), col=mycols, cex.lab=1.5)
dev.off()
