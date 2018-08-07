##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Compute posterior PDF for coin problem with a beta prior for a
##### range of r

n <- 20
alpha.prior <- 10
beta.prior  <- 10
Nsamp <- 200 # no. of points to sample at
pdf("coin2.pdf", 9, 7)
par(mfrow=c(3,4), mgp=c(2,0.8,0), mar=c(3.5,3.5,1.5,1), oma=0.5*c(1,1,1,1))
deltap <- 1/Nsamp # width of rectangles used for numerical integration
p <- seq(from=1/(2*Nsamp), by=1/Nsamp, length.out=Nsamp) # rectangle centres
for(r in seq(from=0, to=20, by=2)) {
  pdense <- dbeta(x=p, shape1=alpha.prior+r, shape2=beta.prior+n-r)
  plot(p, pdense, type="l", lwd=1.5, xlim=c(0,1), ylim=c(0, 6.5), 
       xaxs="i", yaxs="i", xlab="p", ylab="P(p | r,n,M)")
  title(main=paste("r =",r), line=0.3, cex.main=1.2)
  p.mean <- deltap*sum(p*pdense)
  abline(v=p.mean, lty=2)
  # overplot posterior obtained from a uniform prior
  pdense.uniform <- dbinom(x=r, size=n, prob=p)
  lines(p, pdense.uniform/(deltap*sum(pdense.uniform)), lwd=2, 
        col="grey60")
  # Can verify that pdense can also be found by direct calculation
  #pdense2 <- dbinom(x=r, size=n, prob=p) * 
  #           dbeta(x=p, shape1=alpha.prior, shape2=beta.prior)
  #pdense2 <- pdense2/(deltap*sum(pdense2)) # normalize posterior
  #lines(p, pdense2, col="red", lty=2)
}
dev.off()
