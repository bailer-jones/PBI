##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Compute the posterior PDF for coin problem with a uniform prior for a
##### range of r

n <- 20
Nsamp <- 200 # no. of points to sample at
pdf("coin1.pdf", 9, 7)
par(mfrow=c(3,4), mgp=c(2,0.8,0), mar=c(3.5,3.5,1.5,1), oma=0.5*c(1,1,1,1))
deltap <- 1/Nsamp # width of rectangles used for numerical integration
p <- seq(from=1/(2*Nsamp), by=1/Nsamp, length.out=Nsamp) # rectangle centres
for(r in seq(from=0, to=20, by=2)) {
  pdense <- dbinom(x=r, size=n, prob=p)
  pdense <- pdense/(deltap*sum(pdense)) # normalize posterior
  plot(p, pdense, type="l", lwd=1.5, xlim=c(0,1), ylim=c(0,1.1*max(pdense)), 
       xaxs="i", yaxs="i", xlab="p", ylab="P(p | r,n,M)")
  title(main=paste("r =",r), line=0.3, cex.main=1.2)
  p.mean <- deltap*sum(p*pdense)
  abline(v=p.mean, lty=2)
}
dev.off()
