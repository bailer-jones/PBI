##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Plot prior, likelihood, and posterior PDF for coin problem with a
##### beta prior for a range of r with n fixed

n <- 20
alpha.prior <- 10
beta.prior  <- 10
Nsamp <- 200 # no. of points to sample at
pdf("coin3.pdf", 9, 7)
par(mfrow=c(3,4), mgp=c(2,0.8,0), mar=c(3.5,3.5,1.5,1), oma=0.5*c(1,1,1,1))
deltap <- 1/Nsamp # width of rectangles used for numerical integration
p <- seq(from=1/(2*Nsamp), by=1/Nsamp, length.out=Nsamp) # rectangle centres
prior <- dbeta(x=p, shape1=alpha.prior, shape2=beta.prior)
for(r in seq(from=0, to=20, by=2)) {
  like  <- dbinom(x=r, size=n, prob=p)
  like  <- like/(deltap*sum(like)) # for plotting convenience only
  post  <- dbeta(x=p, shape1=alpha.prior+r, shape2=beta.prior+n-r)
  plot(p, prior, type="l", lwd=1.5, lty=2, xlim=c(0,1), ylim=c(0, 6.5),
       xaxs="i", yaxs="i", xlab="p", ylab="density")
  lines(p, like, lwd=1.5, lty=3)
  lines(p, post, lwd=1.5)
  title(main=paste("r =",r), line=0.3, cex.main=1.2)
}
dev.off()
