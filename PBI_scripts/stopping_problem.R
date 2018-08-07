##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### The stopping problem

p <- 0.1 
n <- 102 # no. trials
r <- 5   # no. successes

PA <- 2*sum(dbinom(x=0:r, size=n, prob=p))
# = 2*pbinom(q=r, size=n, prob=p)
PB <- 2*(1 - sum(dnbinom(x=0:(n-r-1), size=r, prob=0.1)))
# = 2*sum(dnbinom(x=(n-r):1e4, size=r, prob=0.1))
# where 1e4 is used as an approximation of infinity.
# Note the definitions adopted in dnbinom():
# x=n-r (no. failures) and size=r. Can check by comparing:
# choose(n-1, r-1)*p^r*(1-p)^(n-r)
# dnbinom(x=n-r, size=r, prob=p)

pdf("stopping_problem.pdf", 4, 4)
par(mfrow=c(1,1), mar=c(3.5,3.5,0.5,1), oma=0.5*c(1,1,1,1), 
    mgp=c(2.2,0.8,0))
Nsamp  <- 1e4
deltap <- 1/Nsamp
pgrid  <- seq(from=1/(2*Nsamp), by=1/Nsamp, length.out=Nsamp)
pdense <- pgrid^r * (1-pgrid)^(n-r)   # with uniform prior
pdense <- pdense/(deltap*sum(pdense)) # normalize posterior
# could instead do analytically for a uniform prior: pdense/beta(r+1, n-r+1)
plot(pgrid, pdense, type="l", lwd=1.5, xaxs="i", yaxs="i", xlim=c(0,0.2), 
     ylim=c(0,19), xlab="p", ylab="P(p | D,M')")
dev.off()
