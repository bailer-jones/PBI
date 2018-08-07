##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Compare data drawn from a Poisson with its theoretical distribution

pdf("dpois2.pdf", 4, 4)
par(mfrow=c(1,1), mgp=c(2.0,0.8,0), mar=c(3.5,3.5,1,1), oma=0.1*c(1,1,1,1))
truelambda <- 5 # = 1/tau
nint <- 100     # number of time intervals
set.seed(200)
ndecay <- rpois(n=nint, lambda=truelambda) # no.decays in each time interval
nobs   <- table(ndecay) # frequency distribution of ndecay
x <- 0:max(ndecay)
# multiply Poisson density by nint to get expected counts
plot(x, dpois(lambda=mean(ndecay), x=x)*nint, xlim=range(x), 
     ylim=c(0,max(nobs)), xlab="no. decays per time interval", 
     ylab="no. observations") 
points(nobs)
dev.off()
