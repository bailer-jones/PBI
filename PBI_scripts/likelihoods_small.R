##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Illustration with the binomial distribution that likelihoods are small

pdf("likelihoods_small.pdf", 8, 4)
par(mfrow=c(1,2), mgp=c(2.2,0.8,0), mar=c(3.5,3.5,0.5,1),oma=0.5*c(1,1,1,0))
expVal <- c(10, 50, 100, 150)
plot(0:1, 0:1, type="n", xlim=c(0.7, 1.4), ylim=c(-8, 0), xlab="r / E[r]", 
     ylab=expression(paste(log," P(r | p=0.5, n=2*E[r])")))
for(i in 1:length(expVal)) {
  rVec <- seq(from=0.7*expVal[i], to=1.3*expVal[i], by=1)  
  loglike <- (1/log(10))*dbinom(rVec, size=2*expVal[i], prob=0.5, log=TRUE)
  points(rVec/expVal[i], loglike, pch=20, cex=1/i)
}
text(1.31, -0.3, "E[r]", pos=4)
text(1.31, c(-1.18, -3.13, -5.25, -7.32), expVal, pos=4)
expVal  <- 2^(1:10)
loglike <- (1/log(10))*dbinom(expVal, size=2*expVal, prob=0.5, log=TRUE)
plot(expVal, loglike, log="x", pch=20, xlab="E[r]", 
     ylab=expression(paste(log," P(r=E[r] | p=0.5, n=2*E[r])")))
dev.off()
