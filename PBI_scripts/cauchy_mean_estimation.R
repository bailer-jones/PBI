##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Is the mean a convergent estimator for the location of a Cauchy?

# Calculate mean and CLT prediction of its standard deviation, of a sample 
# of numbers drawn from (i) Cauchy, (ii) Gaussian, both with mode=0, HWHM=1, 
# with 2^lognsamp[i] no. of samples. Plot these against lognsamp[i]

set.seed(100)
lognsamp <- 1:20
GaussMean  <- vector(length=length(lognsamp))
CauchyMean <- vector(length=length(lognsamp))
GaussSD    <- vector(length=length(lognsamp))
CauchySD   <- vector(length=length(lognsamp))
for(i in 1:length(lognsamp)) {
  s <- rnorm(2^lognsamp[i], mean=0, sd=1/(2*log(2))) # sd=HWHM/2ln2=0.721
  GaussMean[i]  <- mean(s)
  GaussSD[i]    <- sd(s)/sqrt(length(s)) # standard deviation in mean
  s <- rcauchy(2^lognsamp[i], location=0, scale=1) # scale=HWHM
  CauchyMean[i] <- mean(s)
  CauchySD[i]   <- sd(s)/sqrt(length(s)) # standard deviation in mean
}

pdf("cauchy_mean_estimation.pdf", 8, 4)
par(mfrow=c(1,2), mgp=c(2.0,0.8,0), mar=c(3.5,3.5,1,1), oma=0.1*c(1,1,1,1))
plot(  lognsamp, GaussMean, ylim=range(c(GaussMean, CauchyMean)), type="n",
  xlab=expression(paste(log[2], N)), ylab="mean")
lines( lognsamp, GaussMean,  lwd=1.5)
points(lognsamp, GaussMean)
lines( lognsamp, CauchyMean, lwd=1.5, lty=2)
points(lognsamp, CauchyMean, pch=4)
plot(  lognsamp, GaussMean, ylim=range(c(GaussSD, CauchySD)), log="y", 
       type="n", xlab=expression(paste(log[2], N)), 
       ylab="standard deviation of mean")
lines( lognsamp, GaussSD,  lwd=1.5)
points(lognsamp, GaussSD)
lines( lognsamp, CauchySD, lwd=1.5, lty=2)
points(lognsamp, CauchySD, pch=4)
dev.off()
cbind(lognsamp, GaussSD, CauchySD)
