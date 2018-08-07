##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Bayesian inference of a 4-parameter quadratic model to 2D data

library(gplots) # for plotCI
source("metropolis.R")
source("quadraticmodel_functions.R") # provides logpost.quadraticmodel

### Define true model and simulate experimental data from it

set.seed(57)
Ndat <- 20
xrange <- c(0,10)
x <- sort(runif(n=Ndat, min=xrange[1], max=xrange[2]))
sigTrue <- 2
modMat <- c(25, -10, 1) # 1 x P vector: coefficients, b_p, of polynomial sum_{p=0} b_p*x^p
y <- cbind(1,x,x^2) %*% as.matrix(modMat) + rnorm(Ndat, 0, sigTrue)
# Dimensions in matrix multiplication: [Ndat x 1] = [Ndat x P] %*% [P x 1] + [Ndat]
# cbind does the logical thing combining a scalar and vector; then vector addition
y <- drop(y) # converts into a vector
pdf("quadraticmodel_data.pdf", width=4, height=4)
par(mfrow=c(1,1), mar=c(3.5,3.5,0.5,1), oma=0.1*c(1,1,1,1), mgp=c(2.0,0.8,0), cex=1.0)
plot(x, y, xlim=xrange, ylim=c(-6,29), xaxs="i", yaxs="i")
#xsamp <- seq(from=xrange[1], to=xrange[2], length.out=500)
#ysamp <- cbind(1,xsamp,xsamp^2) %*% as.matrix(modMat)
#lines(xsamp, drop(ysamp), col="red", lw=2) # true model
dev.off()
# True parameters, transformed to be conformable with model to be used below
thetaTrue <- c(modMat[1], atan(modMat[2]), modMat[3], log10(sigTrue))
obsdata <- data.frame(cbind(x,y)) # columns must be named "x" and "y"
rm(x,y)

### Define model and infer the posterior PDF over its parameters

# Model to infer: linear regression with Gaussian noise
# Parameters: intercept b_0, gradient b_1, quadratic term b_2; Gaussian noise sigma, ysig.
# MCMC works on: theta=c(b_0, alpha=tan(b_1), b_2, log10(ysig)), a 1x4 vector.
# Prior PDFs:
# b_0:         N(mean=m, sd=s); m,s estimated from global properties of data
# alpha:       Uniform (0 to 2pi)
# b_2:         N(mean=m, sd=s); m,s estimated from global properties of data
# log10(ysig): Uniform (improper)

# Define covariance matrix of MCMC sampling PDF: sigma=c(b_0, alpha, b_2, log10(ysig))
sampleCov <- diag(c(0.1, 0.01, 0.01, 0.01)^2)
# thetaInit <- c(0, 0, 0, log10(1)) # far from a good model
# Need fewer samples if initialize better: do a least squares fit and
# set thetaInit to: lsfit$coefficients, sqrt(mean(lsfit$residuals^2))
# lsfit <- lm(y ~ x + I(x^2), data=obsdata)
thetaInit <- c(27.4, atan(-11.7), 1.18, log10(2.4))
# Run the MCMC to find postSamp, samples of the posterior PDF
set.seed(250)
allSamp <- metrop(func=logpost.quadraticmodel, thetaInit=thetaInit, Nburnin=2e4, Nsamp=2e5,
                   sampleCov=sampleCov, verbose=1e3, obsdata=obsdata)
# 10^(allSamp[,1]+allSamp[,2]) is the unnormalized posterior at each sample
thinSel  <- seq(from=1, to=nrow(allSamp), by=100) # thin by factor 100
postSamp <- allSamp[thinSel,]

# Plot MCMC chains and use density estimation to plot 1D posterior PDFs from these.
# Note that we don't need to do any explicit marginalization to get the 1D PDFs.
pdf("quadraticmodel_mcmc.pdf", width=8, height=9.33)
par(mfrow=c(4,2), mar=c(3.0,3.5,0.5,0.5), oma=0.5*c(1,1,1,1), mgp=c(1.8,0.6,0), cex=0.9)
parnames <- c(expression(b[0]), expression(paste(alpha, " / rad")), expression(b[2]), 
              expression(paste(log, " ", sigma)))
for(j in 3:6) { # columns of postSamp
  plot(1:nrow(postSamp), postSamp[,j], type="l", xlab="iteration", ylab=parnames[j-2])
  postDen <- density(postSamp[,j], n=2^10)
  plot(postDen$x, postDen$y, type="l", lwd=1.5, yaxs="i", ylim=1.05*c(0,max(postDen$y)),
       xlab=parnames[j-2], ylab="density")
  abline(v=thetaTrue[j-2], lwd=1.5, lty=3)
}
dev.off()

# Plot all parameter samples in 2D
pdf("quadraticmodel_parameter_correlations.pdf", width=6, height=6)
par(mfcol=c(3,3), mar=c(3.5,3.5,0.5,0.5), oma=c(0.1,0.1,0.1,0.5), mgp=c(2.0,0.8,0))
for(i in 1:3) {
  for(j in 2:4) {
    if(j<=i) {
        plot.new()
      } else {
        plot(postSamp[,i+2], postSamp[,j+2], xlab=parnames[i], ylab=parnames[j], pch=".")
    }
  }
}
dev.off()

# Find MAP and mean solutions.
# MAP = Maximum A Posteriori, i.e. peak of posterior.
# MAP is not the peak in each 1D PDF, but the peak of the 4D PDF.
# mean is easy, because samples have been drawn from the (unnormalized) posterior.
posMAP    <- which.max(postSamp[,1]+postSamp[,2]) 
thetaMAP  <- postSamp[posMAP, 3:6]
thetaMean <- apply(postSamp[,3:6], 2, mean) # Monte Carlo integration

# Overplot MAP solution with original data
pdf("quadraticmodel_fits.pdf", width=4, height=4)
par(mfrow=c(1,1), mar=c(3.5,3.5,0.5,1), oma=0.1*c(1,1,1,1), mgp=c(2.0,0.8,0), cex=1.0)
plotCI(obsdata$x, obsdata$y, xlim=xrange, ylim=c(-6,29), xaxs="i", yaxs="i", 
       xlab="x", ylab="y", uiw=10^thetaMAP[4], gap=0)
xsamp <- seq(from=xrange[1], to=xrange[2], length.out=500)
ysamp <- cbind(1,xsamp,xsamp^2) %*% as.matrix(modMat)
#lines(xsamp, drop(ysamp), col="red", lwd=2) # true model
#ysamp <- cbind(1,xsamp,xsamp^2) %*% as.matrix(c(thetaMean[1], tan(thetaMean[2]), thetaMean[3])) 
#lines(xsamp, drop(ysamp), col="green", lwd=2) # mean model
ysamp <- cbind(1,xsamp,xsamp^2) %*% as.matrix(c(thetaMAP[1], tan(thetaMAP[2]), thetaMAP[3]))
lines(xsamp, drop(ysamp), lwd=2) # MAP model
dev.off()
