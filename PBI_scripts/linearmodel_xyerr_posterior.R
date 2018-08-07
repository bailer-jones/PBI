##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Bayesian inference of a 2-parameter linear model to 2D data with
##### uncertainties on x and y (one parameter for x, one for y)

library(gplots) # for plotCI
source("metropolis.R")
source("linearmodel_xyerr_functions.R") # provides logpost.linearmodel

########## Define true model, simulate experimental data from it, plot

set.seed(120)
Ndat <- 5
x <- sort(runif(Ndat, 0, 10))
modMat <- c(0,1) # 1 x P vector: coefficients, b_p, of polynomial sum_{p=0} b_p*x^p
y <- drop( cbind(1,x) %*% as.matrix(modMat) )
truedata <- data.frame(cbind(x=x,y=y)) # columns must be named "x" and "y"
thetaTrue <- c(modMat[1], atan(modMat[2]))
sigxTrue <- 3
sigyTrue <- 1
obsdata <- data.frame(cbind(x=x+rnorm(Ndat,0,sigxTrue), y=y+rnorm(Ndat,0,sigyTrue))) # columns must be named "x" and "y"
pdf("linearmodel_xyerr_data.pdf", width=4, height=4)
par(mfrow=c(1,1), mar=c(3.5,3.5,0.5,0.5), oma=0.1*c(1,1,1,1), mgp=c(2.0,0.8,0), cex=1.0)
plotCI(obsdata$x, obsdata$y, uiw=sigxTrue, err="x", gap=0, xlim=c(-1,11), ylim=c(-1,11), 
       xaxs="i", yaxs="i", xlab="x", ylab="y")
plotCI(obsdata$x, obsdata$y, uiw=sigyTrue, err="y", gap=0, add=TRUE)
abline(a=modMat[1], b=modMat[2], col="red", lw=2) # true model
points(truedata$x, truedata$y, pch=8)
dev.off()
# center the data
xMean <- mean(x)
yMean <- mean(y)
# the following are what we use in the model
obsdata$x <- obsdata$x - xMean
obsdata$y <- obsdata$y - yMean
sigx <- sigxTrue
sigy <- sigyTrue
rm(x,y)

### Define model and infer the posterior PDF over its parameters

# Model to infer: linear regression with Gaussian noise in x and y.
# Parameters: intercept b_0, gradient b_1
# MCMC works on: theta=c(b_0, alpha=tan(b_1), a 1x2 vector.
# Prior PDFs:
# b_0:         N(mean=m, sd=s); m,s estimated from global properties of data
# alpha:       Uniform (0 to 2pi)

# Define covariance matrix of MCMC sampling PDF: sigma=c(b_0, alpha)
sampleCov <- diag(c(0.25, 0.05)^2)
# Set starting point as least squares solution
coef <- lm(obsdata$y ~ obsdata$x)$coefficients
thetaInit <- c(coef[1], atan(coef[2]))
#thetaInit <- c(2, pi/8)
# Run the MCMC to get samples from the posterior PDF
set.seed(150)
allSamp <- metrop(func=logpost.linearmodel.xyerr, thetaInit=thetaInit, Nburnin=0, Nsamp=1e4,
                  sampleCov=sampleCov, verbose=10, obsdata=obsdata, sigx=sigx, sigy=sigy)
# 10^(allSamp[,1]+allSamp[,2]) is the unnormalized posterior at each sample
thinSel  <- seq(from=1, to=nrow(allSamp), by=25) # thin by factor 25
postSamp <- allSamp[thinSel,]
#postSamp <- allSamp

# Plot MCMC chains and use density estimation to plot 1D posterior PDFs from these.
# Note that we don't need to do any explicit marginalization to get the 1D PDFs.
pdf("linearmodel_xyerr_mcmc.pdf", width=8, height=7)
par(mfrow=c(2,2), mar=c(3.0,3.5,0.5,0.5), oma=0.5*c(1,1,1,1), mgp=c(1.8,0.6,0), cex=0.9)
parnames <- c(expression(b[0]), expression(paste(alpha, " / rad")))
for(j in 3:4) { # columns of postSamp
  plot(1:nrow(postSamp), postSamp[,j], type="l",
       xlab="iteration", ylab=parnames[j-2])
  postDen <- density(postSamp[,j], n=2^10)
  plot(postDen$x, postDen$y, type="l", lwd=1.5, yaxs="i", ylim=1.05*c(0,max(postDen$y)),
       xlab=parnames[j-2], ylab="density")
#  abline(v=thetaTrue[j-2], lwd=1.5, lty=3)
  if(j==3) { # overplot prior. Must manually set to agree with logprior.linearmodel 
    b0Val <- seq(from=-8, to=8, by=0.01)
    lines(b0Val, dnorm(b0Val, mean=0, sd=5), lty=2)
  }
}
dev.off()

# Plot gradient and intercept of samples in 2D.
# Fix range for b_0 vs alpha to enable comparison to centered data case
pdf("linearmodel_xyerr_parameter_correlations.pdf", width=4, height=4)
par(mfrow=c(1,1), mgp=c(2.0,0.8,0), mar=c(3.0,3.0,0.5,0.5), oma=0.1*c(1,1,1,1))
plot(postSamp[,3], postSamp[,4], xlab=parnames[1], ylab=parnames[2], pch=".",
     xaxs="i", yaxs="i")
dev.off()

# Find MAP and mean solutions.
# MAP is not the peak in each 1D PDF, but the peak of the 3D PDF.
# mean is easy, because samples have been drawn from the (unnormalized) posterior.
posMAP    <- which.max(postSamp[,1]+postSamp[,2]) 
thetaMAP  <- postSamp[posMAP, 3:4]
thetaMean <- apply(postSamp[,3:4], 2, mean) # Monte Carlo integration
cov(postSamp[, 3:4]) # covariance
cor(postSamp[, 3:4]) # correlation

# Overplot MAP solution with original data, in centered data space
# along with 20 draws from posterior, and least squares fit
pdf("linearmodel_xyerr_fits.pdf", width=4, height=4)
par(mfrow=c(1,1), mar=c(3.5,3.5,0.5,0.5), oma=0.1*c(1,1,1,1), mgp=c(2.0,0.8,0), cex=1.0)
plot(obsdata$x, obsdata$y, type="n", xlim=c(-9,9), ylim=c(-7,7), xaxs="i", yaxs="i", 
     xlab="x", ylab="y")
sel <- sample.int(n=nrow(postSamp), size=50)
for(j in sel) {
  abline(a=postSamp[j,3], b=tan(postSamp[j,4]), col="grey") # 50 posterior draws
}
abline(lm(obsdata$y ~ obsdata$x), col="black", lw=1.5, lty=2) # least squares model
abline(a=thetaMAP[1], b=tan(thetaMAP[2]), lw=2) # MAP  model
plotCI(obsdata$x, obsdata$y, uiw=sigxTrue, err="x", gap=0, add=TRUE)
plotCI(obsdata$x, obsdata$y, uiw=sigyTrue, err="y", gap=0, add=TRUE)
#abline(a=modMat[1] - yMean + xMean*modMat[2], b=modMat[2], col="red", lw=2) # true model
dev.off()

# As this can can take tens of minutes to run, save the entire session for later use
save(list=ls(all=TRUE), file="linearmodel_xyerr_posterior.Robj")
