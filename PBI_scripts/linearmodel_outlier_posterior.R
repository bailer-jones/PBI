##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Bayesian inference of a 4-parameter linear model to 2D data,
##### which uses a mixture model for the likelihood to include an outlier model.

library(gplots) # for plotCI
source("metropolis.R")

# Model to infer: linear regression with Gaussian noise + Cauchy noise as a mixture model.
# Parameters: intercept b_0, gradient b_1; Gaussian noise sigma, ysig; mixture fraction, p.
# MCMC works on: theta=c(b_0, alpha=tan(b_1), log10(ysig), p), a 1x4 vector.
# Prior PDFs:
# b_0:         N(mean=m, sd=s), m,s estimated from global properties of data
# alpha:       Uniform (0 to 2pi)
# log10(ysig): Uniform (improper)
# p:           Beta(shape1, shape2); shape1,shape2 set manually

### Define posterior

# Calculate log10(likelihood) for each of the two likelihood components, and log10(prior)
#   for given theta and obsdata.
# If retmode=1 (default), return two-element vector (log10(prior), log10(likelihood).
# If retmode=2 return three-element vector LikeMain, LikeOutlier (not log)
#   where the two likelihood components are not multiplied by p or 1-p.
# theta is vector of parameters; obsdata is 2 column matrix with names [x,y] 
# (works with only one row too).
# The model and the priors and hard-wired into this function.
logpost.linoutmod <- function(theta, obsdata, retmode=1) {
  # convert alpha to b_1 and log10(ysig) to ysig
  theta[2] <- tan(theta[2])
  theta[3] <- 10^theta[3]
  # likelihood
  modPred <- drop( theta[1:2] %*% t(cbind(1,obsdata$x)) )
  # Dimensions in above mixed vector/matrix multiplication: [Ndat] = [P] %*% [P x Ndat] 
  likeMain    <- dnorm(modPred - obsdata$y,   mean=0, sd=theta[3]) # scalar
  likeOutlier <- dcauchy(modPred - obsdata$y, location=0, scale=1) # scalar
  if(theta[4]<0 | theta[4]>1) { # MCMC sampling can provide this
    logLike <- -Inf
  } else {
    logLike <- sum(log10( (1-theta[4])*likeMain + theta[4]*likeOutlier )) # scalar
  }
  # prior (not normalized - doesn't need to be, as we're returning unnormalized posterior)
  b0Prior      <- dnorm(theta[1], mean=0, sd=2)
  alphaPrior   <- 1
  logysigPrior <- 1
  pPrior       <- dbeta(theta[4], shape1=1, shape2=20)
  logPrior <- sum( log10(b0Prior), log10(alphaPrior), log10(logysigPrior), log10(pPrior) ) 
  if(retmode==2) {
    return( c(likeMain, likeOutlier) )
  } else {
    return( c(logPrior, logLike) )
  }
}

### Define basic model, simulate experimental data from it, add an outlier

set.seed(18)
Ndat <- 10
x <- sort(runif(Ndat, 0, 10))
sigTrue <- 1
modMat <- c(0,1) # 1 x P vector: coefficients, b_p, of polynomial sum_{p=0} b_p*x^p
y <- cbind(1,x) %*% as.matrix(modMat) + rnorm(Ndat, 0, sigTrue)
# Dimensions in above matrix multiplication: [Ndat x 1] = [Ndat x P] %*% [P x 1] + [Ndat]
# cbind does the logical thing combining a scalar and vector; then vector addition
#y[2] <- y[2] + 10*sigTrue # add an outlier
y[5] <- y[5] + 10*sigTrue # add an outlier
y <- drop(y) # converts into a vector
#plotCI(x, y, xlim=c(0,10), uiw=sigTrue, gap=0)
#abline(a=modMat[1], b=modMat[2], col="red") # true model
thetaTrue <- c(modMat[1], atan(modMat[2]), log10(sigTrue), NA)
# Above is the true parameters, transformed to be conformable with model to be used below
obsdata <- data.frame(cbind(x,y)) # columns must be named "x" and "y"
rm(x,y)

### Define model and infer the posterior PDF over its parameters

# define covariance matrix of MCMC sampling PDF: sigma=c(b_0, alpha, log10(ysig), p)
sampleCov <- diag(c(0.1, 0.01, 0.05, 0.01)^2)
# set starting point
thetaInit <- c(0, pi/4, log10(1), 0.1)
# run the MCMC to find postSamp, samples of the posterior PDF
set.seed(500)
allSamp <- metrop(func=logpost.linoutmod, thetaInit=thetaInit, Nburnin=1e4, Nsamp=1e5,
                   sampleCov=sampleCov, verbose=1e3, obsdata=obsdata)
# 10^(allSamp[,1]+allSamp[,2]) is the unnormalized posterior at each sample
thinSel  <- seq(from=1, to=nrow(allSamp), by=50) # thin by factor 50
postSamp <- allSamp[thinSel,]

# Plot MCMC chains and use density estimation to plot 1D posterior PDFs from these.
# Note that we don't need to do any explicit marginalization to get the 1D PDFs.
pdf("linearmodel_outlier_mcmc.pdf", 7, 7)
par(mfrow=c(4,2), mar=c(3.0,3.5,0.5,0.5), oma=0.5*c(1,1,1,1), mgp=c(1.8,0.6,0), cex=0.9)
parnames <- c(expression(b[0]), expression(paste(alpha, " / rad")), 
              expression(paste(log, " ", sigma)), "p")
for(j in 3:6) { # columns of postSamp
  plot(1:nrow(postSamp), postSamp[,j], type="l", xlab="iteration", ylab=parnames[j-2])
  postDen <- density(postSamp[,j], n=2^10)
  plot(postDen$x, postDen$y, type="l", lwd=1.5, yaxs="i", ylim=1.05*c(0,max(postDen$y)),
       xlab=parnames[j-2], ylab="density")
  abline(v=thetaTrue[j-2], lwd=1.5, lty=3)
  if(j==3) { # overplot prior
    b0Val <- seq(from=-8, to=8, by=0.01)
    lines(b0Val, dnorm(b0Val, mean=0, sd=2), lty=2)
  }
  if(j==6) { # overplot prior
    pVal <- seq(from=0, to=1, by=0.001)
    lines(pVal, dbeta(pVal, shape1=1, shape2=20), lty=2)
  }
}
dev.off()

# Plot all parameter samples in 2D
pdf("linearmodel_outlier_parameter_correlations.pdf", width=6, height=6)
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

# Find MAP solution and mean solution.
# MAP = Maximum A Posteriori, i.e. peak of posterior.
# MAP is not the peak in each 1D PDF, but the peak of the 3D PDF.
# mean is easy, because samples have been drawn from the (unnormalized) posterior.
posMAP     <- which.max(postSamp[,1]+postSamp[,2]) 
thetaMAP  <- postSamp[posMAP, 3:6]
thetaMean <- apply(postSamp[,3:6], 2, mean) # Monte Carlo integration

# Overplot MAP solution with original data
pdf("linearmodel_outlier_fits.pdf", 4, 4)
par(mfrow=c(1,1), mar=c(3.5,3.5,0.5,1), oma=0.1*c(1,1,1,1), mgp=c(2.0,0.8,0), cex=1.0)
plotCI(obsdata$x, obsdata$y, xlim=c(0,10), ylim=c(-2,17), xaxs="i", yaxs="i",
       xlab="x", ylab="y", uiw=10^thetaMAP[3], gap=0)
#abline(a=modMat[1],    b=modMat[2],         col="red",   lw=2) # true model
#abline(a=thetaMean[1], b=tan(thetaMean[2]), col="green", lw=2) # mean model
abline(a=thetaMAP[1],  b=tan(thetaMAP[2]), lw=2) # MAP  model
# Compare this with the result from ML estimation from lm()
abline(lm(obsdata$y ~ obsdata$x), col="black", lty=2)
dev.off()

# Calculate probability of each point being an outlier under the MAP solution
for(n in 1:Ndat) {
  LikeComp    <- logpost.linoutmod(thetaMAP, obsdata[n,], retmode=2)
  probOutlier <- thetaMAP[4]*LikeComp[2]/(thetaMAP[4]*LikeComp[2] + (1-thetaMAP[4])*LikeComp[1])
  cat(formatC(x=c(n, as.numeric(obsdata[n,]), probOutlier), digits=3, width=7),  "\n")
}
