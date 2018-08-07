##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Bayesian inference of a 3-parameter linear model to 2D data

library(gplots) # for plotCI
source("metropolis.R")
source("linearmodel_functions.R") # provides logpost.linearmodel

########## Define true model and simulate experimental data from it

set.seed(50)
Ndat <- 10
x <- sort(runif(Ndat, 0, 10))
sigTrue <- 1
modMat <- c(0,1) # 1 x P vector: coefficients, b_p, of sum_{p=0} b_p*x^p
y <- cbind(1,x) %*% as.matrix(modMat) + rnorm(Ndat, 0, sigTrue)
# Dimensions in the above: [Ndat x 1] = [Ndat x P] %*% [P x 1] + [Ndat]
# cbind does the logical thing when combining a scalar and vector, 
# then do vector addition
y <- drop(y) # converts into a vector
pdf("linearmodel_data.pdf", width=4, height=4)
par(mfrow=c(1,1), mar=c(3.5,3.5,0.5,0.5), oma=0.1*c(1,1,1,1), 
    mgp=c(2.0,0.8,0), cex=1.0)
plot(x, y, xlim=c(0,9), ylim=c(-1,9), xaxs="i", yaxs="i")
#abline(a=modMat[1], b=modMat[2], col="red") # true model
dev.off()
# True parameters, transformed to conform with model to be used below
thetaTrue <- c(modMat[1], atan(modMat[2]), log10(sigTrue))
obsdata <- data.frame(cbind(x,y)) # columns must be named "x" and "y"
rm(x,y)

### Define model and infer the posterior PDF over its parameters

# Model to infer: linear regression with Gaussian noise
# Parameters: intercept b_0, gradient b_1; Gaussian noise sigma, ysig.
# MCMC works on: theta=c(b_0, alpha=tan(b_1), log10(ysig)), a 1x3 vector.
# Prior PDFs:
# b_0:         N(mean=m, sd=s); m,s estimated from global properties of data
# alpha:       Uniform (0 to 2pi)
# log10(ysig): Uniform (improper)

# Define covariance matrix of MCMC sampling PDF: 
# sigma=c(b_0, alpha, log10(ysig))
sampleCov <- diag(c(0.1, 0.02, 0.1)^2)
# Set starting point
thetaInit <- c(2, pi/8, log10(3))
# Run the MCMC to get samples from the posterior PDF
set.seed(150)
allSamp <- metrop(func=logpost.linearmodel, thetaInit=thetaInit, Nburnin=0, 
                  Nsamp=5e4, sampleCov=sampleCov, verbose=1e3, 
                  obsdata=obsdata)
# 10^(allSamp[,1]+allSamp[,2]) is the unnormalized posterior at each sample
thinSel  <- seq(from=1, to=nrow(allSamp), by=25) # thin by factor 25
postSamp <- allSamp[thinSel,]

# Plot MCMC chains and use density estimation to plot 1D posterior PDFs.
# We don't need to do any explicit marginalization to get the 1D PDFs.
pdf("linearmodel_mcmc.pdf", width=8, height=7)
par(mfrow=c(3,2), mar=c(3.0,3.5,0.5,0.5), oma=0.5*c(1,1,1,1), 
    mgp=c(1.8,0.6,0), cex=0.9)
parname <- c(expression(b[0]), expression(paste(alpha, " / rad")), 
             expression(paste(log, " ", sigma)))
for(j in 3:5) { # columns of postSamp
  plot(1:nrow(postSamp), postSamp[,j], type="l",
       xlab="iteration", ylab=parname[j-2])
  postDen <- density(postSamp[,j], n=2^10)
  plot(postDen$x, postDen$y, type="l", lwd=1.5, yaxs="i", 
       ylim=1.05*c(0,max(postDen$y)), xlab=parname[j-2], ylab="density")
  abline(v=thetaTrue[j-2], lwd=1.5, lty=3)
  if(j==3) { # overplot prior
    b0Val <- seq(from=-8, to=8, by=0.01)
    lines(b0Val, dnorm(b0Val, mean=0, sd=2), lty=2)
  }
}
dev.off()

# Plot gradient and intercept of samples in 2D.
# Fix range for b_0 vs alpha to enable comparison to centered data case
pdf("linearmodel_parameter_correlations.pdf", width=6, height=2)
par(mfrow=c(1,3), mgp=c(2.0,0.8,0), mar=c(3.0,3.0,0.5,0.5), 
    oma=0.1*c(1,1,1,1))
plot(postSamp[,3], postSamp[,4], xlab=parname[1], ylab=parname[2], 
     pch=".", xlim=c(-2,2), ylim=c(0.35,0.95), xaxs="i", yaxs="i")
plot(postSamp[,3], postSamp[,5], xlab=parname[1], ylab=parname[3], pch=".")
plot(postSamp[,4], postSamp[,5], xlab=parname[2], ylab=parname[3], pch=".")
dev.off()

# Find MAP and mean solutions.
# MAP is not the peak in each 1D PDF, but the peak of the 3D PDF.
# Mean is easy, as samples were drawn from the (unnormalized) posterior.
posMAP    <- which.max(postSamp[,1]+postSamp[,2]) 
thetaMAP  <- postSamp[posMAP, 3:5]
thetaMean <- apply(postSamp[,3:5], 2, mean) # Monte Carlo integration
cov(postSamp[, 3:5]) # covariance
cor(postSamp[, 3:5]) # correlation

# Overplot these solutions with original data
pdf("linearmodel_fits.pdf", width=4, height=4)
par(mfrow=c(1,1), mar=c(3.5,3.5,0.5,0.5), oma=0.1*c(1,1,1,1), 
    mgp=c(2.0,0.8,0), cex=1.0)
plotCI(obsdata$x, obsdata$y, xlim=c(0,9), ylim=c(-1,9), xaxs="i", yaxs="i",
       xlab="x", ylab="y", uiw=10^thetaMAP[3], gap=0)
abline(a=thetaMAP[1], b=tan(thetaMAP[2]), lw=2)     # MAP  model
#abline(a=modMat[1],   b=modMat[2], col="red", lw=2) # true model
# Compare this with the result from ML estimation from lm()
#abline(lm(obsdata$y ~ obsdata$x), col="black", lty=2)
dev.off()

### Make prediction: determine PDF(ycand | xnew, obsdata)

# Model and likelihood used here must be consistent with logpost.linearmodel

# Example 1
xnew <- 6
xlim <- c( 0,9) # xlim and ylim for plotting only
ylim <- c(-1,9)
# Example 2
#xnew <- 25
#xlim <- c( 0,27)
#ylim <- c(-1,29)

# Evaluate generative model at posterior samples (from MCMC).
# Dimensions in matrix multiplication: [Nsamp x 1] = [Nsamp x P] %*% [P x 1]
modPred <- cbind(postSamp[,3], tan(postSamp[,4])) %*% t(cbind(1,xnew))

# Direct method
# ycand must span full range of likelihood and posterior
dy    <- 0.01
ymid  <- thetaMAP[1] + xnew*tan(thetaMAP[2]) # to center choice of ycand
ycand <- seq(ymid-10, ymid+10, dy) # uniform grid of y with step size dy
ycandPDF <- vector(mode="numeric", length=length(ycand))
for(k in 1:length(ycand)) {
  like <- dnorm(ycand[k], mean=modPred, sd=10^postSamp[,5]) # [Nsamp x 1]
  ycandPDF[k] <- mean(like) # integration by rectangle rule. Gives a scalar
}
# Note that ycandPDF[k] is normalized, i.e. sum(dy*ycandPDF)=1.
# Find peak and approximate confidence intervals at 1sigma on either side
peak.ind  <- which.max(ycandPDF)
lower.ind <- max( which(cumsum(dy*ycandPDF) < pnorm(-1)) )
upper.ind <- min( which(cumsum(dy*ycandPDF) > pnorm(+1)) )
yPredDirect <- ycand[c(peak.ind, lower.ind, upper.ind)]

# Indirect method. likeSamp is [Nsamp x 1]
likeSamp <- rnorm(n=length(modPred), mean=modPred, sd=10^postSamp[,5])
likeDen  <- density(likeSamp, n=2^10)
# Find peak and confidence intervals
yPredIndirect <- c(likeDen$x[which.max(likeDen$y)], quantile(likeSamp, 
                          probs=c(pnorm(-1), pnorm(+1)), names=FALSE))

# Plot the predictive posterior distribution
pdf("linearmodel_prediction6_PDF.pdf", width=4, height=4)
par(mfrow=c(1,1), mar=c(3.0,3.5,0.5,0.5), oma=0.5*c(1,1,1,1), 
    mgp=c(2.2,0.8,0), cex=1.0)
plot(ycand, ycandPDF, type="l", lwd=1.5, yaxs="i", 
     ylim=1.05*c(0,max(ycandPDF)), xlab=expression(y[p]), 
     ylab=expression(paste("P(", y[p], " | ", x[p], ", D)"))) 
abline(v=yPredDirect, lty=2)
# overplot result from the indirect method
lines(likeDen$x, likeDen$y, type="l", lty=3, lwd=1.5)
dev.off()

# Compare predictions between the two methods
rbind(yPredDirect, yPredIndirect)

# Overplot direct prediction with original data and the MAP model
pdf("linearmodel_prediction6_ondata.pdf", width=4, height=4)
par(mfrow=c(1,1), mar=c(3.5,3.5,0.5,0.5), oma=0.1*c(1,1,1,1), 
    mgp=c(2.0,0.8,0), cex=1.0)
plotCI(obsdata$x, obsdata$y, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i", 
       uiw=10^thetaMAP[3], gap=0, xlab="x", ylab="y")
abline(a=thetaMAP[1], b=tan(thetaMAP[2]), lwd=2) # MAP  model
plotCI(xnew, ycand[peak.ind], li=ycand[lower.ind], ui=ycand[upper.ind],
  gap=0, add=TRUE, lwd=3)
dev.off()
