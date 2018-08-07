##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Calculate the Bayesian evidence for two linear models of 2D data.
##### M1: alpha=0, M2: P(alpha) ~ 1.

library(gplots) # for plotCI

### Define likelihood

# Return log10(likelihood), a scalar.
# theta is the vector of model parameters, here c(b_0, alpha, log10(ysig)).
# data is the two-column matrix [x,y].
# dnorm(..., log=TRUE) returns log base e, so multiply by 1/ln(10) = 0.434
# to get log base 10
log.like <- function(theta, data) {
  # convert alpha to b_1 and log10(ysig) to ysig
  theta[2] <- tan(theta[2])
  theta[3] <- 10^theta[3]
  # likelihood
  modPred <- drop( theta[1:2] %*% t(cbind(1,data$x)) )
  # Dimensions in mixed vector/matrix products: [Ndat] = [P] %*% [P x Ndat] 
  logLike <- (1/log(10))*sum( dnorm(modPred - data$y, mean=0, 
                                    sd=theta[3], log=TRUE) )
  return(logLike)
}

### Define true model and simulate experimental data from it

set.seed(75) # 75 gives data and plots in script
Ndat <- 10
x <- sort(runif(Ndat, -10, 10))
sigTrue <- 1
modMat <- c(0,0.1) # 1 x P vector: coefficients, b_p, of sum_{p=0} b_p*x^p
y <- cbind(1,x) %*% as.matrix(modMat) + rnorm(Ndat, 0, sigTrue) # noisy data
# Dimensions in matrix multiplication: 
# [Ndat x 1] = [Ndat x P] %*% [P x 1] + [Ndat]
# cbind does logical thing combining scalar and vector; then vector addition
y <- drop(y) # converts into a vector
pdf("linearmodel_evidence_data_10.pdf", width=5, height=4)
par(mfrow=c(1,1), mar=c(3.5,3.0,0.5,0.5), oma=0.5*c(1,1,1,1), 
    mgp=c(2.2,0.8,0), cex=1.2)
plot(x, y, xlim=c(-10,10))
#plotCI(x, y, xlim=c(-10, 10), uiw=sigTrue, gap=0) # data and true error bar
#abline(a=modMat[1], b=modMat[2], col="red") # true model
dev.off()

### Sample from prior

# Sample from prior.
# priorSamp is an array with dimensions (Nsamp, 3) containing the 
# samples for b_0, alpha, log10(sigma)
set.seed(100)
Nsamp <- 1e5 # will need to be larger if the priors are broader
priorSamp <- cbind(rnorm(n=Nsamp, mean=0, sd=1), 
                   runif(n=Nsamp, min=-pi/2, max=pi/2),
                   runif(n=Nsamp, min=log10(0.5), max=log10(2)))
sel <- sample.int(n=Nsamp, size=100) # 100 of the prior samples 

# Plot data and overplot 100 prior models from M2 
pdf("linearmodel_evidence_prior_models_M2.pdf", width=5, height=4)
par(mfrow=c(1,1), mar=c(3.5,3.0,0.5,0.5), oma=0.5*c(1,1,1,1), 
    mgp=c(2.2,0.8,0), cex=1.2)
plot(x, y, type="n", xlim=c(-10,10), ylim=c(-10,10), xaxs="i", yaxs="i")
for(j in sel) {
  abline(a=priorSamp[j,1], b=tan(priorSamp[j,2]), col="grey")
}
points(x, y)
dev.off()

# Plot data and overplot prior models from M1
pdf("linearmodel_evidence_prior_models_M1.pdf", width=5, height=4)
par(mfrow=c(1,1), mar=c(3.5,3.0,0.5,0.5), oma=0.5*c(1,1,1,1), 
    mgp=c(2.2,0.8,0), cex=1.2)
plot(x, y, type="n", xlim=c(-10,10), ylim=c(-10,10), xaxs="i", yaxs="i")
for(j in sel) {
  abline(a=priorSamp[j,1], b=0, col="grey")
}
points(x, y)
dev.off()

### Calculate likelihoods, evidences, and BF

data <- data.frame(cbind(x,y))

# Calculate likelihoods and evidence for model M2
logLikeM2 <- rep(NA, Nsamp)
for(j in 1:Nsamp) {
  logLikeM2[j] <- log.like(theta=priorSamp[j,], data)
}
logEvM2 <- log10(mean(10^logLikeM2))
# Calculate likelihoods and evidence for model M1
priorM1   <- cbind(priorSamp[,1], rep(0, Nsamp), priorSamp[,3])
logLikeM1 <- rep(NA, Nsamp)
for(j in 1:Nsamp) {
  logLikeM1[j] <- log.like(theta=priorM1[j,], data)
}
logEvM1 <- log10(mean(10^logLikeM1))
# Print results
cat("log10(Ev_1)  = ", logEvM1, "\n")
cat("log10(Ev_2)  = ", logEvM2, "\n")
cat("log10(BF_12) = ", logEvM1 - logEvM2, "\n")
cat("BF_12        = ", 10^(logEvM1 - logEvM2), "\n")
