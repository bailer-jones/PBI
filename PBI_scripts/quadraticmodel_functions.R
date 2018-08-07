##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Functions to provide evaluations of prior, likelihood and posterior for the
##### quadratic model, plus sampling from the prior

# theta is vector of parameters; obsdata is 2 column dataframe with names [x,y].
# The priors are hard-wired into the functions.

# Return c(log10(prior), log10(likelihood)) (each generally unnormalized) of the quadratic model
logpost.quadraticmodel <- function(theta, obsdata) {
  logprior <- logprior.quadraticmodel(theta)
  if(is.finite(logprior)) { # only evaluate model if parameters are sensible
    return( c(logprior, loglike.quadraticmodel(theta, obsdata)) )
  } else {
    return( c(-Inf, -Inf) )
  }
}

# Return log10(likelihood) for parameters theta and obsdata
# dnorm(..., log=TRUE) returns log base e, so multiply by 1/ln(10) = 0.4342945
# to get log base 10
loglike.quadraticmodel <- function(theta, obsdata) {
  # convert alpha to b_1 and log10(ysig) to ysig
  theta[2] <- tan(theta[2])
  theta[4] <- 10^theta[4]
  modPred <- drop( theta[1:3] %*% t(cbind(1,obsdata$x,obsdata$x^2)) )
  # Dimensions in above mixed vector/matrix multiplication: [Ndat] = [P] %*% [P x Ndat] 
  logLike <- (1/log(10))*sum( dnorm(modPred - obsdata$y, mean=0, sd=theta[4], log=TRUE) )
  return(logLike)
}

# Return log10(unnormalized prior)
logprior.quadraticmodel <- function(theta) {
  b0Prior      <- dnorm(theta[1], mean=0, sd=10)
  alphaPrior   <- 1
  b2Prior      <- dnorm(theta[3], mean=0, sd=5)
  logysigPrior <- 1 
  logPrior <- sum( log10(b0Prior), log10(alphaPrior), log10(b2Prior), log10(logysigPrior) )
  return(logPrior)
}
