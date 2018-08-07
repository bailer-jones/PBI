##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Functions to evaluate the prior, likelihood, and posterior
##### for the linear model, plus to sample from the prior

# theta is vector of parameters
# obsdata is 2 column dataframe with names [x,y]
# The priors are hard-wired into the functions

# Return c(log10(prior), log10(likelihood)) (each generally unnormalized) 
# of the linear model
logpost.linearmodel <- function(theta, obsdata) {
  logprior <- logprior.linearmodel(theta)
  if(is.finite(logprior)) { # only evaluate model if parameters are sensible
    return( c(logprior, loglike.linearmodel(theta, obsdata)) ) 
  } else {
    return( c(-Inf, -Inf) )
  }
}

# Return log10(likelihood) (a scalar) for parameters theta and obsdata
# dnorm(..., log=TRUE) returns log base e, so multiply by 1/ln(10) = 0.434
# to get log base 10
loglike.linearmodel <- function(theta, obsdata) {
  # convert alpha to b_1 and log10(ysig) to ysig
  theta[2] <- tan(theta[2])
  theta[3] <- 10^theta[3]
  modPred <- drop( theta[1:2] %*% t(cbind(1,obsdata$x)) )
  # Dimensions in mixed vector/matrix products: [Ndat] = [P] %*% [P x Ndat] 
  logLike <- (1/log(10))*sum( dnorm(modPred - obsdata$y, mean=0,
                                    sd=theta[3], log=TRUE) )
  return(logLike)
}

# Return log10(unnormalized prior) (a scalar)
logprior.linearmodel <- function(theta) {
  b0Prior      <- dnorm(theta[1], mean=0, sd=2)
  alphaPrior   <- 1
  logysigPrior <- 1 
  logPrior <- sum( log10(b0Prior), log10(alphaPrior), log10(logysigPrior) )
  return(logPrior)
}
