##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Functions to provide evaluations of prior, likelihood and posterior for the
##### linear model, plus sampling from the prior

# theta is vector of parameters; obsdata is 2 column dataframe with names [x,y].
# sigx and sigy are the fixed standard deviations in x and y.
# The priors are hard-wired into the functions.

# Return c(log10(prior), log10(likelihood)) (each generally unnormalized) of the linear model
logpost.linearmodel.xyerr <- function(theta, obsdata, sigx, sigy) {
  logprior <- logprior.linearmodel.xyerr(theta)
  if(is.finite(logprior)) { # only evaluate model if parameters are sensible
    return( c(logprior, loglike.linearmodel.xyerr(theta, obsdata, sigx, sigy)) ) 
  } else {
    return( c(-Inf, -Inf) )
  }
}

# Return log10(likelihood) for parameters theta for obsdata
# dnorm(..., log=TRUE) returns log base e, so multiply by 1/ln(10) = 0.4342945
# to get log base 10
loglike.linearmodel.xyerr <- function(theta, obsdata, sigx, sigy) {
  # convert alpha to b_1
  theta[2] <- tan(theta[2])
  logLpoint <- vector("logical", nrow(obsdata)) # loglike for each data point
  alt <- vector("logical", nrow(obsdata))
  for(i in 1:nrow(obsdata)) {
    xlim <- obsdata$x[i] + 5*sigx*c(-1,+1)
    #xlim <- c(-Inf, +Inf)
    integout <- integrate(f=Vectorize(likeIntegrand, "xprime"), lower=xlim[1], upper=xlim[2], 
                          subdivisions=1e2, rel.tol=1e-4, stop.on.error=FALSE, 
                          x=obsdata$x[i], y=obsdata$y[i], theta=theta, sigx=sigx, sigy=sigy)
    # I think this is normalized to within a factor which is independent of data or theta
    if(integout$message=="OK") {
      logLpoint[i] <- log10(integout$value)
      #cat(i, logLpoint[i], "\n")
    } else {
      stop("Numerical integration failed for point", i, "with theta=", theta, "\n")
    }
    alt[i] <- loglike.onepoint(x=obsdata$x[i], y=obsdata$y[i], sigx=sigx, sigy=sigy, 
                                     b0=theta[1], b1=theta[2])
  }
  logLike <- sum(logLpoint)
  #cat(logLike, sum(alt), "\n")
  return(logLike)
}

# Return log10(likelihood) of a single data point (x,y) 
# for model parameters b0, b1, sigx, sigy.
# Used for the analytic solution.
loglike.onepoint <- function(x, y, sigx, sigy, b0, b1) {
  c0 <-  x^2/sigx^2 + (y-b0)^2/sigy^2
  c1 <- -2*x/sigx^2 - 2*b1*(y-b0)/sigy^2
  c2 <-  1/sigx^2   + b1^2/sigy^2
  return( ((c1^2/(4*c2^2) - c0/2) - log(sqrt(pi)*c2*sigx*sigy))/log(10) )
}

# Return the integrand for the likelihood calculation at xprime,
# for a given data point (x,y) and model parameters theta, sigx, sigy.
# Used for the numerical solution.
likeIntegrand <- function(xprime, x, y, theta, sigx, sigy) {
  yprime <- drop( theta[1:2] %*% t(cbind(1,xprime)) ) # model prediction
  return( dmvnorm(x=c(x,y), mean=c(xprime,yprime), sigma=diag(c(sigx, sigy)^2)) )
}

# Return log10(unnormalized prior)
logprior.linearmodel.xyerr <- function(theta) {
  b0Prior      <- dnorm(theta[1], mean=0, sd=5)
  alphaPrior   <- 1
  logPrior <- sum( log10(b0Prior), log10(alphaPrior)) # scalar
  return(logPrior)
}
