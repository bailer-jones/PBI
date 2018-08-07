##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Demonstration of bootstrapping and gradient ascent

### Define functions

# Prior (normalized)
d.prior3 <- function(r, rlen) ifelse(r>0, (1/(2*rlen^3))*r^2*exp(-r/rlen), 0) 

# Likelihood of one data point
d.like <- function(w, r, wsd) dnorm(x=w, mean=1/r, sd=wsd)

# Posterior (unnormalized) of set of data
d.post <- function(data, r, wsd, rlen) {
  prod(d.like(w=data, r, wsd))*d.prior3(r, rlen)
}

# Gradient of natural logarithm of d.post w.r.t r
grad.log.d.post <- function(r, d.post, data, wsd, rlen) {
  N <- length(data)
  dlnPdr <- 2/r - 1/rlen - (mean(data) - 1/r)*N/(wsd^2*r^2)
  return(dlnPdr)
  #return(dlnPdr * d.post(data, r, wsd, rlen))
}

# Calculate posterior on dense grid, scaled so mode=1.
# Return dataframe of r, posterior
post.grid <- function(data, wsd, rlen) {
  r <- seq(from=0, to=2*rlen, length.out=1e4)  
  post <- numeric(length(r))
  for(i in 1:length(r)) {
    post[i] <- d.post(data=data, r=r[i], wsd, rlen)
  }
  return(data.frame(r=r, post=post/max(post)))
}

# Draw bootstrap sample
draw.sample <- function(data) sample(x=data, size=length(data), 
                                     replace=TRUE)

# Gradient ascent (generic)
# Returns two element list: 
# - par is optimized value of parameter
# - lastIt is iteration number reached
# gradFunc is a function which returns the gradient of the function to 
# maximize. Its first argument must be the (scalar) value of the parameter.
# param is the intitial value of the parameter.
# '...' is used to pass anything else gradFunc needs (data etc).
# Search until relative change in gradient in less than rel.tol
# or maxIts iterations is reached.
grad.ascent <- function(gradFunc, param, ...) {
  rel.tol <- 1e-5
  alpha   <- 1e3
  maxIts  <- 1e4
  for(i in 1:maxIts) {
    paramNext <- param + alpha*gradFunc(param, ...)
    #cat(i, param, paramNext, "\n")
    if(abs(paramNext/param - 1) < rel.tol) {
      break
    } else {
      param <- paramNext
    }
  }
  return(list(par=paramNext, lastIt=i))
}

### Apply method

# Data (parallaxes in arcseconds; so distances are in parsecs)
wsd  <- 1.0e-3
data <- 1e-3*c(0.36, 1.11, 1.16, 1.31, 1.74, 1.94, 2.12, 2.72, 3.56, 4.30)
#data <- rnorm(n=10, mean=1/500, sd=wsd) # Above data were drawn like this

# Set up parameters
rlen  <- 1000 # prior length scale
Kboot <- 1000 # no. bootstrap samples

# Plot some example posteriors
pdf("bootstrap_posteriors.pdf", width=5, height=4)
par(mfrow=c(1,1), mar=c(3.5,4.2,1,1.25), oma=0.1*c(1,1,1,1), 
    mgp=c(2.3,0.9,0), cex=1.15)
z <- post.grid(data, wsd, rlen)
plot(z$r, z$post, type="l", lwd=3, col="grey70",
     xlab="r", ylab=expression(paste(P^symbol("*"), "(r | ", D[k], ")")), 
     ylim=c(0,1.05), xaxs="i", yaxs="i")
set.seed(555)
for(j in 1:7) {
  dataSample <- draw.sample(data)
  z <- post.grid(dataSample, wsd, rlen)  
  lines(z$r, z$post, type="l", lwd=1.2)
}
dev.off()

# Draw bootstrap samples and for each calculate the maximum of the
# log posterior using gradient ascent
rMode <- numeric(Kboot)
for(k in 1:Kboot) {
  dataSample <- draw.sample(data)
  rInit <- 1/median(dataSample) # initial distance for optimization
  opt <- grad.ascent(gradFunc=grad.log.d.post, param=rInit, d.post, 
                     dataSample, wsd, rlen) 
  cat(rInit, opt$par, opt$lastIt, "\n") # print initial, final, no. its
  rMode[k] <- opt$par
}
cat("sd of bootstrap samples =", sd(rMode), "\n")
pdf("bootstrap_histogram.pdf", 5,4)
par(mfrow=c(1,1), mar=c(3.5,3.5,1,1), oma=0.1*c(1,1,1,1), mgp=c(2.2,0.8,0), 
    cex=1.15)
hist(rMode, breaks=25, xlim=c(200, 1400), main="", 
     xlab=expression(hat(r)), ylab="frequency")
segments(x0=1/data, y0=0, x1=1/data, y1=25, lwd=3)
dev.off()

# Estimate distance to the star using the original data set
opt <- grad.ascent(gradFunc=grad.log.d.post, param=rInit, d.post, data, 
                   wsd, rlen) 
cat("Numerical estimate of mode =", opt$par, "\n")

# For comparison, compute sd of posterior from all data (need to normalize 
# it first). r must extend to where post is very small
r <- seq(from=0, to=2*rlen, by=0.1) 
post <- numeric(length(r))
for(i in 1:length(r)) {
  post[i] <- d.post(data=data, r=r[i], wsd, rlen)
}
rMean <- sum(post*r)/sum(post)
plot(r, post, type="l") # ensure we have sampled the full posterior
cat("sd of posterior =", sqrt(sum(post*(r-rMean)^2)/sum(post)), "\n")
