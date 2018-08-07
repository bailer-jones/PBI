##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Infer posterior PDF over amplitude and background parameters

# Define function to return true signal at position x (generative model)
signal <- function(x, a, b, x0, w, t) {
  t*(a*exp(-(x-x0)^2/(2*w^2)) + b)
}

# Define function to return (natural) log posterior over (a,b).
# Prior on a and b: P(a,b) = const if a>0 and b>0, = 0 otherwise.
# Likelihood for one point is Poisson with mean d(x), so total 
# likelihood is their product. Unnormalized posterior is product of these.
# d and x are equal length vectors (or scalars). The rest are scalars.
logupost <- function(d, x, a, b, x0, w, t) {
  if(a<0 || b <0) {return(-Inf)} # the effect of the prior
  sum(dpois(d, lambda=signal(x, a, b, x0, w, t), log=TRUE))
}

# Set model parameters (true and fixed)
x0    <- 0 # centre of peak
w     <- 1 # sd of peak
atrue <- 2 # amplitude
btrue <- 1 # background
t     <- 5 # scale factor (exposure time -> sets SNR)

# Simulate some data (by drawing from the likelihood)
set.seed(205)
xdat  <- seq(from=-7*w, to=7*w, by=0.5*w)
strue <- signal(xdat, atrue, btrue, x0, w, t)
ddat  <- rpois(length(strue), strue)

# Define sampling grid to compute posterior (will be normalized
# over this range too). uniGrid spans the range 0-1 with Nsamp 
# points. This is then scaled to cover the ranges alim and blim.
alim  <- c(0.0, 4.0)
blim  <- c(0.5, 1.5)
Nsamp <- 1e2
uniGrid <- seq(from=1/(2*Nsamp), to=1-1/(2*Nsamp), by=1/Nsamp)
delta_a <- diff(alim)/Nsamp 
delta_b <- diff(blim)/Nsamp
a <- alim[1] + diff(alim)*uniGrid 
b <- blim[1] + diff(blim)*uniGrid 

# Compute log unnormalized posterior, z = ln P^*(a,b|D), on a regular grid
z <- matrix(data=NA, nrow=length(a), ncol=length(b))
for(j in 1:length(a)) {
  for(k in 1:length(b)) {
    z[j,k] <- logupost(ddat, xdat, a[j], b[k], x0, w, t)
  }
}
z <- z - max(z) # set maximum to zero

# Compute normalized marginalized posteriors, P(a|D) and P(b|D)
# by summing over other parameter. Normalize by gridding.
p_a_D <- apply(exp(z), 1, sum)
p_a_D <- p_a_D/(delta_a*sum(p_a_D))
p_b_D <- apply(exp(z), 2, sum)
p_b_D <- p_b_D/(delta_b*sum(p_b_D))

# Compute mean, standard deviation, covariance, correlation, of a and b
mean_a <- delta_a * sum(a * p_a_D)
mean_b <- delta_b * sum(b * p_b_D)
sd_a   <- sqrt( delta_a * sum((a-mean_a)^2 * p_a_D) )
sd_b   <- sqrt( delta_b * sum((b-mean_b)^2 * p_b_D) )
# To calculate the covariance I need to normalize P(a,b|D) = exp(z).
# I do it here by brute force with two loops (there are better ways in R).
# The normalization constant is Z = delta_a*delta_b*sum(exp(z)).
# This is independent of (a,b) so can be calculated outside of the loops.
# The factor delta_a*delta_b will just cancel in the expression for 
# cov_ab, so I omit it entirely.
cov_ab <- 0
for(j in 1:length(a)) {
  for(k in 1:length(b)) {
    cov_ab <- cov_ab + (a[j]-mean_a)*(b[k]-mean_b)*exp(z[j,k])
  }
}
cov_ab <- cov_ab / sum(exp(z))
rho_ab <- cov_ab / (sd_a * sd_b)
cat("  a = ", mean_a, "+/-", sd_a, "\n")
cat("  b = ", mean_b, "+/-", sd_b, "\n")
cat("rho = ", rho_ab, "\n")

# Compute normalized conditional posteriors, P(a|b,D) and P(b|a,D)
# using true values of conditioned parameters. Vectorize(func, par)
# makes a vectorized function out of func in the parameter par.
p_a_bD <- exp(Vectorize(logupost, "a")(ddat, xdat, a, btrue, x0, w, t))
p_a_bD <- p_a_bD/(delta_a*sum(p_a_bD))
p_b_aD <- exp(Vectorize(logupost, "b")(ddat, xdat, atrue, b, x0, w, t))
p_b_aD <- p_b_aD/(delta_b*sum(p_b_aD))

# Make plots

pdf("signal_background_estimation.pdf", 7, 7)
# Plot true model and data
par(mfrow=c(2,2), mgp=c(2,0.8,0), mar=c(3.5,3.5,1,1), oma=0.1*c(1,1,1,1))
xplot <- seq(from=min(xdat), to=max(xdat), by=0.05*w)
splot <- signal(xplot, atrue, btrue, x0, w, t)
plot(xplot, splot, ylim=range(c(splot, ddat)), xlab="x", ylab="s or d", 
     type="l", col="grey", lwd=2)
points(xdat, ddat)
# Plot unnormalized 2D posterior as contours.
# Note that they are labelled by posterior density relative to peak, 
# NOT by how much probabilty they enclose.
contour(a, b, exp(z), nlevels=5, labcex=0.5, lwd=2, xlab="amplitude, a", 
        ylab="background, b")
abline(v=2,h=1,col="grey")
# Plot the 1D marginalized posteriors
plot(b, p_b_D, xlab="background, b", yaxs="i", 
     ylim=1.05*c(0,max(p_b_D, p_b_aD)), ylab="P(b | D)  and  P(b | a,D)", 
     type="l", lwd=2)
lines(b, p_b_aD, lwd=2, lty=2)
abline(v=btrue, col="grey")
plot(a, p_a_D, xlab="amplitude, a", yaxs="i", 
     ylim=1.05*c(0,max(p_a_D, p_a_bD)), ylab="P(a | D)  and  P(a | b,D)", 
     type="l", lwd=2)
lines(a, p_a_bD, lwd=2, lty=2)
abline(v=atrue, col="grey")

dev.off()
