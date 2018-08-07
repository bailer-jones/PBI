##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Demonstration of ridge regression for 1D polynomials

# Function to evaluate polynomial with parameter beta at x.
# poly() as used here excludes x^0 - the constant, beta[1] - so I remove
# it from the matrix multiplication, then add it explicitly afterwards.
polyeval <- function(x, beta) {
  return(poly(x=x, degree=length(beta)-1, raw=TRUE) %*% beta[-1] + beta[1])
}

# Function to solve for ridge regression parameters
# given data vector y of length N, N*J data matrix X, 
# and lambda (scalar).
param.ridge <- function(y, X, lambda) {
  return( solve(t(X)%*%X + diag(lambda, nrow=ncol(X))) %*% t(X) %*% y )
}

pdf("ridge_regression.pdf", 12, 4)
par(mfrow=c(1,3), mar=c(3.5,3.5,0.5,0.5), oma=c(0.5,0.1,0.5,0.5), 
    mgp=c(2.2,0.8,0), cex=1.2)

set.seed(100)

# Simulate data: 3rd order (cubic) polynomial
beta <- c(0, -5, 1, -0.1)
Ndat  <- 15
sigma <- 100
x <-  runif(Ndat, min=-8, max=25)
y <-  polyeval(x, beta) + rnorm(Ndat, 0, sigma)
# xp, yp just for plotting
xp <- seq(from=-10, to=27, length.out=1e3)
yp <- polyeval(xp, beta)

# Plot data and true model
plot(x, y, ylim=c(-580, 240))
lines(xp, yp, lty=3, lwd=2)

# Build matrices, centre y, centre and scale each power of x. "degree" in 
# poly defines the order of the polynomial we use in both solutions. 
ys <- scale(y, scale=FALSE)
X  <- scale(poly(x, degree=6, raw=TRUE)) # Ndat * 6 matrix
XP <- poly(xp, degree=6, raw=TRUE)       # 1e3  * 6 matrix
XP <- t( (t(XP) - attr(X,"scaled:center")) / attr(X,"scaled:scale") ) 
# Solve for OLS and ridge regression parameters and calculate residuals
betaOLS  <- param.ridge(ys, X, lambda=0)
residOLS <- X %*% betaOLS - ys
betaRidge  <- param.ridge(ys, X, lambda=1)
residRidge <- X %*% betaRidge - ys

# Plot raw (non-standardized) data together with fitted curves.
# Latter is done by plotting de-centered model predictions at raw xp.
plot(x, y, ylim=c(-580, 240))
lines(xp, XP %*% betaOLS   + attr(ys, "scaled:center"), lwd=2)
lines(xp, XP %*% betaRidge + attr(ys, "scaled:center"), lwd=2.5, lty=2)

# Plot residuals
plot(x, y, type="n", ylim=range(c(residOLS, residRidge)), 
     ylab="y residuals")
abline(h=0, col="grey")
points(x, residOLS,   pch=20)
points(x, residRidge, pch=4)

dev.off()

# Print coefficients, sum of squares of coefficients, and RMS of residuals
format(data.frame(betaOLS, betaRidge), digits=3)
cat("beta^2 = ", sum(betaOLS^2), sum(betaRidge^2), "\n")
cat("RMS =    ", sqrt(mean(residOLS^2)), sqrt(mean(residRidge^2)), "\n")
