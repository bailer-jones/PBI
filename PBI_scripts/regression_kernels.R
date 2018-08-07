##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

#### Demonstration of regression kernels

library(KernSmooth) # for locpoly

pdf("regression_kernels.pdf", 4, 4)
par(mfrow=c(1,1), mar=c(3.5,3.5,0.5,1), oma=c(0.5,0.5,0.5,0.5), 
    mgp=c(2.2,0.8,0), cex=1.0)

# Simulate data
x <- seq(from=0, to=1, length.out=100)
f <- function(x){0.5 + 0.4*sin(2*pi*x)} 
set.seed(10)
y <- f(x) + rnorm(n=x, mean=0, sd=0.2)
plot(x, y)

# Constant regression using a constant kernel of specified bandwidth.
xp <- seq(from=min(x), to=max(x), length.out=1e3)
lines(xp, ksmooth(x=x, y=y, kernel="box", bandwidth=0.1, x.points=xp)$y, 
      lwd=1.5)

# Local polynomial (oder=degree) regression using a Gaussian kernel
# with specified bandwidth and polynomial degree (here constant)
mod <- locpoly(x=x, y=y, degree=1, bandwidth=0.1, gridsize=1e3)
lines(mod$x, mod$y, lwd=2, lty=2)

# Local polynomial (oder=degree) regression using fraction span of data
# around each prediction point, i.e. span*length(x) nearest neighbours.
# It won't work if span is too small.
xp <- seq(from=min(x), to=max(x), length.out=1e3)
yp <- predict(loess(y ~ x, span=0.2, degree=1), newdata=data.frame(x=xp))
lines(xp, yp, lwd=2, col="grey60")

dev.off()
