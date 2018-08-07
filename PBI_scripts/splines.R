##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Demonstration of 1D splines

library(splines) # for bs

# Simulate data
set.seed(101)
x <- -10:10
y <- sin(2*pi*x/10) + 0.005*x^3 + 1*rnorm(length(x), 0, 1)
xp <- seq(from=min(x), to=max(x), length.out=1e3) # for plotting splines

# Discontinuous splines of orders 0,1,3 with definable knots (here 4)
# I use lm() to explicitly fit within each region.
pdf("splines_discontinuous.pdf", 4, 4)
par(mfrow=c(1,1), mar=c(3.5,3.5,0.5,1), oma=c(0.5,0.5,0.5,0.5), 
    mgp=c(2.2,0.8,0), cex=1.0)
plot(x,y)
knots=c(-3.33, 3.33)
abline(v=c(knots, range(x)), lty=3)
xSeg <- list(which(x<=knots[1]), which(x>knots[1] & x<knots[2]), 
             which(x>=knots[2]))
xpSeg <- list(which(xp<=knots[1]), which(xp>knots[1] & xp<knots[2]), 
              which(xp>=knots[2]))
for(s in 1:3) { # loop over the three subregions
  lines(xp[xpSeg[[s]]], predict(lm(y ~ 1, subset=xSeg[[s]]), 
                                newdata=data.frame(x=xp[xpSeg[[s]]])), 
        lwd=1.5, lty=2)
  lines(xp[xpSeg[[s]]], predict(lm(y ~ x, subset=xSeg[[s]]), 
                                newdata=data.frame(x=xp[xpSeg[[s]]])), 
        lwd=3, col="grey70")
  lines(xp[xpSeg[[s]]], predict(lm(y ~ x + I(x^2) + I(x^3), 
                                   subset=xSeg[[s]]), 
                                newdata=data.frame(x=xp[xpSeg[[s]]])), 
        lwd=1.5, col="black")
}
dev.off()

# Continuous splines of orders 1,3 with definable knots (here 4)
# bs() automatically uses extreme data points as additional knots.
pdf("splines_continuous.pdf", 4, 4)
par(mfrow=c(1,1), mar=c(3.5,3.5,0.5,1), oma=c(0.5,0.5,0.5,0.5), 
    mgp=c(2.2,0.8,0), cex=1.0)
plot(x,y)
knots=c(-3.33, 3.33)
abline(v=c(knots, range(x)), lty=3)
lines(xp, predict(lm(y ~ bs(x, knots=knots, degree=1)), 
                  newdata=data.frame(x=xp)), lwd=3, col="grey70")
lines(xp, predict(lm(y ~ bs(x, knots=knots, degree=3)), 
                  newdata=data.frame(x=xp)), lwd=1.5, col="black")
dev.off()            

# Exact splines (i.e knot at each point) of orders 1,3
pdf("splines_exact.pdf", 4, 4)
par(mfrow=c(1,1), mar=c(3.5,3.5,0.5,1), oma=c(0.5,0.5,0.5,0.5), 
    mgp=c(2.2,0.8,0), cex=1.0)
plot(x,y)
# Plots exact linear spline then exact cubic spline
lines(xp, approxfun(x, y, method="linear")(xp), lwd=3, col="grey70")
lines(xp, splinefun(x, y)(xp), lwd=1.5, col="black")
dev.off()

# Smoothing splines with various degrees of freedom
pdf("splines_smoothing.pdf", 4, 4)
par(mfrow=c(1,1), mar=c(3.5,3.5,0.5,1), oma=c(0.5,0.5,0.5,0.5), 
    mgp=c(2.2,0.8,0), cex=1.0)
plot(x,y)
lines(predict(smooth.spline(x, y, df= 2), xp), lwd=1.5, lty=2, col="black")
lines(predict(smooth.spline(x, y, df= 4), xp), lwd=2,   lty=3, col="black")
lines(predict(smooth.spline(x, y, df= 8), xp), lwd=1.5, lty=1, col="black")
lines(predict(smooth.spline(x, y, df=15), xp), lwd=2,   lty=1, col="grey70")
dev.off()
