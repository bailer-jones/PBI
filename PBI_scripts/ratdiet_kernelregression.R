##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Application of kernel regression to the rat.diet data

library(fields) # for rat.diet
attach(rat.diet)
Ndat <- length(t)

# Calculate RMS using LOO-CV
# loess does not work with just one nearest neighbour, so I just do 2:Ndat.
# It also does not permit extrapolation, so I don't include the two test 
# sets which would each be just the two extreme points. One could set 
# surface="direct" as a loess option, but this gives worse fits for low k. 
rss <- rep(NA, Ndat)
kRange <- 2:Ndat
for(k in kRange) {
  rss[k] <- 0
  for(i in 2:(Ndat-1)){ # RSS will be sum of squares of Ndat-2 residuals
    pred   <- predict(loess(trt ~ t, span=k/Ndat, degree=1, subset=-i), 
                      newdata=data.frame(t=t[i]))
    rss[k] <- rss[k] + (pred - trt[i])^2
  }
}
rms <- sqrt((rss[kRange])/(Ndat-2))

pdf("ratdiet_kernelregression.pdf", 8, 4)
par(mfrow=c(1,2), mar=c(3.5,3.5,0.5,0.5), oma=c(0.5,0.5,0.5,0.5), 
    mgp=c(2.2,0.8,0), cex=1.0)

# Plot RMS vs number of neighbours
plot(kRange, rms, xlab="no. neighbours, k", ylab="RMS", ylim=c(0,4))
bd <- kRange[which.min(rms)]
abline(v=bd, col="grey")
text(bd, 5, bd, pos=4)
cat(bd, rms[which.min(rms)])

# Plot data with three different kernel smoothers
plot(t, trt)
xp <- seq(from=min(t), to=max(t), length.out=1e3)
yp <- predict(loess(trt ~ t, span=16/Ndat, degree=1), 
              newdata=data.frame(t=xp))
lines(xp, yp, lwd=2,   lty=2, col="black")
yp <- predict(loess(trt ~ t, span=32/Ndat, degree=1), 
              newdata=data.frame(t=xp))
lines(xp, yp, lwd=2,   lty=1, col="grey70")
yp <- predict(loess(trt ~ t, span=8/Ndat, degree=1), 
              newdata=data.frame(t=xp))
lines(xp, yp, lwd=1.5, lty=1, col="black")

dev.off()
detach(rat.diet)
