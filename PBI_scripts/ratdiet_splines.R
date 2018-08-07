##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Application of smoothing splines to the rat.diet data

library(fields) # for sreg and rat.diet
Ndat <- nrow(rat.diet)

# Calculate RMS using LOO-CV
dofValues <- seq(2,35,0.2)
Ndof <- length(dofValues)
rss  <- rep(x=0, times=Ndof)
for (k in 1:Ndof) {
  for(i in 1:Ndat) {  
    mod  <- sreg(rat.diet$t[-i], rat.diet$trt[-i], df=dofValues[k])
    pred <- predict(mod, rat.diet$t[i])
    rss[k] <- rss[k] + (pred - rat.diet$trt[i])^2
  }
}
rms <- sqrt(rss/Ndat)

pdf("ratdiet_splines.pdf", 8, 4)
par(mfrow=c(1,2), mar=c(3.5,3.5,0.5,0.5), oma=c(0.5,0.5,0.5,0.5), 
    mgp=c(2.2,0.8,0), cex=1.0)

# Plot RMS vs dof
plot(dofValues, rms, xlab="effective degrees of freedom, dof", ylab="RMS", 
     ylim=c(0,4))
bd <- dofValues[which.min(rms)]
abline(v=bd, col="grey")
text(bd, 5, bd, pos=4)
cat(bd, rms[which.min(rms)])

# Plot data with three different spline smoothers
plot(rat.diet$t, rat.diet$trt, xlab="t", ylab="trt")
xp <- seq(from=min(rat.diet$t), to=max(rat.diet$t), length.out=1e3)
lines(xp, predict(sreg(rat.diet$t, rat.diet$trt, df=bd),   xp), lwd=2,   
      lty=2, col="black")
lines(xp, predict(sreg(rat.diet$t, rat.diet$trt, df=bd/2), xp), lwd=2,   
      lty=1, col="grey70")
lines(xp, predict(sreg(rat.diet$t, rat.diet$trt, df=2*bd), xp), lwd=1.5, 
      lty=1, col="black")

dev.off()
