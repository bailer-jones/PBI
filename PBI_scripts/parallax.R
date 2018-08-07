##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Plot parallax posterior PDFs

# Define functions
d.like <- function(w, r, wsd) dnorm(x=w, mean=1/r, sd=wsd)
d.prior2 <- function(r, rmax) ifelse(r>0 & r<rmax, 1/rmax, 0) # normalized. Vectorized in r and rmax

# Plot likelihood as a function of r together with a Gaussian with
# mean 1/r and standard deviation given by the first order Taylor expansion
pdf("parallax_likelihood_1.pdf", width=6, height=5)
par(mfrow=c(1,1), mar=c(3.5,4.2,1,1.25), oma=c(0.5,0.5,0.5,0.5), mgp=c(2.3,0.9,0), cex.lab=1.4, cex.axis=1.2)
wTrue <- 0.1
w <- seq(from=-0.1, to=0.3, length.out=1e4)
plot(w, w, type="n", xlab="1/r", ylab=expression(paste(P, "(", omega1, " | r)")), ylim=c(0,21), xaxs="i", yaxs="i")
lines(w, d.like(w, r=1/wTrue, wsd=0.02), type="l", lwd=2)
segments(wTrue-2*0.02, 0, wTrue-2*0.02, d.like(wTrue-2*0.02, r=1/wTrue, wsd=0.02), lwd=2)
segments(wTrue+2*0.02, 0, wTrue+2*0.02, d.like(wTrue+2*0.02, r=1/wTrue, wsd=0.02), lwd=2)
lines(w, d.like(w, r=1/wTrue, wsd=0.05), type="l", lwd=2.5, col="grey60")
segments(wTrue-2*0.05, 0, wTrue-2*0.05, d.like(wTrue-2*0.05, r=1/wTrue, wsd=0.05), lwd=2.5, col="grey60")
segments(wTrue+2*0.05, 0, wTrue+2*0.05, d.like(wTrue+2*0.05, r=1/wTrue, wsd=0.05), lwd=2.5, col="grey60")
dev.off()
pdf("parallax_likelihood_2.pdf", width=6, height=5)
par(mfrow=c(1,1), mar=c(3.5,4.2,1,1.25), oma=c(0.5,0.5,0.5,0.5), mgp=c(2.3,0.9,0), cex.lab=1.4, cex.axis=1.2)
rp <- seq(from=-5, to=30, length.out=1e4)
plot(rp, rp, type="n", xlab="r", ylab=expression(paste(P, "(", omega1, " | r)")), ylim=c(0,21), xaxs="i", yaxs="i")
lines(rp, d.like(w=1/rp, r=1/wTrue, wsd=0.02), type="l", lwd=2)
segments(1/(wTrue-2*0.02), 0, 1/(wTrue-2*0.02), d.like(wTrue-2*0.02, r=1/wTrue, wsd=0.02), lwd=2)
segments(1/(wTrue+2*0.02), 0, 1/(wTrue+2*0.02), d.like(wTrue+2*0.02, r=1/wTrue, wsd=0.02), lwd=2)
lines(rp, d.like(w=1/rp, r=1/wTrue, wsd=0.05), type="l", lwd=2.5, col="grey60")
#segments(1/(wTrue-2*0.05), 0, 1/(wTrue-2*0.05), d.like(wTrue-2*0.05, r=1/wTrue, wsd=0.05), col="red")
segments(1/(wTrue+2*0.05), 0, 1/(wTrue+2*0.05), d.like(wTrue+2*0.05, r=1/wTrue, wsd=0.05), lwd=2, col="grey")
# factor of 100 (ratio of standard deviations) required to scale the following to the likelihoods
lines(rp, 1e2*dnorm(x=rp, mean=1/wTrue, sd=0.02/wTrue^2), lwd=2, lty=2)
lines(rp, 1e2*dnorm(x=rp, mean=1/wTrue, sd=0.05/wTrue^2), lwd=2.5, lty=2, col="grey60")
dev.off()

# Plot unnormalised improper posterior for various f_true
pdf("ud.post_runifPrior.pdf", width=6, height=5)
par(mfrow=c(1,1), mar=c(3.5,4.2,1,1.25), oma=c(0.5,0.5,0.5,0.5), mgp=c(2.3,0.9,0), cex.lab=1.4, cex.axis=1.2)
w <- 0.01
fVals <- c(0.1, 0.2, 0.5, 1.0)
rmax <- 310
r <- seq(from=0, to=rmax, length.out=1e4) 
plot(r, r, type="n", xlab="r", ylab=expression(paste(P^symbol("*"), "(r | ", omega1, ")")), ylim=c(0,1.05), xaxs="i", yaxs="i")
for(f in fVals) {
  wsd <- w*f
  udPost <- d.like(w, r, wsd)
  lines(r, udPost/max(udPost), type="l", lwd=1.5)
}
text(x=c(140, 155, 195, 230), y=c(0.15, 0.4, 0.7, 0.9), labels=fVals)
dev.off()

# Plot unnormalised proper (i.e. truncated) posterior for various f_true 
pdf("ud.post_runiftruncPrior.pdf", width=6, height=5)
par(mfrow=c(1,1), mar=c(3.5,4.2,1,1.25), oma=c(0.5,0.5,0.5,0.5), mgp=c(2.3,0.9,0), cex.lab=1.4, cex.axis=1.2)
w <- 0.01
fVals <- c(0.1, 0.2, 0.5, 1.0)
rmax <- 1e3
r <- seq(from=0, to=1.1*rmax, length.out=1e4) 
plot(r, r, type="n", xlab="r", ylab=expression(paste(P^symbol("*"), "(r | ", omega1, ")")), ylim=c(0,1.05), xaxs="i", yaxs="i")
udPost <- d.like(w=-0.01, r, wsd=0.0025)*d.prior2(r, rmax) # not normalized
lines(r, udPost/max(udPost), type="l", lty=2, lwd=2)
for(f in fVals) {
  wsd <- w*f
  udPost <- d.like(w, r, wsd)*d.prior2(r, rmax) # not normalized
  lines(r, udPost/max(udPost), type="l", lwd=1.5)
}
text(x=c(140, 175, 210, 230), y=c(0.1, 0.4, 0.7, 0.9), labels=fVals)
dev.off()

# Plot normalised proper (i.e. truncated) posterior for various f_true 
pdf("d.post_runiftruncPrior.pdf", width=6, height=5)
par(mfrow=c(1,1), mar=c(3.5,4.2,1,1.25), oma=c(0.5,0.5,0.5,0.5), mgp=c(2.3,0.9,0), cex.lab=1.4, cex.axis=1.2)
w <- 0.01
fVals <- c(0.1, 0.2, 0.5, 1.0)
rmax <- 1e3
r <- seq(from=0, to=1.1*rmax, length.out=1e4) 
plot(r, r, type="n", xlab="r", ylab=expression(paste(P, "(r | ", omega1, ")")), ylim=c(0,4e-2), xaxs="i", yaxs="i")
udPost <- d.like(w=-0.01, r, wsd=0.0025)*d.prior2(r, rmax) # not normalized
lines(r, udPost/(r[length(r)]*mean(udPost)), type="l", lty=2, lwd=2) # see below for normalization
for(f in fVals) {
  wsd <- w*f
  udPost <- d.like(w, r, wsd)*d.prior2(r, rmax) # not normalized
  # Normalization constant = sum(p*dr), where dr is bin size
  dPost  <- udPost/(r[length(r)]*mean(udPost)) # correct only if r sampling is uniform, starts at 0, and covers full range of prior
  lines(r, dPost, type="l", lwd=1.5)
}
dev.off()
