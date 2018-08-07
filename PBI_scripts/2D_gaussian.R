##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Plot bivariate Gaussian as a 3D mesh plot
##### and as contours of constant probability density

library(mvtnorm) # for dmvnorm
sigma.x <- 1
sigma.y <- 1
rho <- 0.5 # correlation coefficient
Cov <- matrix(data=c(sigma.x^2, rho*sigma.x*sigma.y, rho*sigma.x*sigma.y,
                     sigma.y^2), nrow=2, ncol=2)
Nsig  <- 3.5
Nsamp <- 100
x <- seq(from=-Nsig*sigma.x, to=Nsig*sigma.x, length.out=Nsamp)
y <- seq(from=-Nsig*sigma.y, to=Nsig*sigma.y, length.out=Nsamp)
z <- matrix(dmvnorm(x=expand.grid(x,y), mean=c(0,0), sigma=Cov), 
            nrow=length(x), ncol=length(y))
z <- z/max(z)

pdf("2D_gaussian_3Dmesh.pdf", 4, 4)
par(mfrow=c(1,1), mar=c(1,1,1,1), oma=c(0,0,0,0), mgp=c(2.2,0.8,0), cex=1.0) 
persp(x=x, y=y, z=z, phi=20, theta=20, d=5, zlab="density")
dev.off()

pdf("2D_gaussian_contours.pdf", 4, 4)
par(mfrow=c(1,1), mgp=c(2.0,0.8,0), mar=c(3,3,1,1), oma=0.1*c(1,1,1,1))
contour(x, y, z, asp=1, xlim=c(-3,3), ylim=c(-3,3), xlab="x", ylab="y")
dev.off()
