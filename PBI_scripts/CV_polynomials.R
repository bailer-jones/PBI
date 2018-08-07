##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Demonstration of cross-validation model selection using polynomials

# Function to evaluate polynomial with parameter beta at x.
# poly() as used here excludes x^0 - the constant, beta[1] - so I remove
# it from the matrix multiplication, then add it explicitly afterwards.
polyeval <- function(x, beta) {
  return(poly(x=x, degree=length(beta)-1, raw=TRUE) %*% beta[-1] + beta[1])
}

pdf("CV_polynomials.pdf", 10, 12)
par(mfrow=c(4,3), mar=c(3.5,3.5,1.5,0.5), oma=c(0.5,0.5,1,0.5), 
    mgp=c(2.2,0.8,0), cex=1.0)
set.seed(63)

# Simulate data: 4th order polynomial
beta <- c(0, 5, 1, -2, 0.1)
Ndat  <- 25
sigma <- 250
x <-  runif(Ndat, min=-8, max=25)
y <-  polyeval(x, beta) + rnorm(Ndat, 0, sigma)
# xp, yp just for plotting
xp <- seq(from=-10, to=30, length.out=1e3)
yp <- polyeval(xp, beta)
plot(x,y)
#lines(xp, yp, col="red")

# Do CV. Plot all fits for each j in a separate panel
jmax <- 9 # evaluate all polynomials with order from 0 to jmax
# Do j=0 separately
rss0 <- 0
plot(x, y, main=c("J = 0"))
for(i in 1:Ndat) {
  pred <- mean(y[-i])
  rss0 <- rss0 + (pred-y[i])^2
  abline(h=pred)
}
# Do j=1:jmax
rss <- vector(mode="numeric", length=jmax)
for(j in 1:jmax) {
  rss[j] <- 0
  plot(x, y, main=paste("J =", j))
  for(i in 1:Ndat) {
    # poly() used in lm() does include x^0 term
    mod  <- lm(y ~ poly(x, j, raw=TRUE), subset=-i)
    pred <- predict(mod, newdata=data.frame(x=x[i]))    
    rss[j] <- rss[j] + (pred - y[i])^2
    lines(xp, predict(mod, newdata=data.frame(x=xp)))
  }
}
plot(0:jmax, log10(sqrt(c(rss0, rss)/Ndat)), xlab="J", ylab="log10 RMS")
dev.off()
