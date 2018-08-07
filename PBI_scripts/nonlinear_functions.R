##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Demonstration of linear regression with nonlinear functions of the
##### data using lm()

x <- c(0.9, 2.1, 2.9, 3.8, 5.3, 6.0, 7.0, 8.4, 9.1, 9.8, 11.1, 11.9, 13.2, 
       13.8, 14.8)
y <- c(19.5, 13.8, 16.6, 11.6, 8.3, 9.9, 7.6, 6.6, 6.2, 5.4, 5.9, 5.5, 5.6, 
       4.4, 4.7)
par(mfrow=c(2,2), mgp=c(2.0,0.8,0), mar=c(3.5,3.5,1,1), oma=0.1*c(1,1,1,1))
plot(x,y) # looks nonlinear...

# Straight line model
model1 <- lm(y ~ x)
summary(model1) # r^2 is a bit low
abline(model1, col="blue", lw=2)
plot(x, model1$residuals, pch=20, col="blue", ylim=c(-5,5)) 
# we see structure in the residuals
abline(h=0)

# Quadratic model
model2 <- lm(y ~ x + I(x^2)) # note meaning of "+" and I()
summary(model2) # looks better than linear
xv  <- seq(from=0, to=30, by=0.1) # generate data for prediction, 
                                  # as abline doesn't do curves
yv2 <- predict(model2, data.frame(x=xv))
plot(x,y)
lines(xv, yv2, col="red", lw=2)  
plot(x, model2$residuals, pch=20, col="red", ylim=c(-5,5)) 
# we see less structure in residuals
abline(h=0)

# Cubic model
model3 <- lm(y ~ x + I(x^2) + I(x^3))
summary(model3) 
# no evidence for cubic term; r^2 hardly drops compared to quadratic
yv3 <- predict(model3, data.frame(x=xv))
plot(x,y)
lines(xv, yv3, col="magenta", lw=2)  
plot(x, model3$residuals, pch=20, col="magenta", ylim=c(-5,5)) 
# this looks no better than quadratic
abline(h=0)

# Quadratic model without the linear term
model2b <- lm(y ~ I(x^2)) # can we drop linear term in the quadratic model?
summary(model2b)
yv2b <- predict(model2b, data.frame(x=xv))
plot(x,y)
lines(xv, yv2b, col="brown", lw=2) # no we can't!
plot(x, model2b$residuals, pch=20, col="brown", ylim=c(-5,5)) 
abline(h=0)

# Plot linear and quadratic results together
pdf("nonlinear_functions.pdf", 8, 4)
par(mfrow=c(1,2), mgp=c(2.0,0.8,0), mar=c(3.5,3.5,1,1), oma=0.1*c(1,1,1,1))
plot(x,y)
abline(model1, lw=2)
lines(xv, yv2, lw=2, lty=2)
plot(x,   model1$residuals, pch=20, ylim=c(-4.5,4.5), ylab="residuals")
points(x, model2$residuals, pch=4)
abline(h=0, col="grey60")
dev.off()
