##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Demonstration of lm

x <- c(10.9, 12.4, 13.5, 14.6, 14.8, 15.6, 16.2, 17.5, 18.3, 18.6)
y <- c(24.8, 30.0, 31.0, 29.3, 35.9, 36.9, 42.5, 37.9, 38.9, 40.5)
model1 <- lm(y ~ x)
model1$fitted.values
predict(model1) # lm has a predict method
xv <- 10:20
#yv <- predict(model1, xv) # does not work: need to name the variable 'x'
yv <- predict(model1, data.frame(x=xv))
pdf("linear_regression.pdf", 8, 4)
par(mfrow=c(1,2), mgp=c(2.0,0.8,0), mar=c(3.5,3.5,1,1), oma=0.1*c(1,1,1,1))
plot(x, y)
lines(xv, yv)
# lm has an abline method so you can instead do abline(model1)
plot(x, model1$residuals, ylab="residuals", pch=16)
abline(h=0, col="grey60")
dev.off()
