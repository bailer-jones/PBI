##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Illustration with the binomial distribution that p-values depend on
##### sample size

pdf("pvalues_binomial.pdf", 4, 4)
par(mfrow=c(1,1), mgp=c(2.0,0.8,0), mar=c(3.5,3.5,1,1), oma=0.1*c(1,1,1,1))
sampSize <- seq(from=10, to=200, by=10)
fracLim  <- c(0.4, 0.3, 0.2)
plot(range(sampSize), c(-10,0), type="n", xlab="sample size, n", 
     ylab="log(p-value)")
for(f in fracLim) {
  pVal <- pbinom(q=f*sampSize, size=sampSize, prob=0.5)
  points(sampSize, log10(pVal), pch=20)
  lines(sampSize,  log10(pVal), lty=2)
}
text(c(170, 170, 90), c(-1.8, -6.7, -8.2), fracLim, pos=4)
dev.off()
