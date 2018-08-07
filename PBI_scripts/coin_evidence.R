##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Calculate the evidence for the coin problem

phi <- 1/2
n   <- c(5,10,25,50,100)
pdf("coin_evidence.pdf", 4, 4)
par(mfrow=c(1,1), mgp=c(2.0,0.8,0), mar=c(3.5,3.5,1,1), oma=0.1*c(1,1,1,1))
plot(c(0,1), c(0,1), type="n", xlim=c(0,1), ylim=c(-2,1), xlab="f = r/n", 
     ylab=expression(paste(log, B[12][u])))
for(i in 1:length(n)) {
  r  <- 0:n[i]
  f  <- r/n[i]
  BF <- (n[i]+1)*dbinom(x=r, size=n[i], prob=phi)
  points(f, log10(BF), cex=0.5, pch=13+i)
  lines(f, log10(BF), lty=6-i)
}
dev.off()
