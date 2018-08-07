##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Plot the binomial distribution

# Plot P vs. r for fixed n for a range of p
n <- 10
r <- 0:n
pseq <- c(0.1, 0.2, 0.5, 0.8, 0.9)
pdf("dbinom1.pdf", 4, 4)
par(mfrow=c(1,1), mgp=c(2.0,0.8,0), mar=c(3.5,3.5,1,1), oma=0.1*c(1,1,1,1))
plot(r, r, type="n", xlim=c(0,max(r)), ylim=c(0,0.4), xlab="r", 
     ylab="P(r | p,n)")
for (p in pseq) {
  points(r, dbinom(x=r, size=n, prob=p), pch=20)
  lines(r,  dbinom(x=r, size=n, prob=p), lty=2)
}
text(c(1,2,5,8,9), c(0.39,0.30,0.285,0.30,0.39), pseq, pos=c(4,4,1,2,2))
dev.off()

# Plot P vs. p for fixed n for a range of r
p    <- seq(from=0, to=1, by=0.001)
rseq <- c(0,1,3,5)
n    <- 10
pdf("dbinom2.pdf", 4, 4)
par(mfrow=c(1,1), mgp=c(2.0,0.8,0), mar=c(3.5,3.5,1,1), oma=0.1*c(1,1,1,1))
plot(p, p, type="n", xlim=range(p), ylim=c(0,0.5), xlab="p", 
     ylab="P(r | p,n)")
for (r in rseq) {
  lines(p,  dbinom(x=r, size=n, prob=p))
}
text(c(0.08,0.15,0.25,0.45), c(0.45,0.35,0.29,0.27), rseq, pos=4)
dev.off()
