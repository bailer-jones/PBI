##### C.A.L. Bailer-Jones, Practical Bayesian Inference, 2017
##### Version 2017-03-24. CC-BY-4.0 licence (see README file)

##### Apply Bayes theorem to infer cancer probability given test result

# Return P(M|D) given P(D|M), P(D|M'), P(M)
post <- function(p_d_m, p_d_nm, p_m) {
  p_nm  <- 1-p_m
  oddsr <- (p_d_m * p_m) / (p_d_nm * p_nm)
  p_m_d <- 1/(1+1/oddsr)
  return(p_m_d)
}

# Vary reliability of test P(D|M)
p_d_m  <- seq(0.0, 1.0, 0.01) # prob. true positive
p_d_nm <- 0.07                # prob. false positive
p_m    <- 8/1000              # prior probability of m
p_m_d  <- post(p_d_m, p_d_nm, p_m)
pdf("cancer_test_1.pdf", 4, 4)
par(mfrow=c(1,1), mgp=c(2.0,0.8,0), mar=c(3.5,3.5,1,1), oma=0.1*c(1,1,1,1))
plot(p_d_m, p_m_d, type="l", lwd=2, xaxs="i", yaxs="i",
     xlab="P(positive result | Cancer)", ylab="P(Cancer | positive result)")
dev.off()

# Vary false positive probability P(D|!M)
p_d_m  <- 0.9               # prob. true positive
p_d_nm <- 10^seq(-4,0,0.02) # prob. false positive
p_m    <- 8/1000            # prior probability of m
p_m_d  <- post(p_d_m, p_d_nm, p_m)
pdf("cancer_test_2.pdf", 4, 4)
par(mfrow=c(1,1), mgp=c(2.0,0.8,0), mar=c(3.5,3.5,1,1), oma=0.1*c(1,1,1,1))
plot(log10(p_d_nm), p_m_d, type="l", lwd=2, ylim=c(0,1), yaxs="i",
     xlab="log[P(positive result | no Cancer)]", 
     ylab="P(Cancer | positive result)")
lines(log10(p_d_nm), post(p_d_m=1.0, p_d_nm, p_m), lty=2)
dev.off()
