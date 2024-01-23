# Figures for Bonafede & al, 2015

rm(list=ls())

# Parms
K = 4
alpha = 2^(1:K)
lambda = 30
u.list = seq(0, 4, by=.01)
y.list = (0:100)
Lwd = 3

# par(mfrow=c(2, 1), mex=.5)
pdf('BPR15-Biometrics-MixtGamma.pdf')
plot(u.list, dgamma(u.list, alpha[1], alpha[K]), col=0, 
     xlab='', ylab='', ylim=c(0, 2))
for(k in 1:K){
  lines(u.list, dgamma(u.list, alpha[k], alpha[k]), col=k, lwd=Lwd)
}
abline(v=1, col=8, lty=3, lwd=Lwd)
dev.off()

pdf('BPR15-Biometrics-MixtNegBinom.pdf')
plot(y.list, dnbinom(y.list, mu=lambda, size=alpha[K]), col=0, xlab='', ylab='')
for(k in 1:K){
  lines(y.list, dnbinom(y.list, mu=lambda, size=alpha[k]), col=k, lwd=Lwd, type='s')
}
abline(v=lambda, col=8, lty=3, lwd=Lwd)
dev.off()
# B = 1e4
# alpha = 5; U = rgamma(B, alpha, alpha); mean(U); var(U)
# lambda = 10
# Y1 = rpois(B, lambda*U); Y2 = rnbinom(B, mu=lambda, size=alpha)
# mean(Y1); mean(Y2); var(Y1); var(Y2)
