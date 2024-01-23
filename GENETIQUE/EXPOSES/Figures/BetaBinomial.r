# Bayesian inference in the BetaBinomial model

# Parms & data
n = 30
x = 4
#a = .5; b = .5
a = 1; b = 1
a = 10; b = 5
Lwd1 = 3
Lwd2 = 4

# Exact posterior and prior
PlotName = paste('BetaBinom-n', n, '-x', x, '-a', a, '-b', b, sep='')
pp = seq(0, 1, by=1e-3)
Prior = dbeta(pp, a, b)
Posterior = dbeta(pp, a+x, b+n-x)

# Plots
postscript(paste(PlotName, "-Prior.ps", sep=''))
plot(pp, Posterior, type='l', col=0, lwd=Lwd2, lty=1, main='', xlab='', ylab='')
lines(pp, Prior, col=4, lwd=Lwd2)
dev.off()

postscript(paste(PlotName, "-Data.ps", sep=''))
plot(pp, Posterior, type='l', col=0, lwd=Lwd2, lty=1, main='', xlab='', ylab='')
lines(pp, Prior, col=4, lwd=Lwd2)
abline(v=x/n, col=1, lwd=Lwd2, lty=2)
dev.off()

postscript(paste(PlotName, "-Posterior.ps", sep=''))
plot(pp, Posterior, type='l', col=0, lwd=Lwd2, lty=1, main='', xlab='', ylab='')
lines(pp, Prior, col=4, lwd=Lwd2)
lines(pp, Posterior, col=2, lwd=Lwd2)
abline(v=x/n, col=1, lwd=Lwd2, lty=2)
dev.off()

postscript(paste(PlotName, "-Credibility.ps", sep=''))
plot(pp, Posterior, type='l', col=0, lwd=Lwd2, lty=1, main='', xlab='', ylab='')
lines(pp, Prior, col=4, lwd=Lwd2)
lines(pp, Posterior, col=2, lwd=Lwd2)
abline(v=x/n, col=1, lwd=Lwd2, lty=2)
abline(v=(a+x)/(a+b+n), col=2, lwd=Lwd2, lty=2)
abline(v=qbeta(0.025, a+x, b+n-x), col=2, lwd=Lwd1, lty=2)
abline(v=qbeta(0.975, a+x, b+n-x), col=2, lwd=Lwd1, lty=2)
dev.off()

# Sampled posterior from prior
B = 1e4
P = rbeta(B, a, b)
PX = dbinom(rep(x, B), n, P)
wP = PX/sum(PX)
pB = seq(0, 1, by=1/sqrt(B))
SampPost0 = rep(0, 1+length(pB))
for (i in (1:length(pB)))
  {
  SampPost0[i] = sum(wP[which(P < pB[i])])
  }
# Sampled posterior from true posterior
B = 1e4
P = rbeta(B, a+x, b+n-x)
PX = dbinom(rep(x, B), n, P)*dbeta(P, a, b)/dbeta(P, a+x, b+n-x)
wP = PX/sum(PX)
pB = seq(0, 1, by=1/sqrt(B))
SampPost1 = rep(0, 1+length(pB))
for (i in (1:length(pB)))
  {
  SampPost1[i] = sum(wP[which(P < pB[i])])
  }

# Plots
Max = max(max(sqrt(B)*diff(SampPost0)), max(sqrt(B)*diff(SampPost1)))
postscript(paste(PlotName, "-SampPost-0.ps", sep=''))
plot(pB, sqrt(B)*diff(SampPost0), type='s', ylim=c(0, 1.1*Max), lwd=Lwd1, col=2, main='', xlab='', ylab='')
lines(pp, Prior, col=4, lwd=Lwd2)
lines(pp, Posterior, col=2, lwd=Lwd1, lty=2)
abline(v=x/n, col=1, lwd=Lwd2, lty=2)
#abline(v=P%*%wP, col=2, lwd=Lwd2, lty=2)
dev.off()

postscript(paste(PlotName, "-SampPost-1.ps", sep=''))
plot(pB, sqrt(B)*diff(SampPost1), type='s', ylim=c(0, 1.1*Max), lwd=Lwd1, col=2, main='', xlab='', ylab='')
lines(pp, Prior, col=4, lwd=Lwd2)
lines(pp, Posterior, col=2, lwd=Lwd1, lty=2)
abline(v=x/n, col=1, lwd=Lwd2, lty=2)
#abline(v=P%*%wP, col=2, lwd=Lwd2, lty=2)
dev.off()

