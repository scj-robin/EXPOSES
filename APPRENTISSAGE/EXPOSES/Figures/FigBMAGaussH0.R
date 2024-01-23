source('/media/donnees/RECHERCHE/APPRENTISSAGE/PGM/F_FitMixtGaussNull.R')


par(mfrow=c(3, 4))
n = 1000
X = c(); xx = seq(-5, 15, by=0.1); f = rep(0, length(xx)); f1 = f

m = 0; s = 1; p0 = 0.5
fg = p0*dnorm(xx); plot(xx, fg, col=2, lwd=2, type='l')
f = f+fg; X = c(X, rnorm(round(p0*n), mean=m, sd=s))

umin = 1; umax = 7; p = 0.1
fg = p*dunif(xx, min=umin, max=umax); lines(xx, fg, col=3, lwd=2)
f = f+fg; f1 = f1+fg; X = c(X, runif(round(p*n), min=umin, max=umax))

m = 5; s = 1.5; p = 0.1
fg = p*dnorm(xx, mean=m, sd=s); lines(xx, fg, col=4, lwd=2)
f = f+fg; f1 = f1+fg; X = c(X, rnorm(round(p*n), mean=m, sd=s))

m = 2; s = .5; p = 0.1
fg = p*dnorm(xx, mean=m, sd=s); lines(xx, fg, col=5, lwd=2)
f = f+fg; f1 = f1+fg; X = c(X, rnorm(round(p*n), mean=m, sd=s))

a = 3; b = 1; p=0.1
fg = p*dgamma(xx, shape=a, rate=b); lines(xx, fg, col=6, lwd=2)
f = f+fg; f1 = f1+fg; X = c(X, rgamma(round(p*n), shape=a, rate=b))

f1 = f1
lines(xx, f, col=1, lwd=2)
lines(xx, f1, col=1, lwd=2, lty=2)


n = length(X)
Rinit = 10
Gmax = 12

for (G in (2:Gmax)){
  cat('G = ', G, ': ')
  # Mixture
  Res = F_FitMixtGaussNull(X, G, 'E')
  logL = Res$logL
  cat('R = ')
  for (R in (1:Rinit))
    {
    cat(R, '')
    tau = matrix(0, n, G)
    for (i in (1:n)){
      prTmp = rgamma(G, shape=Res$tau[i,])
      prTmp = prTmp / sum(prTmp)
      tau[i,] = rmultinom(1, 1, Res$tau[i,])
      }
    #tau = t(rmultinom(n, 1, c(.5, rep(1/2/(G-1), G-1))))
    ResTmp = F_FitMixtGaussNull(X, G, 'E', tau)
    if (ResTmp$logL > logL){
      Res = ResTmp
      logL = Res$logL
      }
    }
  cat('\n')

  # Plot
  H = hist(X, breaks=ceiling(sqrt(n)), main='', xlab='', ylab='')
  w = mean(diff(H$mids))
  Xmin = min(X)
  Xmax = max(X)
  xx = seq(1.1*Xmin-0.1*Xmax, 1.1*Xmax-0.1*Xmin, length.out=500)
  f = rep(0, length(xx))
  for (g in (1:G))
  {
    fg = dnorm(xx, mean=Res$mean[g], sd=sqrt(Res$variance[g]))
    f = f + Res$prop[g]*fg
    lines(xx, n*w*Res$prop[g]*fg, lwd=2, col=(1+g))
  }
  lines(xx, n*w*f, lwd=2, col=1)
  
  #readline('...')
  
}


