# Figure segmentation -> Budapest
rm(list=ls())

Dir = '/media/donnees/RECHERCHE/RUPTURES/Exposes/Figures/'

# Parms
n = 10
S = 3
K = 2
Export = T
gamma = 1

# Simul
U = gamma*rnorm(n)
m = matrix(0, S, K)
nk = matrix(0, S, K)
t = matrix(0, S, (K+1))
mu = matrix(0, S, n)
E = mu; Y = mu
for (s in (1:S))
  {
  m[s,] = 10+2*rnorm(K)
  t[s, ] = c(0, sort(sample((1:(n-1)), (K-1))), n)
  nk[s, ] = diff(t[s,])
  mutmp = rep(m[s, 1], nk[s, 1])
  for (k in (2:K))
    {mutmp = c(mutmp, rep(m[s, k], nk[s, k]))}
  mu[s,] = mutmp
  E[s,] = rnorm(n)
  Y[s,] = mu[s,]+U+E[s,]
  }

# plot
setEPS()
if (Export==T){postscript(file=paste(Dir, "FigSegMult-Budapest-1.eps", sep=""))}
par(mfrow=c((S+1),1))
plot((1:n), U, col=3, xlab='', ylab='', pch=20, cex=2)
#for (s in (1:S))
#  {
#  plot((1:n), Y[s, ], pch=20, cex=1.5, col=0, xlab='', ylab='', col.axis=0)
#}
if (Export==T){dev.off()}

if (Export==T){postscript(file=paste(Dir, "FigSegMult-Budapest-2.eps", sep=""))}
par(mfrow=c((S+1),1))
plot((1:n), U, col=3, xlab='', ylab='', pch=20, cex=2)
for (s in (1:S))
  {
  plot((1:n), Y[s, ], pch=20, cex=1.5, col=0, xlab='', ylab='')
  for (k in (1:K))
    {lines(c((t[s, k]+.5), (t[s, k+1]+.5)), rep(m[s, k], 2), lwd=3, col=2)}
}
if (Export==T){dev.off()}

if (Export==T){postscript(file=paste(Dir, "FigSegMult-Budapest-3.eps", sep=""))}
par(mfrow=c((S+1),1))
plot((1:n), U, col=3, xlab='', ylab='', pch=20, cex=2)
for (s in (1:S))
  {
  plot((1:n), Y[s, ], pch=20, cex=2, col=0, xlab='', ylab='')
  for (k in (1:K))
    {lines(c((t[s, k]+.5), (t[s, k+1]+.5)), rep(m[s, k], 2), lwd=3, col=2)}
  points((1:n), mu[s,] + U, pch=20, col=3, cex=2)
}
if (Export==T){dev.off()}
    
if (Export==T){postscript(file=paste(Dir, "FigSegMult-Budapest-4.eps", sep=""))}
par(mfrow=c((S+1),1))
plot((1:n), U, col=3, xlab='', ylab='', pch=20, cex=2)
for (s in (1:S))
  {
  plot((1:n), Y[s, ], pch=20, cex=1.5, col=0, xlab='', ylab='')
  for (k in (1:K))
    {lines(c((t[s, k]+.5), (t[s, k+1]+.5)), rep(m[s, k], 2), lwd=3, col=2)}
  points((1:n), mu[s,] + U, pch=20, col=3, cex=2)
  points((1:n), Y[s, ], pch=20, cex=2, col=1, xlab='', ylab='')
}
if (Export==T){dev.off()}
