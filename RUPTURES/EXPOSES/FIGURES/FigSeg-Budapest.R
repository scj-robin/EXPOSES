# Figure segmentation -> Budapest
# rm(list=ls())
Dir = '/home/robin/RECHERCHE/RUPTURES/EXPOSES/FIGURES/'

# Parms
n = 100
K = 5
cex = 1.5
Export = T

# Simul
s = s+1
# s = 19, 31, 49, 58, 70, 90
s = 90
set.seed(s)
m = 20+4*rnorm(K)
t = c(0, sort(sample((1:(n-1)), (K-1))), n)
while(min(diff(t))<2){t = c(0, sort(sample((1:(n-1)), (K-1))), n)}
nk = diff(t)
mu = rep(m[1], nk[1])
for (k in (2:K))
  {mu = c(mu, rep(m[k], nk[k]))}
E = rnorm(n)
Y = mu+E

# plot
if (Export==T){
   pdf(file=paste(Dir, "FigSeg-Budapest-0.pdf", sep=""))
}
plot((1:n), Y, col=1, pch=20, cex=cex, xlab='time', ylab='signal')
if (Export==T){dev.off()}

if (Export==T){
   pdf(file=paste(Dir, "FigSeg-Budapest-1.pdf", sep=""))
}
plot((1:n), Y, col=0, xlab='time', ylab='signal')
if (Export==T){dev.off()}

if (Export==T){
   pdf(file=paste(Dir, "FigSeg-Budapest-2.pdf", sep=""))
}
plot((1:n), Y, col=0, xlab='time', ylab='signal')
abline(v=t+.5, lty=2, lwd=5, col=4)
if (Export==T){dev.off()}

if (Export==T){
   pdf(file=paste(Dir, "FigSeg-Budapest-3.pdf", sep=""))
}
plot((1:n), Y, col=0, xlab='time', ylab='signal', col.axis=0)
abline(v=t+.5, lty=2, lwd=5, col=4)
for (k in (1:K))
  {lines(((t[k]+1):t[k+1]), rep(m[k], nk[k]), lwd=5, col=2)}
if (Export==T){dev.off()}

if (Export==T){
   pdf(file=paste(Dir, "FigSeg-Budapest-4.pdf", sep=""))
}
plot((1:n), Y, col=0, xlab='time', ylab='signal')
abline(v=t+.5, lty=2, lwd=5, col=4)
for (k in (1:K))
  {lines(((t[k]+1):t[k+1]), rep(m[k], nk[k]), lwd=5, col=2)}
points((1:n), Y, pch=20, col=1, cex=cex)
if (Export==T){dev.off()}
