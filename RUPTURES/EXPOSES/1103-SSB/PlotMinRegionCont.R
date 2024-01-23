# Plot for 'MinRegionCont' slides

rm(list=ls())

# Sim 1 profile
n = 10
lambda = 5
mu = 75
L = -log(runif(n))/lambda
M = -log(runif(n))/mu
D = c(0, cumsum(as.vector(t(cbind(L, M)))))
#postscript('../Figures/SingleProfile.ps')
plot(D, rep(0, (2*n+1)), col=0, xlab='', ylab='', xlim=c(0, floor(10*max(D))/10), axes=F)
for (i  in (1:(2*n)))
  {
  s =  .5*(1+(-1)^i)
  lines(c(D[i], D[i+1]), rep(s, 2), col=3-s, lwd=5)
  lines(c(D[i+1], D[i+1]), c(0, 1), col=1, lwd=2, lty=2)
  }
D
#dev.off()

#postscript('../Figures/SingleProfileBis.ps')
plot(D, rep(0, (2*n+1)), col=0, xlab='', ylab='', xlim=c(0, floor(10*max(D))/10), axes=F)
for (i  in (1:(2*n)))
  {
  s =  .5*(1+(-1)^i)
  lines(c(D[i], D[i+1]), rep(0, 2), col=3-s, lwd=5)
  #lines(c(D[i+1], D[i+1]), c(0, 1), col=1, lwd=2, lty=2)
  }
D
#dev.off()

# Sim X_+
m = 100
M = ceiling(sqrt(m))
N = 300
D = -log(runif(N))
X = rep(0, length(D))
Time = X
for (i in (1:length(D)))
  {
  Lambda = (m-X[i])*lambda
  Mu = X[i]*mu
  U = runif(1)
  if (U < Lambda/(Lambda+Mu))
    {
    X[i+1] = X[i]+1
    Time[i+1] = Time[i] + D[i]/Lambda
    }
  if (U > Lambda/(Lambda+Mu))
    {
    X[i+1] = X[i]-1
    Time[i+1] = Time[i] + D[i]/mu
    }
  }
#write.table(cbind(Time, X), file='SimXplus.txt')
Sim = read.table('SimXplus.txt')
X = Sim$X
Time = Sim$Time
postscript('../Figures/SumProcess1.ps')
plot(Time, X, type='s', ylab='X+(t)', xlab='t', xlim=c(0, 1.1), lwd=2)
dev.off()

postscript('../Figures/SumProcess2.ps')
plot(Time, X, type='s', ylab='X+(t)', xlab='t', xlim=c(0, 1.1), lwd=2)
XM = which(X==M)
XM = XM[1:5]
XM = XM[-which(X[XM-1]>M)]
XM_1 = rep(0, length(XM))
for (j in (1:length(XM)))
  {
  XM_1[j] = XM[j] + min(which(X[(XM[j]):length(Time)]==(M-1))) -1
  }
abline(h=M, lwd=3, col=2, lty=1)
for (j in XM[1])
  {
  abline(v=Time[j], lwd=4, col=4, lty=2)
  }
dev.off()

postscript('../Figures/SumProcess3.ps')
l = .05
plot(Time, X, type='s', ylab='X+(t)', xlab='t', xlim=c(0, 1.1), lwd=2)
abline(h=M, lwd=3, col=2, lty=1)
abline(v=Time[XM[1]], lwd=4, col=4, lty=2)
abline(v=Time[XM_1[1]], lwd=4, col=5, lty=2)
lines(c(Time[XM[1]], Time[XM[1]]+l), c(M, M), lwd=4, col=4)
dev.off()

postscript('../Figures/SumProcess4.ps')
plot(Time, X, type='s', ylab='X+(t)', xlab='t', xlim=c(0, 1.1), lwd=2)
abline(h=M, lwd=3, col=2, lty=1)
for (j in c(1, 2))
  {
  abline(v=Time[XM[j]], lwd=4, col=4, lty=2)
  }
abline(v=Time[XM_1[1]], lwd=4, col=5, lty=2)
dev.off()

postscript('../Figures/SumProcess5.ps')
plot(Time, X, type='s', ylab='X+(t)', xlab='t', xlim=c(0, 1.1), lwd=2)
abline(h=M, lwd=3, col=2, lty=1)
for (j in (1:length(XM)))
  {
  abline(v=Time[XM[j]], lwd=4, col=4, lty=2)
  abline(v=Time[XM_1[j]], lwd=4, col=5, lty=2)
  lines(c(Time[XM[j]], Time[XM[j]]+l), c(M, M), lwd=4, col=4)
  }
dev.off()


