# Figures segmentation -> Cambridge
rm(list=ls())

# Dir = '/media/donnees/RECHERCHE/RUPTURES/Exposes/Figures/'
Export = F
library(cghseg)
library(Segmentor3IsBack)

# # Simul CGH
# n = 5e2; K = round(log(n)); sigma = 0.5
# m = c(2, 3, 1, 2, 4, 1) #m = rpois(K, 2) 
# t = c(0, 251, 313, 354, 369, 434, 500) #t = c(0, sort(sample((1:(n-1)), (K-1))), n)
# nk = diff(t)
# mu = rep(m, times=nk)
# E = sigma*rnorm(n)
# Y = mu+E

# # Plot CGH
# pdf('FigSeg-Cambridge-CGH.pdf', width=10, height=5)
# Y = mu+E
# plot(Y, pch=20, col=1, xlab='', ylab='')
# dev.off()
# # lines(mu, col='blue', lwd=2)
# CGHo = new('CGHoptions')
# CGHd = new("CGHdata", Y = Y)
# CGHr  = uniseg(CGHd,CGHo)
# t.hat = which(getbp(CGHr)$Y[, 2] ==1)
# mu.hat = getsegprofiles(CGHr)[, 1]; 
# pdf('FigSeg-Cambridge-CGH-seg.pdf', width=10, height=5)
# Y = mu+E
# plot(Y, pch=20, col=1, xlab='', ylab='')
# abline(v=t.hat, col='red', lty=2)
# lines(mu.hat, lwd=3, col='red')
# dev.off()

# Simul RNAseq
n = 1e5; K = n/1e3; sigma = .75
phi = .5
odd = ((1:K)%%2 == 1)
p = rep(c(.99, .33), ceiling(K/2)); p = p[1:K]
t = c(0, sort(sample((1:(n-1)), (K-1))), n)
#t = c(0, cumsum(rpois(K, n/K)))
nk = diff(t); 
nk[odd] = 4*nk[odd]
t = c(0, cumsum(nk)); n = max(t)
pi = rep(p, times=nk)
Y = rnbinom(n, size=rep(phi, n), prob=pi)

# Plot RNAseq
pdf('FigSeg-Cambridge-RNAseq.pdf', width=10, height=5)
plot(Y, pch=20, cex=.5, type='l', col=1, xlab='', ylab='')
dev.off()
SEG = Segmentor(Y, model=3, keep=T, Kmax = 2*K)
K.hat = SelectModel(SEG)
bestSEG = BestSegmentation(SEG, K.hat)
t.hat = c(0, getBreaks(SEG)[K.hat,1:(K.hat-1)])
pdf('FigSeg-Cambridge-RNAseq-seg.pdf', width=10, height=5)
Y = getData(SEG)
plot(Y, pch=20, cex=.5, type='l', col=1, xlab='', ylab='')
abline(v=t.hat, col='red', lwd=.5, lty=2)
lines(Y, pch=20, cex=.5, col=1)
dev.off()
