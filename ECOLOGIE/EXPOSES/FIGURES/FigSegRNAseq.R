# Fig Ecole NGS -> réannotation RNAseq

rm(list=ls())

# fonctions
source('/media/donnees/RECHERCHE/RUPTURES/Meteo/SegCorrSpat/Prog/functions.segmentation.R')

# param
phi = 1e-5 # sur-dispersion

#DataName = 'RpBgeneYAL054C-428R1_L1'
#NGS = read.table(paste(DataName, ".txt", sep='')
#names(NGS) = c('position', 'signal')
#head(NGS)
#attach(NGS)
#Begin = 42600; End = 43800

# Data
#DataName = 'RpBgeneYAR002C-A-repdel'; Begin = 1500; End = 2800
DataName = 'geneYAL030W-ypdrep'; Begin = 1400; End = 2500
signal = read.table(paste(DataName, ".txt", sep=''), h=T)$x
position = (1:length(signal))
#Begin = 1; End = length(signal)

# truncature
signal[which.max(signal)] = 100

postscript(paste(DataName, "-",Begin, "-", End, "-raw.ps", sep=''), width=20, height=8, horizontal=F)
SubSet = which((position>Begin) & (position<End))
par(mai=rep(.4, 4))
plot(position[SubSet], log(1+signal[SubSet]), 
     type='b', xlab='', ylab='log(1+signal)', lwd=2)
dev.off()

# Cost matrix (Poisson)
Y = signal[SubSet]
n = length(Y)
cat('n=', n, '\n')
C = matrix(Inf, n, n)
cat('i=')
for (i in (1:(n-1)))
{
  if (100*floor(i / 100)==i){cat('', i)}
  Yk = Y[i]
  nk = 1
  for (j in ((i+1):n))
  {
    Yk = Yk + Y[j]
    nk = nk + 1
    mk = Yk/nk
#    if (mk > 0)    {C[i, j] = }
    if (mk > 0)    {C[i, j] = - Yk*log(mk) + nk*mk}
    if (mk == 0)    {C[i, j] = 0}    
  }
 }
cat('', n, '\n')
#image((1:n), (1:n), C)

# Segmentation
Kmax = 20
DP = DynProg(C, Kmax)
Klist = c(5, 7)
Kcol = c('red', 'blue')
postscript(paste(DataName, "-", Begin, "-", End, "-J.ps", sep=''), width=20, height=8, horizontal=F)
plot(DP$J.est, xlab='', ylab='', pch=20, cex=2, cex.axis=2)
dev.off()

Pen = (1:Kmax) * log(n/(1:Kmax))
Kcoef = c(50, 20)
#par(mfrow=c(2,1))
postscript(paste(DataName, "-", Begin, "-", End, "-J-K.ps", sep=''), width=20, height=8, horizontal=F)
plot(DP$J.est, xlab='', ylab='', pch=20, cex=2, cex.axis=2)
for (i in (1:length(Klist)))
{
  points(DP$J.est+Kcoef[i]*Pen, col=Kcol[i], pch=20, cex=2)
  abline(v = Klist[i], col=Kcol[i], lwd=3, lty=2)
}
dev.off()

postscript(paste(DataName, "-", Begin, "-", End, "-seg.ps", sep=''), width=20, height=8, horizontal=F)
plot((1:n), log(1+Y), type='p', xlab='', ylab='', lwd=2, pch=20, cex=1.5)
for (i in (1:length(Klist)))
{
  K = Klist[i]
  tau = c(0, DP$t.est[K,1:K])
  for (k in (1:K))
  {
    Ik = ((tau[k]+1):tau[k+1])
    nk = tau[k+1] - tau[k]
    lines(Ik, rep(log(1+mean(Y[Ik])), nk), col=Kcol[i], lwd=4)
    abline(v= tau[k+1]+.5, col=Kcol[i], lwd=3, lty=2)
  }
}
dev.off()
