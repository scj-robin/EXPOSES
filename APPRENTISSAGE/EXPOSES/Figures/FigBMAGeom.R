source('/media/donnees/RECHERCHE/APPRENTISSAGE/PGM/F_FitMixtGeom.R')
DirFig = '/media/donnees/RECHERCHE/ECOLOGIE/EXPOSES/FIGURES/'

###########################################################################
# Simul
library(triangle)
n = 1000
Gtrue = 10
Support = rgeom(Gtrue, .05)
Ptrue = rgamma(Gtrue, rep(1, Gtrue), rep(1, Gtrue))
Ptrue = Ptrue/sum(Ptrue)
Z = (1:Gtrue) %*% rmultinom(n, 1, Ptrue)
X = rep(0, n)
for (i in (1:n))
{
  X[i] = rtriangle(1, a=0, c=0, b=Support[Z[i]])
}
X = floor(X)

###########################################################################
# Tap data
tap = read.csv2('~/Dropbox/ConvexDiscret/articleEcolo/Notes/tap_et_al/tap_et_al.data.csv')
X <- tap$poids
XX = matrix(0, max(X), 2)
for (j in (1:max(X))){
  XX[j, 1] = j
  XX[j, 2] = sum(X==j)
}
XX = as.data.frame(XX)
names(XX) = c('j', 'n_j')
DataName = 'Tap'; coefy = 4

###########################################################################
# Species data
library(SPECIES)
#DataName = 'Butterfly'; data(butterfly); XX = butterfly; XX = XX[XX$j < 25,]; coefy = 3
#DataName = 'Cottontail'; data(cottontail); XX = cottontail; coefy = 3
#DataName = 'EST'; data(EST); XX = EST; XX = as.data.frame(XX); coefy = 6
#DataName = 'Traffic'; data(traffic); XX = traffic; XX = as.data.frame(XX); coefy = 3
#DataName = 'Insects'; data(insects); XX = insects; coefy = 3
#DataName = 'Microbial'; data(microbial); XX = microbial; names(XX) = c('j', 'n_j'); XX = XX[XX$j < 25,]; coefy = 3
X = c()
for(i in (1:dim(XX)[1])){X = c(X, rep(XX$j[i], XX$n_j[i]))}
X = X-1

###########################################################################
# Fit mixture
n = length(X)
Rinit = 10
Gmax = 5
#par(mfrow=c(ceiling(sqrt(Gmax)), ceiling(sqrt(Gmax))))
p0 = rep(0, Gmax); logL = p0

postscript(paste(DirFig, DataName, '.eps' ,sep=''))
xx = (0:(max(XX$j)))
xlimsup = max(XX$j)
xlimsup = 50
ylimsup = coefy*max(XX$n_j)
plot(XX$j, XX$n_j, pch=20, xlim=c(0, xlimsup), ylim=c(0, ylimsup), 
     xlab='', ylab='', main='', cex=2, cex.axis=2)
dev.off()

for (G in (1:Gmax)){
  cat('G =', G)
  # Mixture
  Res = F_FitMixtGeom(X, G)
  logL[G] = Res$logL
  
  cat(': R = ')
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
    ResTmp = F_FitMixtGeom(X, G, tau)
    if (ResTmp$logL > logL){
      Res = ResTmp
      logL = Res$logL
      }
    }
  cat('\n')
  
  # Pred N0
  p0g = Res$prob / (1 - Res$prob)
  p0[G] = sum(Res$prop * p0g)

  # Plot
  postscript(paste(DirFig, DataName, '-MixtGeom-K', G, '.eps' ,sep=''))
  xx = (0:(max(XX$j)))
  plot(XX$j, XX$n_j, pch=20, xlim=c(0, xlimsup), ylim=c(0, ylimsup), 
       xlab='', ylab='', main='', cex=2, cex.axis=2)
  #hist(X, breaks=(-.5:(max(X)+.5)), main='', xlab='', ylab='')
  fg = matrix(NaN, length(xx), G)
  for (g in (1:G))
  {
    fg[(2:length(xx)),g] = dgeom(xx[xx>0], prob=Res$prob[g]) / (1-Res$prob[g])
    lines(xx, n*Res$prop[g]*fg[,g], lwd=4, col='blue')
    lines(xx[1:2], n*Res$prop[g]*c(p0g[g], fg[2, g]), lwd=4, lty=2, col='blue')
  }
  f = fg %*% Res$prop
  lines(xx, n*f, lwd=4, col=1)
    
  lines(xx[1:2], n*c(p0[G], f[2]), lwd=4, lty=2, col='red')
  points(0, n*p0[G], col='red', cex=2)
  dev.off()
  
  #readline('...')
  
}

###########################################################################
BIC = logL - .5*log(n)*(2*(1:Gmax)-1)
plot(BIC)
plot(p0)
