# Fig for hidden metric space model
# Inspired from Handcock & al

library(lattice)
library(mnormt)
library(igraph)

# Parms
n = 20
prop = c(.5, .3, .2)
prob = .8*matrix(c(.9, .2, .5, .2, .1, .3, .5, .3, .8), 3, 3)
K = length(prop)
mu = matrix(c(0, sqrt(3)/2, sqrt(3)/2, 0, 1/sqrt(2), -1/sqrt(2)), 3, 2)
sigma = .2
grid.size = 100

####################################################
F_PlotNet <- function(Y, X, Col){
   # Y = position, X = adjacency, Col = colors
   dimY = dim(Y)[1]
   plot(Y, xlab='', ylab='', pch=20, 
        cex=1.5, col.axis=0, xaxt='n', yaxt='n')
   points(Y, pch=20, cex=2, col=Col)
   for(i in (1:(dimY-1))){
      for (j in ((i+1):dimY)){
         if (X[i, j]==1){lines(Y[c(i, j), 1], Y[c(i, j),2], lwd=1)
         }
      }
   }
#   points(Y, pch=20, cex=2, col=Col)
   if (dimY < 100){
      text(Y, label=(1:n), pos=3, offset=.3)
   }
}

####################################################

####################################################
# Hoff latent position model : LPM
####################################################
# Sim locations
ssigma = 3*sigma
Y = ssigma*matrix(rnorm(2*n), n, 2)
x.grid = seq(-3*ssigma, 3*ssigma, length.out=grid.size)
y.grid = x.grid
dY = matrix(0, length(x.grid), length(y.grid))
for (i in (1:length(x.grid))){
   for (j in (1:length(y.grid))){
      dY[i, j] = dmnorm(c(x.grid[i], y.grid[j]), varcov=diag(rep(ssigma^2, 2)))
   }
}
png('FigRSyst-LPM-density.png')
wireframe(dY, xlab='', ylab='', zlab='', drape=T); dev.off()
png('FigRSyst-LPM-contour.png')
contour(dY, xlab='', ylab='', col=(1:8), lwd=3); dev.off()

# Sim edges
D = as.matrix(dist(Y))
X = matrix(rbinom(n^2, 1, exp(-2*D)), n, n)
X[lower.tri(X, diag=T)] = 0

# # Plots
png('FigRSyst-LPM-Y.png')
F_PlotNet(Y, matrix(0, n, n), rep(6, n)); dev.off()
png('FigRSyst-LPM-XY.png')
F_PlotNet(Y, X, rep(6, n)); dev.off()
png('FigRSyst-LPM-X1.png')
F_PlotNet(cbind(runif(n), runif(n)), X, rep(6, n)); dev.off()
png('FigRSyst-LPM-X2.png')
F_PlotNet(cbind(runif(n), runif(n)), X, rep(6, n)); dev.off()
print(X[1:5, 1:5])

# # Larger network
# N = 5*n
# Y = ssigma*matrix(rnorm(2*N), N, 2)
# D = as.matrix(dist(Y))
# X = matrix(rbinom(N^2, 1, exp(-2*D)), N, N)
# X[lower.tri(X, diag=T)] = 0
# hist(colSums(X), breaks=sqrt(N))

####################################################
# Handcock latent position cluster model : LPCM
####################################################
# Sim locations
mmu = 1.5*mu
ssigma = 2*sigma
#ssigma = 1.5*sigma
x.grid = seq(min(mmu[, 1])-3*ssigma, max(mmu[, 1])+3*ssigma, length.out=grid.size)
y.grid = seq(min(mmu[, 2])-3*ssigma, max(mmu[, 2])+3*ssigma, length.out=grid.size)
dY = matrix(0, length(x.grid), length(y.grid))
for (i in (1:length(x.grid))){
   for (j in (1:length(y.grid))){
      for (k in (1:K)){
         dY[i, j] = dY[i, j] + prop[k]*dmnorm(c(x.grid[i], y.grid[j]), 
                                              mean = mmu[k, ], 
                                              varcov=diag(rep(ssigma^2, 2)))
      }
   }
}
png('FigRSyst-LPCM-density.png')
wireframe(dY, xlab='', ylab='', zlab='', drape=T);  dev.off()
png('FigRSyst-LPCM-contour.png')
contour(dY, xlab='', ylab='', col=(1:8), lwd=3); dev.off()

Z = as.vector((1:K) %*% rmultinom(n, 1, prop))
Y = mmu[Z, ] + ssigma*matrix(rnorm(2*n), n, 2)

# Sim edges
D = as.matrix(dist(Y))
X = matrix(rbinom(n^2, 1, exp(-2*D)), n, n)
X[lower.tri(X)] = 0
# Plots
png('FigRSyst-LPCM-Y.png')
F_PlotNet(Y, matrix(0, n, n), 1); dev.off()
png('FigRSyst-LPCM-YZ.png')
F_PlotNet(Y, matrix(0, n, n), 1+Z); dev.off()
png('FigRSyst-LPCM-XYZ.png')
F_PlotNet(Y, X, 1+Z); dev.off()
png('FigRSyst-LPCM-XY.png')
F_PlotNet(Y, X, rep(6, n)); dev.off()
# #png('FigRSyst-LPCM-X1.png')
# F_PlotNet(cbind(runif(n), runif(n)), X, rep(6, n)); #dev.off()
# #png('FigRSyst-LPCM-X2.png')
# F_PlotNet(cbind(runif(n), runif(n)), X, rep(6, n)); #dev.off()

####################################################
# Channarond latent position cluster model : KerNet
####################################################
# K = 3
prop = c(.5, .3, .2); K = length(prop); mu = matrix(c(0, sqrt(3)/2, sqrt(3)/2, 0, 1/sqrt(2), -1/sqrt(2)), K, 2); 
coef=5; mmu = coef*mu; ssigma = 4/3*coef*sigma
# K = 2
prop = c(.8, .2); K = length(prop); mu = matrix(c(sqrt(3)/2, sqrt(3)/2, 1/sqrt(2), -1/sqrt(2)), K, 2)
coef=5; mmu = coef*mu; ssigma = 5/3*coef*sigma
# K = 5
prop = c(.3, .2, .1, .3, .1); K = length(prop); 
mu = matrix(c(0, sqrt(3)/2, sqrt(3)/2, -sqrt(3)/2, -sqrt(3)/2, 0, 1/sqrt(2), -1/sqrt(2), 1/sqrt(2), -1/sqrt(2)), K, 2); 
coef=5; mmu = coef*mu; ssigma = 4/3*coef*sigma

N = 100*n
x.grid = seq(min(mmu[, 1])-3*ssigma, max(mmu[, 1])+3*ssigma, length.out=grid.size)
y.grid = seq(min(mmu[, 2])-3*ssigma, max(mmu[, 2])+3*ssigma, length.out=grid.size)
dY = matrix(0, length(x.grid), length(y.grid))
for (i in (1:length(x.grid))){
   for (j in (1:length(y.grid))){
      for (k in (1:K)){
         dY[i, j] = dY[i, j] + prop[k]*dmnorm(c(x.grid[i], y.grid[j]), 
                                              mean = mmu[k, ], 
                                              varcov=diag(rep(ssigma^2, 2)))
      }
   }
}
png(paste('FigRSyst-KerNet-K', K, '-density.png', sep=''))
wireframe(dY, xlab='', ylab='', zlab='', drape=T);  dev.off()
png(paste('FigRSyst-KerNet-K', K, '-contour.png', sep=''))
contour(x.grid, y.grid, dY, xlab='', ylab='', col=(1:8), lwd=3); dev.off()

#Z = as.vector((1:K) %*% rmultinom(N, 1, prop))
Z = as.vector((1:K) %*% rmultinom(N, 1, rep(1/K, K)))
Y = mmu[Z, ] + ssigma*matrix(rnorm(2*N), N, 2)
D = as.matrix(dist(Y))
Prob = exp(-2*D)
Prob[D > 2] = 0
X = matrix(rbinom(N^2, 1, Prob), N, N)
X[lower.tri(X, diag=T)] = 0
X = X + t(X)

TrueDensity = rep(0, N)
for (i in (1:N)){
   for (k in (1:K)){
      TrueDensity[i] = TrueDensity[i] + prop[k]*dmnorm(Y[i,] , 
                                                         mean = mmu[k, ], 
                                                         varcov=diag(rep(ssigma^2, 2)))
   }
}   

Deg = colSums(X) 
DegNorm = Deg/(N-1)
plot(TrueDensity, Deg, pch=20, cex=.5)
Order = order(Deg)
Y = Y[Order, ]
X = X[Order, Order]
Deg = Deg[Order]
DegNorm = DegNorm[Order]

png(paste('FigRSyst-KerNet-K', K, '-degree.png', sep=''))
hist(Deg, breaks=sqrt(N), main='', xlab='', ylab=''); dev.off()
png(paste('FigRSyst-KerNet-K', K, '-degnorm.png', sep=''))
hist(DegNorm, breaks=sqrt(N), main='', xlab='', ylab=''); dev.off()
png(paste('FigRSyst-KerNet-K', K, '.png', sep=''))
F_PlotNet(Y, X, 1); dev.off()

t.nb = 10; t.list = seq(0, max(DegNorm), length.out=t.nb)
t.list = unique(sort(DegNorm)); t.nb = length(t.list)
Ct = rep(0, t.nb); 
Nt = Ct
Xt = X
for (i in (1:length(t.list))){
   t = t.list[i]
   Nt[i] = sum(DegNorm >= t)
   cat(i, t, Nt[i])
   if (Nt[i] > 1){
      G = graph.adjacency(X[(DegNorm >= t), (DegNorm >= t)])
      Ct[i] = no.clusters(G, mode='weak')
      Xt = X
      Xt[(DegNorm < t), ] = 0
      Xt[, (DegNorm < t)] = 0
      if (i == 1){
         png(paste('FigRSyst-KerNet-K', K, '-t', i, '.png', sep=''))
         F_PlotNet(Y, Xt, 6*(DegNorm >= t)+8*(DegNorm < t)); dev.off()
#          plot(G); readline()
      }
      if((i > 1)){
         if(Ct[i] != Ct[i-1]){
#            png(paste('FigRSyst-KerNet-K', K, '-C', Ct[i], '.png', sep=''))
            png(paste('FigRSyst-KerNet-K', K, '-t', i, '.png', sep=''))
            F_PlotNet(Y, Xt, 6*(DegNorm >= t)+8*(DegNorm < t)); dev.off()
#             plot(G); readline()
         }
      }
      #readline()
#      points(Y[DegNorm < t, ], col=8, pch=20, cex=1.5)
      cat('', Ct[i])
   }
   cat('\n')
}
png(paste('FigRSyst-KerNet-K', K, '-C-N.png', sep=''))
plot(Nt/N, Ct, xlab='N(t)/n', ylab='C(t)', pch=20, cex=1.5, type='b'); dev.off()
for (i in (1:100)){dev.off()}
