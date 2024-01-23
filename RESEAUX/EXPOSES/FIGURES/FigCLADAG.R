# Fig for hidden metric space model
# Inspired from Handcock & al

library(lattice)

# Parms
n = 20
prop = c(.5, .3, .2)
prob = .8*matrix(c(.9, .2, .5, .2, .1, .3, .5, .3, .8), 3, 3)
K = length(prop)
mu = matrix(c(0, sqrt(3)/2, sqrt(3)/2, 0, 1/sqrt(2), -1/sqrt(2)), 3, 2)
sigma = .2
kappa = 3 # distance decrease

####################################################
F_PlotNet <- function(Y, X, Col){
   # Y = position, X = adjacency, Col = colors
   n = dim(Y)[1]
   plot(Y, xlab='', ylab='', pch=20, 
        cex=1.5, col.axis=0, xaxt='n', yaxt='n')
   for(i in (1:(n-1))){
      for (j in ((i+1):n)){
         if (X[i, j]==1){lines(Y[c(i, j), 1], Y[c(i, j),2], lwd=2)
         }
      }
   }
   points(Y, pch=20, cex=2, col=Col)
   text(Y, label=(1:n), pos=3, offset=.3)
}
####################################################

####################################################
# Hoff latent position model : LPM
####################################################
# Sim locations
Y = 3*sigma*matrix(rnorm(2*n), n, 2)
# Sim edges
D = as.matrix(dist(Y))
X = matrix(rbinom(n^2, 1, exp(-2*D)), n, n)
X[lower.tri(X, diag=T)] = 0

# # Plots
png('FigCLADAG-LPM-Y.png')
F_PlotNet(Y, matrix(0, n, n), rep(6, n)); dev.off()
png('FigCLADAG-LPM-XY.png')
F_PlotNet(Y, X, rep(6, n)); dev.off()
png('FigCLADAG-LPM-X1.png')
F_PlotNet(cbind(runif(n), runif(n)), X, rep(6, n)); dev.off()
png('FigCLADAG-LPM-X2.png')
F_PlotNet(cbind(runif(n), runif(n)), X, rep(6, n)); dev.off()
print(X[1:5, 1:5])

####################################################
# Handcock latent position cluster model : LPCM
####################################################
# Sim locations
Z = as.vector((1:K) %*% rmultinom(n, 1, prop))
Y = mu[Z, ] + sigma*matrix(rnorm(2*n), n, 2)
# Sim edges
D = as.matrix(dist(Y))
X = matrix(rbinom(n^2, 1, exp(-2*D)), n, n)
X[lower.tri(X)] = 0
# # Plots
# png('FigCLADAG-LPCM-YZ.png')
# F_PlotNet(Y, matrix(0, n, n), 1+Z); dev.off()
# png('FigCLADAG-LPCM-XYZ.png')
# F_PlotNet(Y, X, 1+Z); dev.off()
# png('FigCLADAG-LPCM-XY.png')
# F_PlotNet(Y, X, rep(6, n)); dev.off()
# png('FigCLADAG-LPCM-X1.png')
# F_PlotNet(cbind(runif(n), runif(n)), X, rep(6, n)); dev.off()
# png('FigCLADAG-LPCM-X2.png')
# F_PlotNet(cbind(runif(n), runif(n)), X, rep(6, n)); dev.off()

####################################################
# SBM
####################################################
# Sim groups
Z = as.vector((1:K) %*% rmultinom(n, 1, prop))
Y = cbind(runif(n), runif(n))
# Sim edges
X = matrix(rbinom(n^2, 1, prob[Z, Z]), n, n)
X[lower.tri(X)] = 0
# # Plots
# png('FigCLADAG-SBM-Z.png')
# F_PlotNet(Y, matrix(0, n, n), 1+Z); dev.off()
# png('FigCLADAG-SBM-XZ.png')
# F_PlotNet(Y, X, 1+Z); dev.off()
# png('FigCLADAG-SBM-X1.png')
# F_PlotNet(Y, X, rep(6, n)); dev.off()
# png('FigCLADAG-SBM-X2.png')
# F_PlotNet(cbind(runif(n), runif(n)), X, rep(6, n)); dev.off()
h = 1e-2
u = seq(0, 1, by=h)
grid = (u) %*% t(u)
size = length(u)
cumprop = cumsum(prop)
gammaSBM = matrix(0, size, size)
prob = matrix(c(.1, .2, .5, .2, .4, .7, .5, .7, .9), K, K)
for (i in (1:size)){
   for (j in (1:size)){
      k = 1+sum(cumprop < u[i])
      l = 1+sum(cumprop < u[j])
      gammaSBM[i, j] = prob[k, l]
   }
}
png('FigCLADAG-SBM-graphon.png')
wireframe(gammaSBM, xlab='k', ylab="l", zlab='gamma', drape=T)
dev.off()

####################################################
# Smoothed SBM graphon
####################################################
w = 5 # smoothing
W = (-w:w) %*% t(-w:w)
W = exp(-abs(W))
W = W / sum(W)
gammaSBMs = gammaSBM
for (i in (1:size)){
   for (j in (1:size)){
      gammaSBMs[i, j] = sum(W * gammaSBM[(i-w):(i+w), (j-w):(j+w)])   
   }
}
SBMs = smooth.spline(grid, gammaSBM, spar=.1)
gammaSBMs = predict(SBMs, grid)$y
wireframe(gammaSBMs, xlab='k', ylab="l", zlab='gamma', drape=T)

####################################################
# Graphon
####################################################
gamma = matrix(0, size, size)
rho = .1
lambda = 4
for (i in (1:size)){
   for (j in (1:size)){
      gamma[i, j] = rho*lambda^2*(u[i]*u[j])^(lambda-1)
   }
}
png('FigCLADAG-W-graphon.png')
wireframe(gamma, xlab='z', ylab="z'", drape=T)
dev.off()
png('FigCLADAG-W-graphon2.png')
wireframe(gamma, xlab='u', ylab="u'", zlab='w', drape=T)
dev.off()

######################################################################################
# Universal thresholding
load('~/Dropbox/vbemapp/gof-network/pgm/exemples/Blog.Rdata')
Y = Blog$Net; SVD = svd(Y); n = ncol(Y); 
plot(sort(SVD$d), pch=20, type='b'); 
t = 2.01*sqrt(n); 
t = 2.01*sqrt(sum(Y/2)/n); 
# t = 7.5
abline(h=t, col=2, lwd=2)
sel = which(SVD$d >= t)
P.hat = SVD$u[, sel] %*% diag(SVD$d[sel]) %*% t(SVD$v[, sel])
D.hat = rowSums(P.hat)
library(lattice)
plot(D.hat[order(D.hat)], pch=20, type='b')
image(Y[order(D.hat), order(D.hat)])
image(P.hat[order(D.hat), order(D.hat)])
wireframe(P.hat[order(D.hat), order(D.hat)], xlab='u', ylab="u'", zlab='w', drape=T)
