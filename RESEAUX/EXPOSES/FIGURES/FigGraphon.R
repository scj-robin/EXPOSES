# Fig for graphon model

library(lattice)
library(mixer);
source("/home/robin/Bureau/BMA/Dropbox-vbemapp/PgmSim/functionsGraphon.R")
source("/home/robin/Bureau/BMA/Dropbox-vbemapp/PgmSim/functionsGraphonInference.R")
source("/home/robin/Bureau/BMA/Dropbox-vbemapp/PgmSim/functionsMotifsInference.R")

# Parms
n = 5e2
prop = c(.5, .3, .2)
K = length(prop)
prob = matrix(c(.1, .2, .5, .2, .4, .7, .5, .7, .9), K, K)

####################################################
# SBM
####################################################
# Sim groups
Zmat = rmultinom(n, 1, prop)
Z = as.vector((1:K) %*% Zmat)
Y = cbind(runif(n), runif(n))
# Sim edges
X = matrix(rbinom(n^2, 1, prob[Z, Z]), n, n)
X = X - diag(diag(X))
X[lower.tri(X)] = 0
X = X + t(X)
# # Plots
# png('FigGraphon-SBM-Z.png')
# F_PlotNet(Y, matrix(0, n, n), 1+Z); dev.off()
# png('FigGraphon-SBM-XZ.png')
# F_PlotNet(Y, X, 1+Z); dev.off()
# png('FigGraphon-SBM-X1.png')
# F_PlotNet(Y, X, rep(6, n)); dev.off()
# png('FigGraphon-SBM-X2.png')
# F_PlotNet(cbind(runif(n), runif(n)), X, rep(6, n)); dev.off()
h = 1e-2
u = seq(0, 1, by=h)
grid = (u) %*% t(u)
size = length(u)
cumprop = cumsum(prop)
gammaSBM = matrix(0, size, size)
for (i in (1:size)){
   for (j in (1:size)){
      k = 1+sum(cumprop < u[i])
      l = 1+sum(cumprop < u[j])
      gammaSBM[i, j] = prob[k, l]
   }
}
png('FigGraphon-SBM-graphon.png')
wireframe(gammaSBM, xlab='k', ylab="l", zlab='gamma', drape=T)
dev.off()

pdf('FigGraphon-SBM-graphon-alpha.pdf')
# wireframe(gammaSBM, xlab='k', ylab="l", zlab='alpha')
trellis.par.set("axis.line",list(col=NA,lty=1,lwd=1))
wireframe(qlogis(gammaSBM), shade = TRUE, light.source = c(3,3,3), aspect = c(61/87, 0.4), 
          row.values = u, column.values = u, 
          xlab='', ylab='', zlab='', default.scales=list(arrows=F))
dev.off()

pdf('FigGraphon-SBM-graphon-alpha-largeticks.pdf')
# wireframe(gammaSBM, xlab='k', ylab="l", zlab='alpha')
trellis.par.set("axis.line",list(col=NA,lty=1,lwd=1))
wireframe(qlogis(gammaSBM), shade = TRUE, light.source = c(3,3,3), aspect = c(61/87, 0.4), 
          row.values = u, column.values = u, 
          xlab='', ylab='', zlab='', default.scales=list(arrows=F), scales=list(cex=1.5))
dev.off()

####################################################
# Averaged SBM graphon
####################################################
S = 1e2
p0 = rep(1, K)
Nemp = rowSums(Zmat)
g0 = matrix(1, K, K)
g1 = matrix(1, K, K)
M1emp = Zmat %*% X %*% t(Zmat)
M0emp = Zmat %*% (1-X) %*% t(Zmat)
pave = rep(0, K)
gave = matrix(0, K, K)
surfave = matrix(0, size, size)
for (s in (1:S)){
   psim = rep(0, K)
   gsim = matrix(0, K, K)
   for (k in (1:K)){
      psim[k] = rgamma(1, p0+Nemp, 1)
      for (l in (k:K)){
         gsim[k, l] = rbeta(1, g1[k, l]+M1emp[k, l], g0[k, l]+M0emp[k, l])
      }
   }
   psim = psim / sum(psim)
   surfave = surfave + surfSBM(size, psim, gsim)
}
png('FigGraphon-SBM-average.png')
wireframe(surfave, xlab='z', ylab="z'", zlab='gamma', drape=T)
dev.off()

####################################################
# VBEM inference
####################################################
Qmin = 1; Qmax = K; nbrepeat = 10; L = 100
# Loop of repeated VBEM
lout<-list()
tab_criteria<-matrix(0, nbrepeat, Qmax)
for(start in 1:nbrepeat) {
   cat(start, '')
   xout<-mixer(X, qmin=1, qmax=Qmax, method="bayesian", verbose=F)
   lout[[start]]<-xout
   tab_criteria[start, ]<-unlist(lapply(xout$output, function(x) x$criterion))
}

# Selection of best starts for each Q
bestStart = apply(tab_criteria, 2, which.max)
bestCrit = apply(tab_criteria, 2, max)
SBM = list()
for (Q in Qmin:Qmax){
   SBM[[Q]] = lout[[bestStart[Q]]]$output[Q][[1]]
   ord = order(SBM[[Q]]$alphas %*% SBM[[Q]]$Pis)
   SBM[[Q]]$alphas = SBM[[Q]]$alphas[ord]
   SBM[[Q]]$a = SBM[[Q]]$a[ord]
   SBM[[Q]]$Pis = SBM[[Q]]$Pis[ord, ord]
   SBM[[Q]]$eta = SBM[[Q]]$eta[ord, ord]
   SBM[[Q]]$zeta = SBM[[Q]]$zeta[ord, ord]
   SBM[[Q]]$Taus = SBM[[Q]]$Taus[ord, ]
}
phiBestMAP <- surfSBM(L, SBM[[K]]$alphas, SBM[[K]]$Pis)
png('FigGraphon-SBM-VBEM.png')
wireframe(phiBestMAP, xlab='k', ylab="l", zlab='gamma', drape=T)
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
SBMs = smooth.spline(grid, gammaSBM, spar=1e-6)
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
png('FigGraphon-W-graphon.png')
wireframe(gamma, xlab='z', ylab="z'", drape=T)
dev.off()
