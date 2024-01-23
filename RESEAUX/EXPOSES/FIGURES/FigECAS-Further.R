# Figures for ECAS / SFdS / NAIM : reminder on prob. dist

library(sna); library(mvtnorm); library(huge); library(LITree); library(ROCR); 
library(mixtools); library(latex2exp)
source('FigECAS-Functions.R')
source('/home/robin/PgmR/General/FunctionsMatVec.R')
source('/home/robin/RECHERCHE/RESEAUX/S-Founas/NoisyNetworkInference/Pgm/Functions/SimulFunctions.R')
source('/home/robin/RECHERCHE/RESEAUX/S-Founas/NoisyNetworkInference/Pgm/Functions/ResultFunctions.R')

# Missing actor
p = 5; G = matrix(c(0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0), p, p)
pdf('FigECAS-GGM-5nodes.pdf'); par(mex=0.01)
pos = gplot(G, gmode='graph', label=1:p, label.pos=5, vertex.col=c(0, 0, 8, 0, 0), vertex.cex=2, label.cex=2)
dev.off()

# Missing actor : node 3 missing
iMiss = 3; Gmiss = G; 
sapply(1:(p-1), function(j){sapply((j+1):p, function(k){if(G[iMiss, j]*G[iMiss, k]==1){Gmiss[j, k] <<- Gmiss[k, j] <<- 1}})})
Gmiss = Gmiss[-iMiss, ]; Gmiss = Gmiss[, -iMiss]
pdf(paste0('FigECAS-GGM-5nodes-', iMiss, 'missing.pdf')); par(mex=0.01)
gplot(Gmiss, gmode='graph', label=(1:p)[-iMiss], label.pos=5, vertex.col=0, vertex.cex=2, label.cex=2, coord=pos[-iMiss, ])
dev.off()

# Complete covariance
Gsign = F_Vec2Sym(F_Sym2Vec(G * matrix(2*rbinom(p^2, 1, .5)-1, p, p)))
lambda = 1.1; Omega = lambda*diag(rowSums(G)) - Gsign
while(min(eigen(Omega)$values) < 1e-10){lambda = 1.1*lambda; Omega = lambda*diag(rowSums(G)) - G}
sigma = 2*rgamma(p, 1, 1); corSigma = cov2cor(solve(Omega)); Sigma = diag(sigma)%*%corSigma%*%diag(sigma)
Omega = solve(Sigma); Omega = Omega * (abs(Omega) > 1e-10); corOmega = -cov2cor(Omega); diag(corOmega) = 1

# Marginal covariance
OmegaMiss = Omega[-iMiss, -iMiss] - Omega[-iMiss, iMiss]%o%Omega[iMiss, -iMiss]/Omega[iMiss, iMiss]
SigmaMiss = solve(OmegaMiss); corSigmaMiss = cov2cor(SigmaMiss); 
corOmegaMiss = -cov2cor(OmegaMiss); diag(corOmegaMiss) = 1

par(mfrow=c(2, 2))
image(corSigma, axes=FALSE); abline(h=(iMiss)/(p+1), v=(iMiss)/(p+1)); image(corOmega, axes=FALSE); abline(h=(iMiss)/(p+1), v=(iMiss)/(p+1))
image(corSigmaMiss, axes=FALSE); image(corOmegaMiss, axes=FALSE)

corSigma[lower.tri(corSigma, diag=TRUE)] = NA; Ftab(round(corSigma, 2))
corOmega[lower.tri(corOmega, diag=TRUE)] = NA; Ftab(round(corOmega, 2))
corSigmaMiss[lower.tri(corSigmaMiss, diag=TRUE)] = NA; Ftab(round(corSigmaMiss, 2))
corOmegaMiss[lower.tri(corOmegaMiss, diag=TRUE)] = NA; Ftab(round(corOmegaMiss, 2))

# Simulated GGM 4 nodes
seed = 3; set.seed(seed)
n = 25; Y = rmvnorm(n, sigma=Sigma)
SigmaHat = cov(Y); corSigmaHat = cor(Y)
OmegaHat = solve(SigmaHat); corOmegaHat = -cov2cor(OmegaHat); diag(corOmegaHat) = 1
Ftab(cbind(round(OmegaHat, 2))); 
# Ftab(round(corOmegaHat, 2))
Ftab(cbind(round(SigmaHat, 2))); 
# Ftab(round(corSigmaHat, 2))

# L0, L1, L2
alphaList = c(.8, .6, .4)
beta = c(1.5, 3); Sigma = matrix(c(1.1, 1.3, 1.3, 2), 2, 2)
plus.cex = 4; beta.cex = 3; lwd = 4

PlotBetaHat <- function(){
   plot(c(0, 0), xlim=c(-1.5, 4), ylim=c(-1.5, 5), col=0, axes=FALSE, cex.lab=2, xlab='', ylab='')
        # xlab=TeX('$\\beta_1'), ylab=TeX('$\\beta_2')); 
   abline(h=0, v=0)
   points(beta[1], beta[2], pch='+', cex=plus.cex); 
   text(beta[1], beta[2], pos=4, labels=TeX('$\\widehat{\\beta}'), cex=beta.cex)
   text(3.5, 0.15, pos=4, labels=TeX('$\\beta_1'), cex=beta.cex)
   text(0, 4.5, pos=4, labels=TeX('$\\beta_2'), cex=beta.cex)
   sapply(alphaList, function(alpha){
      ellipse(beta, Sigma, alpha=alpha, col=8, lwd=lwd, lty=2, type='l', xlab='', ylab='')
   })
}

pdf('FigECAS-RegularizationL0.pdf'); par(mex=.01)
PlotBetaHat()
abline(h=0, v=0, col=4, lwd=lwd)
ellipse(beta, Sigma, alpha=.37, col=8, lwd=lwd, type='l', xlab='', ylab='')
points(0, 1.2, col=2, pch='+', cex=plus.cex); 
dev.off()

pdf('FigECAS-RegularizationL1.pdf'); par(mex=.01)
PlotBetaHat()
l1Ball = matrix(c(1, 0, -1, 0, 1, 0, 1, 0, -1, 0), 5, 2)
lines(l1Ball, col=4, lwd=lwd)
ellipse(beta, Sigma, alpha=.34, col=8, lwd=lwd, type='l', xlab='', ylab='')
points(0, 1, col=2, pch='+', cex=plus.cex); 
dev.off()

pdf('FigECAS-RegularizationL2.pdf'); par(mex=.01)
PlotBetaHat()
alpha = 2*acos(-1)*seq(0, 1, by=.01); 
l2Ball = cbind(cos(alpha), sin(alpha))
lines(l2Ball, col=4, lwd=lwd)
ellipse(beta, Sigma, alpha=.37, col=8, lwd=lwd, type='l', xlab='', ylab='')
points(0.2, .95, col=2, pch='+', cex=plus.cex); 
dev.off()
