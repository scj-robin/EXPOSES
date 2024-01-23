# Figure for noisy network inference

rm(list=ls())
source('/home/robin/PgmR/General/FunctionsMatVec.R')
source('FigECAS-Functions.R')
source('/home/robin/RECHERCHE/RESEAUX/S-Founas/NoisyNetworkInference/Pgm/Functions/SimulFunctions.R')
library(mvtnorm); library(EMtree); library(huge); library(sna); library(PLNmodels)

# Parms
seed = 3; set.seed(seed)
 
###############################################################################
# Scores boxplots
p = 30; n = 89; density = log(p)/p; Q = 3
pi = (1:Q); pi = pi / sum(pi)
gamma = (Q:1)^2%o%(Q:1)^2; gamma = density*gamma/(pi%*%gamma%*%pi)[1, 1]
Z = t(rmultinom(p, 1, pi))
G = F_Vec2Sym(F_Sym2Vec(matrix(rbinom(p^2, 1, Z%*%gamma%*%t(Z)), p, p)))
while(min(colSums(G))==0){
   Z = t(rmultinom(p, 1, pi))
   G = F_Vec2Sym(F_Sym2Vec(matrix(rbinom(p^2, 1, Z%*%gamma%*%t(Z)), p, p)))
}
gplot(G, gmode='graph', vertex.col=1+Z%*%(1:Q));
Y = rmvnorm(n, sigma=MakeOmegaSigma(G)$Sigma)
scoreGlasso = F_Sym2Vec(fitHuge(Y, method='glasso'))
scoreMB = F_Sym2Vec(fitHuge(Y, method='mb'))
scoreTree = F_Sym2Vec(EdgeProba(-n/2*log(1- cor(Y)^2)))
par(mfrow=c(1, 3))
Gvec = F_Sym2Vec(G)
pdf('FigNoisyNet-BoxplotMB.pdf'); par(mex=.5, pch=20)
boxplot(scoreMB ~ Gvec, xlab='', ylab='', main='', cex.axis=1.5, cex.main=2, lwd=2)
dev.off()
pdf('FigNoisyNet-BoxplotGlasso.pdf'); par(mex=.5, pch=20)
boxplot(scoreGlasso ~ Gvec, xlab='', ylab='', main='', cex.axis=1.5, cex.main=2, lwd=2)
dev.off()
pdf('FigNoisyNet-BoxplotTree.pdf'); par(mex=.5, pch=20)
boxplot(scoreTree ~ Gvec, xlab='', ylab='', main='', cex.axis=1.5, cex.main=2, lwd=2)
dev.off()

###############################################################################
# Schematic model
# Dim
seed = 2; set.seed(seed); par(mfrow=c(3, 3), mex=.5)
p = 30; n = 89; Q = 3; density = 4*log(p)/p; speciesName = paste0('sp', (1:p))
# Parms
pi = (1:Q); pi = pi / sum(Q)
gamma = (Q:1)^1 %o% (Q:1)^1; gamma = density * gamma / (t(pi)%*%gamma%*%pi)[1, 1]
Z = t(rmultinom(p, 1, pi))
G = F_Vec2Sym(F_Sym2Vec(matrix(rbinom(p^2, 1, Z%*%gamma%*%t(Z)), p, p)))
while(min(colSums(G))==0){
   Z = t(rmultinom(p, 1, pi))
   F_Vec2Sym(F_Sym2Vec(matrix(rbinom(p^2, 1, Z%*%gamma%*%t(Z)), p, p)))
   }
Sigma = cov2cor(MakeOmegaSigma(G)$Sigma)
# Data
Ylatent = rmvnorm(n, mean=1+rnorm(p), sigma=Sigma)
Y = matrix(rpois(n*p, exp(Ylatent)), n, p); colnames(Y) = speciesName
# Edge scores
pln = PLN(Y ~ 1)
weightTree = -n/2*log(1 - cov2cor(pln$model_par$Sigma)^2)
S = EdgeProba(weightTree); colnames(S) = rownames(S) = speciesName
# Plots
pdf('FigNoisyNet-Graph.pdf'); par(mex=.01)
nodePos = gplot(G, gmode='graph', vertex.col=0, vertex.cex=1.5, 
      label=speciesName, label.cex=2, label.pos=5)
dev.off()
pdf('FigNoisyNet-NodesGroup.pdf'); par(mex=.01)
gplot(matrix(0, p, p), gmode='graph', vertex.col=1+(Z%*%(1:Q)), vertex.cex=1.5, coord=nodePos, 
      label=speciesName, label.cex=2, label.pos=5)
dev.off()
# Tables
specSel = 1:5
Fvec(speciesName[specSel])
Ftab(Y[1:6, specSel])
Ftab(round(100*S[specSel, specSel], 1))
