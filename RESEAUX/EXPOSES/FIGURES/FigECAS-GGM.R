# Figures for ECAS / SFdS / NAIM : reminder on prob. dist

library(sna); library(mvtnorm); library(huge); library(LITree); library(ROCR); 
library(mixtools); library(latex2exp)
source('FigECAS-Functions.R')
source('/home/robin/PgmR/General/FunctionsMatVec.R')
source('/home/robin/RECHERCHE/RESEAUX/S-Founas/NoisyNetworkInference/Pgm/Functions/SimulFunctions.R')
source('/home/robin/RECHERCHE/RESEAUX/S-Founas/NoisyNetworkInference/Pgm/Functions/ResultFunctions.R')

###############################################################################
# Small example
seed = 3; set.seed(seed)
p = 4; G = matrix(c(0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0), 4, 4)
pdf('FigECAS-GGM-4nodes.pdf')
par(mex=0.01)
gplot(G, gmode='graph', label=1:p, label.pos=5, vertex.col=0, vertex.cex=2, label.cex=2)
dev.off()
OmegaSigma = MakeOmegaSigma(G); Omega = OmegaSigma$Omega; Sigma = OmegaSigma$Sigma
corSigma = cov2cor(Sigma); corOmega = -cov2cor(Omega); diag(corOmega) = 1

Ftab(cbind(round(Sigma, 1))); 
Ftab(round(corSigma, 2))
Ftab(cbind(round(Omega, 1))); 
Ftab(round(corOmega, 2))

# Simulated GGM 4 nodes
n = 25; Y = rmvnorm(n, sigma=Sigma)
SigmaHat = cov(Y); corSigmaHat = cor(Y)
OmegaHat = solve(SigmaHat); corOmegaHat = -cov2cor(OmegaHat); diag(corOmegaHat) = 1
Ftab(cbind(round(SigmaHat, 2))); 
# Ftab(round(corSigmaHat, 2))
Ftab(cbind(round(OmegaHat, 2))); 
# Ftab(round(corOmegaHat, 2))

###############################################################################
# Missing node
seed = 1; set.seed(seed); p = 15; density = log(p)/p
# Graph generation
G = F_Vec2Sym(F_Sym2Vec(matrix(rbinom(p^2, 1, density), p, p)))
while(min(colSums(G))==0){G = F_Vec2Sym(F_Sym2Vec(matrix(rbinom(p^2, 1, density), p, p)))}
gplot(G, gmode='graph')
missNode = which.max(colSums(G)); 
# Making the missing more central
missEdge = rbinom(p, 1, 2*density); G[missNode, ] = G[, missNode] = G[missNode, ]+missEdge
G[which(G > 1)] = 1; G = F_Vec2Sym(F_Sym2Vec(G))
missNeighbors = sum(G[missNode, ])
nodeOrder = unique(c(missNode, which(G[missNode, ]==1), which(G[missNode, ]==0)))
G = G[nodeOrder, nodeOrder]; missNode = 1
gplot(G, gmode='graph')
# Create Omega & Sigma
OmegaSigma = MakeOmegaSigma(G); OmegaAll = OmegaSigma$Omega ; SigmaAll = OmegaSigma$Sigma
OmegaAll = OmegaAll; SigmaAll = SigmaAll 
corSigmaAll = cov2cor(SigmaAll); corOmegaAll = -cov2cor(OmegaAll); diag(corOmegaAll) = 1
# Create Omega & Sigma with missing
OmegaMiss = OmegaAll[-missNode, -missNode] - 
   OmegaAll[-missNode, missNode]%o%OmegaAll[missNode, -missNode]/(OmegaAll[missNode, missNode])
Gmiss = 1*(OmegaMiss!=0)
SigmaMiss = solve(OmegaMiss); corSigmaMiss = cov2cor(SigmaMiss); 
corOmegaMiss = -cov2cor(OmegaMiss); diag(corOmegaMiss) = 1

pdf('FigECAS-MissingNode.pdf')
par(mfrow=c(2, 3), mex=.6)
alpha = 1/3
nodeCol = rep(0, p); nodeCol[missNode] = 2
nodePos = gplot(G, gmode='graph', label=(1:p), label.pos=5, vertex.col=nodeCol, vertex.cex=2.5, main='complete GM')
image(1:p, 1:p, (abs(corSigmaAll)^alpha), zlim=c(0, 1), main='complete Sigma'); abline(v=.5+c(1, 1+missNeighbors), h=.5+c(1, 1+missNeighbors))
image(1:p, 1:p, (abs(corOmegaAll)^alpha), zlim=c(0, 1), main='complete Omega'); abline(v=.5+c(1, 1+missNeighbors), h=.5+c(1, 1+missNeighbors))
gplot(Gmiss, gmode='graph', coord=nodePos[-missNode, ], label=(2:p), label.pos=5, vertex.col=0, vertex.cex=2.5, main='marginal GM')
image(2:p, 2:p, (abs(corSigmaMiss)^alpha), zlim=c(0, 1), main='marginal Sigma'); abline(v=.5+1+missNeighbors, h=.5+1+missNeighbors)
image(2:p, 2:p, (abs(corOmegaMiss)^alpha), zlim=c(0, 1), main='marginal Omega'); abline(v=.5+1+missNeighbors, h=.5+1+missNeighbors)
dev.off()