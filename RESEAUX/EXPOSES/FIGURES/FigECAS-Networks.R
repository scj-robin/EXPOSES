# Figures for ECAS / SFdS / NAIM : networks

rm(list=ls())
source('FigECAS-Functions.R')
source('/home/robin/PgmR/General/FunctionsMatVec.R')

# Karate
library(blockmodels); library(igraph); library(igraphdata); library(sna)
seed = 1; set.seed(seed)
data("karate"); G = as.matrix(get.adjacency(karate)); p = nrow(G)
pdf('FigECAS-KarateUnordered.pdf')
coord = matrix(runif(2*p), p, 2); gplot(G, gmode='graph', coord=coord, vertex.col=1, edge.col=8)
dev.off()
Zhat = c(rep(1, 3), rep(2, 11), rep(3, 18), rep(4, 2))
BM = BM_bernoulli('SBM_sym', G); BM$estimate(); Zhat = apply(BM$memberships[[4]]$Z, 1, which.max)
pdf('FigECAS-KarateOrdered.pdf')
gplot(G, gmode='graph', vertex.col=Zhat, edge.col=8)
dev.off()

# Species interactions
library(vegan); library(ade4)
data(aravo); varSel = c(1, 2, 4); p = nrow(aravo$traits)
Fvec(names(aravo$traits)[varSel]) 
Ftab(cbind(1:5, as.matrix(aravo$traits[1:5, varSel])))
Dspe = as.matrix(dist(scale(t(aravo$spe)), diag=T))
# Dtraits = F_Sym2Vec(as.matrix(dist(scale(aravo$traits), diag=T)))
Inter = 1*(Dspe < 10)
invisible(sapply(1:7, function(r){
   i = sample(1:p, 1); j = sample(1:p, 1); cat(i, '&', j, '&', Inter[i, j], '\\\\ \n')
   }))
