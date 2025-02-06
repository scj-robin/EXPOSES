# Pictures for VBEM-IS talk

rm(list=ls())
palette('R3')
library(sna)
library(weights)
library(blockmodels)
par(cex=5)
exportFig <- TRUE

######################################################################################
# Tree data
source('/home/robin/PgmR/General/FunctionsMatVec.R')
DataDir <- '/home/robin/RECHERCHE/BAYES/VBEM-IS/VBEM-IS.git/Pgm/TEST_TREE_DATA/'
load(paste(DataDir, 'Tree-all-data.Rdata', sep=''))
n = nrow(data$Y.mat)

# Network
if(exportFig){png('FigUVSQ-Tree-Network.png')}
gplot(1*(data$Y.mat>0), gmode='graph')
if(exportFig){dev.off()}
if(exportFig){png('FigUVSQ-Tree-Adjacency.png')}
image((1:n), (1:n), log10(1+data$Y.mat), xlab='species', ylab='species')
if(exportFig){dev.off()}

# Covariates
if(exportFig){png('FigUVSQ-Tree-GeneticDistance.png')}
image((1:n), (1:n), data$Xo.list[[1]], xlab='species', ylab='species')
if(exportFig){dev.off()}
if(exportFig){png('FigUVSQ-Tree-GeographicDistance.png')}
image((1:n), (1:n), data$Xo.list[[2]], xlab='species', ylab='species')
if(exportFig){dev.off()}
if(exportFig){png('FigUVSQ-Tree-TaxonomicDistance.png')}
image((1:n), (1:n), data$Xo.list[[3]], xlab='species', ylab='species')
if(exportFig){dev.off()}

######################################################################################
# VEM SBM
resFile <- 'FigUVSQ-TreeSBM.Rdata'
if(!file.exists(resFile)){
  sbm0 <- BM_poisson(membership_type='SBM_sym', adj=data$Y.mat)
  sbm0$estimate()
  sbm1 <- BM_poisson_covariates(membership_type='SBM_sym', adj=data$Y.mat, covariates=data$Xo.list)
  sbm1$estimate()
  save(sbm0, sbm1, file=resFile)
}else{load(resFile)}

######################################################################################
# VEM Plots
PlotAdjClust <- function(sbm){
  Q <- which.max(sbm$ICL)
  order <- order(sbm$memberships[[Q]]$Z%*%(1:Q))
  image((1:n), (1:n), log10(1+data$Y.mat[order, order]), xlab='species', ylab='species')
  abline(v=0.5+cumsum(colMeans(n*sbm$memberships[[Q]]$Z)), 
         h=0.5+cumsum(colMeans(n*sbm$memberships[[Q]]$Z)), col=4, lwd=2)
}

if(exportFig){png('FigUVSQ-Tree-ClustAdjacency-noCovar.png')}
PlotAdjClust(sbm0)
if(exportFig){dev.off()}
if(exportFig){png('FigUVSQ-Tree-ClustAdjacency-Covar.png')}
PlotAdjClust(sbm1)
if(exportFig){dev.off()}

round(sbm1$model_parameters[[which.max(sbm1$ICL)]]$beta, 3)
