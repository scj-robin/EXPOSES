# Pictures for VBEM-IS / PLN sampling talk at UVSQ (march 2025)

rm(list=ls())
palette('R3')
library(sna)
library(weights)
library(PLNmodels)
library(blockmodels)
library(mclust)
par(cex=5)
cexLab <- 1.75
cexTick <- 1.25
exportFig <- FALSE

######################################################################################
# Barents data
ordering <- TRUE
# dev.off()
DataDir <- '/home/robin/RECHERCHE/ECOLOGIE/PLN/DataSR/'
load(paste(DataDir, 'BarentsFish.Rdata', sep=''))
X <- as.matrix(scale(Data$covariates)); Y <- as.matrix(Data$count)
n <- nrow(Y); p <- ncol(Y); d <- ncol(X)
pln <- PLN(Y ~ X)
Ypred <- exp(cbind(rep(1, n), X)%*%pln$model_par$B + rep(1, n)%o%diag(pln$model_par$Sigma)/2)
# plot(1+Y, 1+Ypred, log='xy', pch=20); abline(0, 1)
# Ordering species
gmm <- Mclust(t(pln$model_par$B[-1, ]), modelNames=c('VEI'))
if(ordering){speciesOrder <- order(gmm$z%*%(1:gmm$G))}else{speciesOrder <- 1:p}
if(exportFig){png(paste0('FigUVSQ-BarentsFish-coeffAll-woIntercept-specOrder', ordering, '.png'))}
image(1:d, 1:p, pln$model_par$B[-1, speciesOrder], xlab='covariates', ylab='species', cex.lab=cexLab, col.axis=0)
axis(1, at=1:d, labels=colnames(X), cex.axis=cexTick)
axis(2, at=5*(1:(p/5)), labels=5*(1:(p/5)), cex.axis=cexTick)
if(exportFig){dev.off()}
if(exportFig){png(paste0('FigUVSQ-BarentsFish-corrPred-specOrder', ordering, '.png'))}
image(1:p, 1:p, cor(Ypred)[speciesOrder, speciesOrder], xlab='species', ylab='species', cex.lab=cexLab)
if(exportFig){dev.off()}
if(exportFig){png(paste0('FigUVSQ-BarentsFish-corrAll-specOrder', ordering, '.png'))}
image(1:p, 1:p, cov2cor(pln$model_par$Sigma)[speciesOrder, speciesOrder], xlab='species', ylab='species', cex.lab=cexLab)
if(exportFig){dev.off()}

hist(pln$model_par$B[-1, ], breaks=sqrt(p*d))

######################################################################################
# Tree data
source('/home/robin/PgmR/General/FunctionsMatVec.R')
DataDir <- '/home/robin/RECHERCHE/BAYES/VBEM-IS/VBEM-IS.git/Pgm/TEST_TREE_DATA/'
load(paste(DataDir, 'Tree-all-data.Rdata', sep=''))
n = nrow(data$Y.mat)
# coord <- cbind(cos(2*acos(-1)*order(rowSums(data$Y.mat))/n), 
               # sin(2*acos(-1)*order(rowSums(data$Y.mat))/n))

# Network
if(exportFig){png('FigUVSQ-Tree-Network.png')}
gplot(1*(data$Y.mat>0), gmode='graph')
if(exportFig){dev.off()}
if(exportFig){png('FigUVSQ-Tree-WeightedNetwork.png')}
gplot(1*(data$Y.mat), edge.lwd=(data$Y.mat), edge.lty=1, edge.col=ceiling(10*data$Y.mat/max(data$Y.mat)), gmode='graph')
if(exportFig){dev.off()}
if(exportFig){png('FigUVSQ-Tree-Adjacency.png')}
image((1:n), (1:n), log10(1+data$Y.mat), xlab='species', ylab='species', cex.lab=cexLab)
if(exportFig){dev.off()}

# Covariates
if(exportFig){png('FigUVSQ-Tree-GeneticDistance.png')}
image((1:n), (1:n), data$Xo.list[[1]], xlab='species', ylab='species', cex.lab=cexLab)
if(exportFig){dev.off()}
if(exportFig){png('FigUVSQ-Tree-GeographicDistance.png')}
image((1:n), (1:n), data$Xo.list[[2]], xlab='species', ylab='species', cex.lab=cexLab)
if(exportFig){dev.off()}
if(exportFig){png('FigUVSQ-Tree-TaxonomicDistance.png')}
image((1:n), (1:n), data$Xo.list[[3]], xlab='species', ylab='species', cex.lab=cexLab)
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
  image((1:n), (1:n), log10(1+data$Y.mat[order, order]), xlab='species', ylab='species', cex.lab=cexLab)
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

