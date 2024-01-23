# Example from Barents fish data

# Parms
rm(list=ls())
exportPlot <- TRUE; mex <- .9
par(mfrow=c(1, 1), mex=mex)
source('/home/robin/PgmR/General/FunctionsMatVec.R')
library(PLNmodels); library(EMtree)

# Dirs
dataDir <- '/home/robin/RECHERCHE/ECOLOGIE/CountPCA/pln/Data/BBVA/'
dataName <- 'BarentsFish'
figDir <- '/home/robin/RECHERCHE/ECOLOGIE/EXPOSES/FIGURES/'

# Functions
F_SignedLog <- function(x){sign(x)*log(abs(x))}
F_SignedPower <- function(x, alpha=.5){sign(x)*abs(x)^alpha}
F_PlotCrit <- function(plnFit, lwd=2){
   plot(plnFit$criteria$loglik, type='b', lwd=lwd, pch=20, xlab='', ylab='criteria', cex.lab=2, cex.axis=1.5)
   lines(plnFit$criteria$BIC, type='b', lwd=lwd, col=4, pch=20); 
   abline(v=which.max(plnFit$criteria$BIC), lty=2, lwd=lwd, col=4)
   lines(plnFit$criteria$ICL, type='b', lwd=lwd, col=2, pch=20); 
   abline(v=which.max(plnFit$criteria$ICL), lty=2, lwd=lwd, col=2)
   legend(x=.75*length(plnFit$criteria$loglik), y=mean(range(plnFit$criteria$loglik)), 
          c('ELBO', 'BIC', 'ICL'), col=c(1,4, 2), lty=1, lwd=lwd, cex=1.5)
}

# Data
load(paste0(dataDir, dataName, '.Rdata'))
counts <- Data$count; covariatesOrg <- Data$covariates; offset <- Data$offset
covariates <- scale(covariatesOrg); 
n <- nrow(counts); p <- ncol(counts); d <- ncol(covariates)
head(covariates)

###############################################################################
# PLN
plnNull <- PLN(counts ~ 1)
plnAll <- PLN(counts ~ covariates)

# Regression coefficients
if(exportPlot){pdf(paste0(figDir, dataName, '-coeffAll.pdf'))}
image(0:d, 1:p, t(plnAll$model_par$Theta), xlab='covariates', ylab='species', cex.lab=2)
if(exportPlot){dev.off()}

if(exportPlot){pdf(paste0(figDir, dataName, '-coeffAll-woIntercept.pdf'))}
image(1:d, 1:p, t(plnAll$model_par$Theta[, -1]), xlab='covariates', ylab='species', cex.lab=2)
# image(1:p, 1:d, plnAll$model_par$Theta[, -1], xlab='species', ylab='covariates', cex.lab=2)
if(exportPlot){dev.off()}

predAll <- cbind(rep(1, n), covariates)%*%t(plnAll$model_par$Theta)

# Covariance
# par(mfrow=c(2, 2), mex=.6, pch=20)
# image(1:p, 1:p, plnNull$model_par$Sigma, xlab='species', ylab='species')
# image(1:p, 1:p, cov(predAll), xlab='species', ylab='species')
# image(1:p, 1:p, plnAll$model_par$Sigma, xlab='species', ylab='species')
# image(1:p, 1:p, plnAll$model_par$Sigma+cov(predAll), xlab='species', ylab='species')

# par(mfrow=c(1, 1), mex=.6, pch=20)
# plot(F_SignedPower(plnNull$model_par$Sigma), F_SignedPower(cov(predAll)))
# points(F_SignedPower(plnNull$model_par$Sigma), F_SignedPower(plnAll$model_par$Sigma), col=4)
# points(F_SignedPower(plnNull$model_par$Sigma), F_SignedPower(cov(predAll) + plnAll$model_par$Sigma), col=2)
# abline(0, 1)

if(exportPlot){pdf(paste0(figDir, dataName, '-corrNull.pdf'))}
image(1:p, 1:p, cov2cor(plnNull$model_par$Sigma), xlab='species', ylab='species', cex.lab=2)
if(exportPlot){dev.off()}
if(exportPlot){pdf(paste0(figDir, dataName, '-corrPred.pdf'))}
image(1:p, 1:p, cor(predAll), xlab='species', ylab='species', cex.lab=2)
if(exportPlot){dev.off()}
if(exportPlot){pdf(paste0(figDir, dataName, '-corrAll.pdf'))}
image(1:p, 1:p, cov2cor(plnAll$model_par$Sigma), xlab='species', ylab='species', cex.lab=2)
if(exportPlot){dev.off()}

#################################N#############################################
# PCA
# plnPcaNull <- PLNPCA(counts ~ 1, ranks=1:12)
# plnPcaAll <- PLNPCA(counts ~ covariates, ranks=1:12)
if(exportPlot){pdf(paste0(figDir, dataName, '-pcaCritNull.pdf'))}
# plot(plnPcaNull)
F_PlotCrit(plnPcaNull)
if(exportPlot){dev.off()}
if(exportPlot){pdf(paste0(figDir, dataName, '-pcaCritAll.pdf'))}
# plot(plnPcaAll)
F_PlotCrit(plnPcaAll)
if(exportPlot){dev.off()}

plnPcaNullSel <- getModel(plnPcaNull, 2); plnPcaNullInd <- plnPcaNullSel$plot_individual_map()
PCnull <- cbind(plnPcaNullInd$data$a1, plnPcaNullInd$data$a2)
plnPcaAllSel <- getModel(plnPcaAll, 2); plnPcaAllInd <- plnPcaAllSel$plot_individual_map()
PCall <- cbind(plnPcaAllInd$data$a1, plnPcaAllInd$data$a2)
h <- 1; k <- 4; 
if(exportPlot){pdf(paste0(figDir, dataName, '-compPcaNull', h, '-', colnames(covariatesOrg)[k], '.pdf'))}
plot(PCnull[, h], covariatesOrg[, k], xlab=paste0('PCnull', h), ylab=colnames(covariatesOrg)[k], cex.lab=2, pch=20)
if(exportPlot){dev.off()}
h <- 1; k <- 1; 
if(exportPlot){pdf(paste0(figDir, dataName, '-compPcaNull', h, '-', colnames(covariatesOrg)[k], '.pdf'))}
plot(PCnull[, h], covariatesOrg[, k], xlab=paste0('PCnull', h), ylab=colnames(covariatesOrg)[k], cex.lab=2, pch=20)
if(exportPlot){dev.off()}
h <- 2; k <- 3; 
if(exportPlot){pdf(paste0(figDir, dataName, '-compPcaNull', h, '-', colnames(covariatesOrg)[k], '.pdf'))}
plot(PCnull[, h], covariatesOrg[, k], xlab=paste0('PCnull', h), ylab=colnames(covariatesOrg)[k], cex.lab=2, pch=20)
if(exportPlot){dev.off()}

if(exportPlot){pdf(paste0(figDir, dataName, '-pcaIndNull.pdf'))}
plnPcaNullSel$plot_individual_map()
if(exportPlot){dev.off()}
if(exportPlot){pdf(paste0(figDir, dataName, '-pcaIndAll.pdf'))}
plnPcaAllSel$plot_individual_map()
if(exportPlot){dev.off()}
if(exportPlot){pdf(paste0(figDir, dataName, '-pcaCorNull.pdf'))}
plnPcaNullSel$plot_correlation_circle()
if(exportPlot){dev.off()}
if(exportPlot){pdf(paste0(figDir, dataName, '-pcaCorAll.pdf'))}
plnPcaAllSel$plot_correlation_circle()
if(exportPlot){dev.off()}
# Dnull <- as.matrix(dist(PCnull[, 1:2])); Dall <- as.matrix(dist(PCall[, 1:2]))
# if(exportPlot){pdf(paste0(figDir, dataName, '-compDist-pcaNull-pcaAll.pdf'))}
# plot(Dnull, Dall, xlab='distance PCnull', ylab='distance PCall', pch=20, cex.lab=2)
# if(exportPlot){dev.off()}

###############################################################################
# Network : glasso
plnNetNull <- PLNnetwork(counts ~ 1)
plnNetNull$plot()
plnNetAll <- PLNnetwork(counts ~ covariates, control_init=list(inception=plnAll)); 
plnNetAll$plot()
lambda <- plnNetAll$penalties; lambda <- exp(seq(min(log(lambda))-max(log(lambda)), max(log(lambda)), length.out=2*length(lambda)-1))
plnNetAll <- PLNnetwork(counts ~ covariates, penalties=lambda, control_init=list(inception=plnAll)); 
if(exportPlot){pdf(paste0(figDir, dataName, '-netCrit.pdf'))}
plnNetAll$plot()
if(exportPlot){dev.off()}

# Inferred network
plnGraphNull <- plnNetNull$getBestModel()
plnGraphNull$plot_network()
plnGraphAll <- plnNetAll$getBestModel()
plnGraphAll$plot_network()

###############################################################################
# Network : tree-based
plnEMtreeNull <- EMtree(plnNull)
plnEMtreeAll <- EMtree(plnAll)
if(exportPlot){pdf(paste0(figDir, dataName, '-EMtreeNull.pdf'))}
image(1:p, 1:p, plnEMtreeNull$edges_prob, xlab='species', ylab='species')
if(exportPlot){dev.off()}
if(exportPlot){pdf(paste0(figDir, dataName, '-EMtreeAll.pdf'))}
image(1:p, 1:p, plnEMtreeAll$edges_prob, xlab='species', ylab='species')
if(exportPlot){dev.off()}

###############################################################################
# Sampling from the conditional tree distribution