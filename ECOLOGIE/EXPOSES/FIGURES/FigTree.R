# Example from tree network

library(sna); library(blockmodels)

exportPlot = TRUE; mex = .9
par(mfrow=c(1, 1), mex=mex)

# Data
dataDir = '/home/robin/RECHERCHE/BAYES/VBEM-IS/VBEM-IS.git/Data/Tree/'
dataName = 'Tree'
figDir = '/home/robin/RECHERCHE/RESEAUX/EXPOSES/FIGURES/'
load(paste0(dataDir, dataName, '.Rdata'))
p = ncol(interaction)

# Network
if(exportPlot){pdf(paste0(figDir, dataName, '-adjMat.pdf')); par(mex=mex)}
image(1:p, 1:p, log(1+interaction), xlab='species', ylab='species', cex.lab=2)
if(exportPlot){dev.off()}

# Blockmodel null
# BMnull = BM_poisson('SBM_sym', adj=interaction, plotting='')
# BMnull$estimate()
# save(BMnull, file=paste0(dataDir, dataName, '-SBMnull.Rdata'))
load(paste0(dataDir, dataName, '-SBMnull.Rdata'))
Q = which.max(BMnull$ICL); 
if(exportPlot){pdf(paste0(figDir, dataName, '-ICL-SBMnull.pdf')); par(mex=mex)}
plot(BMnull$ICL, pch=20, type='b', ylab='', xlab='K', lwd=3, cex.axis=2, cex.lab=2); 
abline(v = c(Q, 7), lwd=3, lty=c(1, 2))
if(exportPlot){dev.off()}
Q = 7
Tau = BMnull$memberships[[Q]]$Z; Tau = Tau[, order(colMeans(Tau))]
R = Tau%*%(1:Q); Z = apply(Tau, 1, which.max); Znull = Z
if(exportPlot){pdf(paste0(figDir, dataName, '-adjMat-SBMnull.pdf')); par(mex=mex)}
image(1:p, 1:p, log(1+interaction)[order(R), order(R)], , xlab='species', ylab='species', cex.lab=2)
abline(v=.5+cumsum(table(Z)), h=.5+cumsum(table(Z)), col=4, lwd=3)
if(exportPlot){dev.off()}

# Blockmodel taxo
# BMtaxo = BM_poisson_covariates('SBM_sym', adj=interaction, covariates=d.taxo, plotting='')
# BMtaxo$estimate()
# save(BMtaxo, file=paste0(dataDir, dataName, '-SBMtaxo.Rdata'))
load(paste0(dataDir, dataName, '-SBMtaxo.Rdata'))
Q = which.max(BMtaxo$ICL); 
if(exportPlot){pdf(paste0(figDir, dataName, '-ICL-SBMtaxo.pdf')); par(mex=mex)}
plot(BMtaxo$ICL, pch=20, type='b', ylab='', , xlab='K', lwd=3, cex.axis=2, cex.lab=2); 
abline(v = c(Q, 7), lwd=3, lty=c(1, 2))
if(exportPlot){dev.off()}
Tau = BMtaxo$memberships[[Q]]$Z; Tau = Tau[, order(colMeans(Tau))]
R = Tau%*%(1:Q); Z = apply(Tau, 1, which.max); Ztaxo = Z
if(exportPlot){pdf(paste0(figDir, dataName, '-adjMat-SBMtaxo.pdf')); par(mex=mex)}
image(1:p, 1:p, log(1+interaction)[order(R), order(R)], , xlab='species', ylab='species', cex.lab=2)
abline(v=.5+cumsum(table(Z)), h=.5+cumsum(table(Z)), col=4, lwd=3)
if(exportPlot){dev.off()}
print(BMtaxo$model_parameters[[Q]]$beta)

# Blockmodel taxo-geo
Xtaxogeo = list(d.taxo=d.taxo, d.geo=d.geo)
# BMtaxogeo = BM_poisson_covariates('SBM_sym', adj=interaction, covariates=Xtaxogeo, plotting='')
# BMtaxogeo$estimate()
# save(BMtaxogeo, file=paste0(dataDir, dataName, '-SBMtaxogeo.Rdata'))
load(paste0(dataDir, dataName, '-SBMtaxogeo.Rdata'))
Q = which.max(BMtaxogeo$ICL); 
# plot(BMtaxogeo$ICL, pch=20, type='b'); abline(v = Q, col=4, lwd=2)
Tau = BMtaxogeo$memberships[[Q]]$Z; Tau = Tau[, order(colMeans(Tau))]
R = Tau%*%(1:Q); Z = apply(Tau, 1, which.max); Ztaxogeo = Z
if(exportPlot){pdf(paste0(figDir, dataName, '-adjMat-SBMtaxogeo.pdf')); par(mex=mex)}
image(1:p, 1:p, log(1+interaction)[order(R), order(R)], , xlab='species', ylab='species', cex.lab=2)
abline(v=.5+cumsum(table(Z)), h=.5+cumsum(table(Z)), col=4, lwd=3)
if(exportPlot){dev.off()}
print(as.vector(BMtaxogeo$model_parameters[[Q]]$beta))

# Blockmodel full
Xfull = list(d.taxo=d.taxo, d.geo=d.geo, d.gene=log.d.gene)
# BMfull = BM_poisson_covariates('SBM_sym', adj=interaction, covariates=Xfull, plotting='')
# BMfull$estimate()
# save(BMfull, file=paste0(dataDir, dataName, '-SBMfull.Rdata'))
load(paste0(dataDir, dataName, '-SBMfull.Rdata'))
Q = which.max(BMfull$ICL); 
# plot(BMfull$ICL, pch=20, type='b'); abline(v = Q, col=4, lwd=2)
Tau = BMfull$memberships[[Q]]$Z; Tau = Tau[, order(colMeans(Tau))]
R = Tau%*%(1:Q); Z = apply(Tau, 1, which.max); Zfull = Z
if(exportPlot){pdf(paste0(figDir, dataName, '-adjMat-SBMfull.pdf')); par(mex=mex)}
image(1:p, 1:p, log(1+interaction)[order(R), order(R)], , xlab='species', ylab='species', cex.lab=2)
abline(v=.5+cumsum(table(Z)), h=.5+cumsum(table(Z)), col=4, lwd=3)
if(exportPlot){dev.off()}
print(as.vector(BMfull$model_parameters[[Q]]$beta))


table(Znull, Ztaxo)
table(Ztaxo, Ztaxogeo)
table(Ztaxogeo, Zfull)

if(exportPlot){pdf(paste0(figDir, dataName, '-ICL-SBMall.pdf')); par(mex=mex)}
plot(BMnull$ICL, pch=20, type='b', ylab='', xlab='K', lwd=3, cex.axis=2, cex.lab=2, ylim=c(min(BMnull$ICL), max(BMtaxo$ICL))); abline(v=which.max(BMnull$ICL), lwd=3)
points(BMtaxo$ICL, pch=20, type='b', col=2, lwd=3); abline(v=which.max(BMtaxo$ICL), col=2, lwd=3)
points(BMtaxogeo$ICL, pch=20, type='b', col=3, lwd=3); abline(v=which.max(BMtaxogeo$ICL), col=3, lwd=3)
points(BMfull$ICL, pch=20, type='b', col=4, lwd=3); abline(v=which.max(BMfull$ICL), col=4, lwd=3)
if(exportPlot){dev.off()}
