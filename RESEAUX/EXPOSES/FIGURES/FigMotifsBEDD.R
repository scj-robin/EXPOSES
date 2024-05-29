# Figure pour Motif & B-EDD

rm(list=ls())
library(gtools); library(sna); library(bmotif); library(polyaAeppli)

# Dirs
dirFig <- './'
dirCodeGOF <- '/home/robin/RECHERCHE/RESEAUX/GOF-Network/GOFnetwork2018/Pgm/EDD/Functions/'
dirCodeSBM <- '/home/robin/ENSEIGN/COURS/RESEAUX/EcoStat/Materiel pour diffusion/'
source(paste0(dirCodeGOF, 'MotifCount.R'))
source(paste0(dirCodeGOF, 'MotifSimul.R'))
source(paste0(dirCodeGOF, 'MotifAnalysis.R'))
# source(paste0(dirCodeSBM, 'function_for_blockmodels.R'))

# Figs
mex <- 1; cex <- 10; lwdEdge <- 4; lwdCurve <- 8; cex.main <- 8; 
topCh <- 1; bottomCh <- 0; edgeCol <- 1
exportFig <- TRUE
figName <- paste0(dirFig, 'FigMotifsBEDD')

# Motif characteristics
sixNode <- TRUE
ifelse(sixNode, load(paste0(dirCodeGOF, 'MotifsCharacteristics-6.RData')), load(paste0(dirCodeGOF, 'MotifsCharacteristics-5.RData')))

###############################################################################
# Modele B-EDD
###############################################################################
# Dims and Parms
m <- 30; n <- 50; rho <- .1; lambda.gList <- c(1, 2.5); lambda.hList <- c(1, 4); 
seed <- 3; set.seed(seed)
F_Glambda <- function(u, lambda){lambda*u^(lambda-1)}
F_PlotNet <- function(lambda.g, lambda.h, cex.main=cex.main, lwd=lwdEdge){
   net <- SimulB_EDD(m, n, rho, lambda.g, lambda.h)$net
   if(lambda.g!=1){net <- net[order(rowSums(net)), ]}
   if(lambda.h!=1){net <- net[, order(colSums(net))]}
   PlotBNet(net, cex.main=cex.main, lwd=lwdEdge, topch=topCh, bottomch=bottomCh, edgecol=edgeCol)
}
uList <- seq(0, 1, length.out=100)
par(mfrow=c(3, 3), mex=mex); 
plot(0, 0, col=0, axes=0, xlab='', ylab='')

if(exportFig){png(paste0(figName, '-dist-g', round(10*lambda.gList[1]), '.png')); par(mfrow=c(1, 1), mex=mex)}
plot(uList, F_Glambda(uList, lambda.gList[1]), xlab='', ylab='', type='l', lwd=lwdCurve, col.axis=0, col=4)
if(exportFig){dev.off()}
if(exportFig){png(paste0(figName, '-dist-g', round(10*lambda.gList[2]), '.png')); par(mfrow=c(1, 1), mex=mex)}
plot(uList, F_Glambda(uList, lambda.gList[2]), xlab='', ylab='', type='l', lwd=lwdCurve, col.axis=0, col=4)
if(exportFig){dev.off()}

if(exportFig){png(paste0(figName, '-dist-h', round(10*lambda.hList[1]), '.png')); par(mfrow=c(1, 1), mex=mex)}
plot(uList, F_Glambda(uList, lambda.hList[1]), xlab='', ylab='', type='l', lwd=lwdCurve, col.axis=0, col=4)
if(exportFig){dev.off()}
if(exportFig){png(paste0(figName, '-adj-g', round(10*lambda.gList[1]), '-h', round(10*lambda.hList[1]), '.png')); par(mfrow=c(1, 1), mex=mex)}
net <- SimulB_EDD(m, n, rho, lambda.gList[1], lambda.hList[1]); net <- net$net[order(net$u), order(net$v)]
image(1:n, 1:m, t(net), xlab='', ylab='', axes=0)
if(exportFig){dev.off()}
if(exportFig){png(paste0(figName, '-adj-g', round(10*lambda.gList[2]), '-h', round(10*lambda.hList[1]), '.png')); par(mfrow=c(1, 1), mex=mex)}
net <- SimulB_EDD(m, n, rho, lambda.gList[2], lambda.hList[1]); net <- net$net[order(net$u), order(net$v)]
image(1:n, 1:m, t(net), xlab='', ylab='', axes=0)
if(exportFig){dev.off()}

if(exportFig){png(paste0(figName, '-dist-h', round(10*lambda.hList[2]), '.png')); par(mfrow=c(1, 1), mex=mex)}
plot(uList, F_Glambda(uList, lambda.hList[2]), xlab='', ylab='', type='l', lwd=lwdCurve, col.axis=0, col=4)
if(exportFig){dev.off()}
if(exportFig){png(paste0(figName, '-adj-g', round(10*lambda.gList[1]), '-h', round(10*lambda.hList[2]), '.png')); par(mfrow=c(1, 1), mex=mex)}
net <- SimulB_EDD(m, n, rho, lambda.gList[1], lambda.hList[2]); net <- net$net[order(net$u), order(net$v)]
image(1:n, 1:m, t(net), xlab='', ylab='', axes=0)
if(exportFig){dev.off()}
if(exportFig){png(paste0(figName, '-adj-g', round(10*lambda.gList[2]), '-h', round(10*lambda.hList[2]), '.png')); par(mfrow=c(1, 1), mex=mex)}
net <- SimulB_EDD(m, n, rho, lambda.gList[2], lambda.hList[2]); net <- net$net[order(net$u), order(net$v)]
image(1:n, 1:m, t(net), xlab='', ylab='', axes=0)
if(exportFig){dev.off()}

###############################################################################
# Proprietes des motifs : automorphismes
###############################################################################
s <- 10
motif <- motifList[[s]]
sapply(1:max(length(motif$Rset), 6), function(r){
   if(exportFig){png(paste0(figName, '-motif', s, '-automorphism', r, '.png')); par(mfrow=c(1, 1), mex=mex)}
   PlotBNet(motif$Rset[[r]], cex=cex, lwd=lwdEdge, topch=topCh, bottomch=bottomCh, edgecol=edgeCol)
   if(exportFig){dev.off()}
})

###############################################################################
# Proprietes des motifs : probabilite
###############################################################################
adj <- motif$Rset[[1]]
if(exportFig){png(paste0(figName, '-motif', s, '.png')); par(mfrow=c(1, 1), mex=mex)}
PlotBNet(adj, cex=cex, lwd=lwdEdge, topch=topCh, bottomch=bottomCh, edgecol=edgeCol)
if(exportFig){dev.off()}
sapply(unique(motif$topdegree), function(u){
   if(exportFig){png(paste0(figName, '-motif', s, '-top', u, '.png')); par(mfrow=c(1, 1), mex=mex)}
   PlotBNet(motifList[[starNum$top[u]]]$A, cex=cex, lwd=lwdEdge, topch=topCh, bottomch=bottomCh, edgecol=edgeCol)
   if(exportFig){dev.off()}
})
sapply(unique(motif$bottomdegree), function(v){
   if(exportFig){png(paste0(figName, '-motif', s, '-bottom', v, '.png')); par(mfrow=c(1, 1), mex=mex)}
   PlotBNet(motifList[[starNum$bottom[v]]]$A, cex=cex, lwd=lwdEdge, topch=topCh, bottomch=bottomCh, edgecol=edgeCol)
   if(exportFig){dev.off()}
})

###############################################################################
# Proprietes des motifs : super-motifs
###############################################################################
superMotif <- superMotifList[[s]][[s]]; 
sapply(1:16, function(r){
   if(exportFig){png(paste0(figName, '-motif', s, '-supermotif', r, '.png')); par(mfrow=c(1, 1), mex=mex)}
   PlotBNet(superMotif[[r]], cex=cex, lwd=lwdEdge, topch=topCh, bottomch=bottomCh, edgecol=edgeCol)
   if(exportFig){dev.off()}
})

###############################################################################
# Normalité asymptotique
###############################################################################
dirPgm <- '/home/robin/RECHERCHE/RESEAUX/GOF-Network/GOFnetwork2018/Pgm/EDD/'
source(paste0(dirPgm, 'Functions/MotifCount.R'))
source(paste0(dirPgm, 'Functions/MotifAnalysis.R'))
source(paste0(dirPgm, 'Functions/MotifSimul.R'))
dirSim <- '/home/robin/RECHERCHE/RESEAUX/GOF-Network/GOFnetwork2018/Pgm/EDD/Simul/'
mList <- c(50, 100, 200)
rho <- 0.05; lambda.g <- 2; lambda.h <- 3; B <- 5e2; 
if(sixNode){
   load(paste0(dirPgm, 'Functions/MotifsCharacteristics-6.RData'))
}else{
   load(paste0(dirPgm, 'Functions/MotifsCharacteristics-5.RData'))
}
x <- seq(-5, 5, length.out=200); phi <- dnorm(x)

# Plots: motifs counts
par(mfrow=c(4, 4), mex=.6, lwd=2);
for(m in mList){
   # Simulated data
   n <- round(2*m/3); 
   simName <- paste0('SimMotifBEDD-m', m, '-n', n, '-rho', round(100*rho), '-g', lambda.g, '-h', lambda.h, '-6nodes', sixNode, '-B', B)
   simFile <- paste0(dirSim, simName, '.RData')
   load(simFile)
   for (s in c(2, nonStarNb)){
      # Motifs
      if(exportFig){png(paste0(figName, '-adjMatMotif', nonStarNum[s], '.png')); par(mfrow=c(1, 1), mex=mex, lwd=2)}
      PlotBNet(motifList[[nonStarNum[s]]]$A, cex=cex, lwd=lwdEdge, topch=topCh, bottomch=bottomCh, edgecol=edgeCol)
      if(exportFig){dev.off()}
      # Count distribution
      if(exportFig){png(paste0(figName, '-', simName, '-distCountMotif', nonStarNum[s], '.png')); par(mfrow=c(1, 1), mex=mex, lwd=2)}
      H <- hist(N[, s], breaks=sqrt(B), main='', xlab='', ylab='', freq=FALSE)
      theoMeanN <- theoMoments$mean[nonStarNum]; theoSdN <- theoMoments$sd[nonStarNum]
      lines(sort(N[, s]), dnorm(sort(N[, s]), mean=theoMeanN[s],
                                sd=theoSdN[s]), col=4, lwd=lwdCurve)
      if(theoMeanN[s] <= 1e4){
         compPoisTheo <- CompPoissonParms(theoMeanN[s], theoSdN[s]^2)
         lines(sort(N[, s]), dPolyaAeppli(sort(N[, s]),
                                          compPoisTheo$lambda, compPoisTheo$prob), col=6, lwd=lwdCurve)
      }
      if(exportFig){dev.off()}
      # Raw statistics
      if(exportFig){png(paste0(figName, '-', simName, '-distRawStat', nonStarNum[s], '.png')); par(mfrow=c(1, 1), mex=mex, lwd=2)}
      H <- hist(statN_Nbar[, s], breaks=sqrt(B), main='', xlab='', ylab='', freq=FALSE)
      lines(x, phi, lwd=lwdCurve, col=4)
      if(exportFig){dev.off()}
      # Normalized statistics
      if(exportFig){png(paste0(figName, '-', simName, '-distNormStat', nonStarNum[s], '.png')); par(mfrow=c(1, 1), mex=mex, lwd=2)}
      H <- hist(normStatN_Nbar[, s], breaks=sqrt(B), main='', xlab='', ylab='', freq=FALSE)
      lines(x, phi, lwd=lwdCurve, col=2)
      if(exportFig){dev.off()}
      # Raw statistics
      if(exportFig){png(paste0(figName, '-', simName, '-distRawStat', nonStarNum[s], '.png')); par(mfrow=c(1, 1), mex=mex, lwd=2)}
      H <- hist(statN_Nbar[, s], breaks=sqrt(B), main='', xlab='', ylab='', freq=FALSE)
      lines(x, phi, lwd=lwdCurve, col=4)
      if(exportFig){dev.off()}
   }
}
