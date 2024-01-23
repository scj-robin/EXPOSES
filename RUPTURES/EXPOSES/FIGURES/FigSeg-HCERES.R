# Figures segmentation -> Cambridge
rm(list=ls())

seed <- 5; set.seed(seed)
exportFig = TRUE
figName <- 'FigSeg-HCERES'

################################################################################
# Discrete time
################################################################################
# Parms
n <- 150; K <- 5
tau <- c(0, sort(sample(2:(n-1), (K-1))), n)
mu <- 3*rnorm(K)
muVec <- unlist(sapply(1:K, function(k){rep(mu[k], tau[k+1]-tau[k])}))
Y <- rnorm(n) + muVec
 
# Plot
if(exportFig){png(paste0(figName, '-discrete-seg.png'))}
plot(Y, xlab='t', ylab='Yt', pch=20, cex=1.5)
sapply(1:K, function(k){lines(.5+c(tau[k], tau[k+1]), 
                              rep(mean(Y[(tau[k]+1):tau[k+1]]), 2), col=2, lwd=2)})
abline(v=tau+.5, lwd=2, lty=2, col=2)
if(exportFig){dev.off()}

################################################################################
# Continuous time: Vulcano
################################################################################
# Data
library(CptPointProcess)
dataDir <- '/home/robin/RECHERCHE/RUPTURES/Hawkes/segHP/data/vulcano/'
dataName <- 'Kilauea';  Kmax <- 4; # Kopt <- 4
# dataName <- 'MaunaLoa'; Kmax <- 10; Kopt <- 6
times <- as.vector(read.table(paste0(dataDir, dataName, '-PP.csv'), sep=';')$x)
selList <- list(); selList[[1]] <- c(1, 3); selList[[2]] <- c(2, 4); 

dataDir <- '/home/robin/RECHERCHE/RUPTURES/Hawkes/DataSR/chiroptères/'
dataName <- 'chiro'; load(paste0(dataDir, 'sequences.Rdata'))
seqNum <- 2295; seq <- seqList[[seqNum]]; Kmax <- 5
times <- seq$times
selList <- list(); selList[[1]] <- c(1, 3); selList[[2]] <- c(2, 5);

pp <- as.data.frame(times)
cpt <- CptPointProcess(pp, Kmax=Kmax)
if(exportFig){png(paste0(figName, '-', dataName, '-seg.png'))}
plotCptPP(pp, cpt$SegK)
if(exportFig){dev.off()}
cpt$K.est
cpt$SegK
for(sel in selList){
  cpt$SegK$lambda[sel] <- cpt$SegK$DN[sel]%*%cpt$SegK$lambda[sel]/sum(cpt$SegK$DN[sel])
}
cpt$SegK

if(exportFig){png(paste0(figName, '-', dataName, '-clustSeg.png'))}
plotCptPP(pp, cpt$SegK)
if(exportFig){dev.off()}
