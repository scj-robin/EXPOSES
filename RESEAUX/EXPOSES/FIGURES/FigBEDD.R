# Figure pour Motif & B-EDD

rm(list=ls())
source('/home/robin/Bureau/GOF-Network/GOFnetwork2018/Pgm/EDD/Functions/MotifAnalysis.R')
library(lattice); library(gtools)
figName <- 'FigBEDD'
exportFig <- TRUE

################################################################################
# Examples
################################################################################
# Dirs
dirFig <- './'
## Bipartite
# dirData <- '/home/robin/Bureau/RECHERCHE/RESEAUX/BipartiteNetwork/Data/ClassifyingBipartiteNetworks/'
# load(paste0(dirData, 'ecologicalinteractions-binaryNetworks.Rdata'))
# length(networks); netNames <- unlist(lapply(networks, function(net){net$name}))
# netNums <- c(544, 554)
# netList <- list(Robertson=networks[[netNums[1]]], Silva=networks[[netNums[2]]])
## Web of life
# dirData <- '/home/robin/Bureau/RECHERCHE/RESEAUX/BipartiteNetwork/Data/WebOfLife/'
# load(paste0(dirData, 'web-of-life-binaryNetworks.Rdata'))
# length(networks); netNames <- unlist(lapply(networks, function(net){net$name}))
# netNums <- 12; net <- networks[[netNums]]
## Zackenberg
dataName <- 'Zackenberg'
dirData <- '/home/robin/Bureau/GOF-Network/GOFnetwork2018/Data/Zackenberg/'
load(paste0(dirData, 'Zackenberg.RData'))
netNums <- 3; net <- netList[[netNums]]
dataName <- paste0(dataName, '-', net$name)

# Reduction
net$adjMat <- net$adjMat[-which(rowSums(net$adjMat)==0), -which(colSums(net$adjMat)==0)]
dim(net$adjMat)

################################################################################
# Excheangeability
################################################################################
# Permutation
seed <- 1; set.seed(seed)
netPerm <- list()
netPerm[[1]] <- net$adjMat
netPerm[[2]] <- net$adjMat[order(runif(nrow(net$adjMat))), order(runif(ncol(net$adjMat)))]
netPerm[[3]] <- net$adjMat[order(rowSums(net$adjMat)), order(colSums(net$adjMat))]

# Plots
for(perm in 1:3){
  if(exportFig){png(paste0(dirFig, figName, '-', dataName, '-Net', perm, '.png'))}
  PlotBNet(netPerm[[perm]], cex=2, topch=1, bottomch=0)
  if(exportFig){dev.off()}
  if(exportFig){png(paste0(dirFig, figName, '-', dataName, '-Adj', perm, '.png'))}
  image(1:ncol(netPerm[[perm]]), 1:nrow(netPerm[[perm]]), t(netPerm[[perm]]), xlab='', ylab='')
  if(exportFig){dev.off()}
}

################################################################################
# Graphon
################################################################################
uGrid <- vGrid <- seq(0, 1, length.out=100)
phiGraphon <- function(u, v, lambda=1, mu=1){(u^{lambda-1} + v^{mu-1}) * (1 - (u+v-1)^2)}
seed <- 7; set.seed(seed); K <- 3; L <- 4
pi <- rdirichlet(1, rep(1, K)); rho <- rdirichlet(1, rep(1, L)); alpha <-matrix(runif(K*L), K, L)
alpha <- alpha[order(rowSums(alpha)), order(colSums(alpha))]
phiLBM <- function(u, v){
  alpha[(1+max(c(0, which(cumsum(pi)[-K] < u)))), (1+max(c(0, which(cumsum(rho)[-L] < v))))]
}

if(exportFig){png(paste0(dirFig, figName, '-graphon.png'))}
phiGrid <- t(sapply(uGrid, function(u){sapply(vGrid, function(v){phiGraphon(u, v, 2, 3)})}))
phiGridMax <- 2*max(phiGrid); phiGrid <- phiGrid / phiGridMax
plot(wireframe(phiGrid, shade = TRUE, light.source = c(3,3,3), aspect = c(61/87, 0.4), 
     xlab='u', ylab='v', zlab='phi'), xlim=c(0, 1), ylim=c(0, 1), zlim=c(0, 1))
if(exportFig){dev.off()}
# plot(wireframe(phiGrid[order(rowSums(phiGrid)), order(rowSums(phiGrid))], 
#                shade = TRUE, light.source = c(3,3,3), aspect = c(61/87, 0.4), 
#                xlab='u', ylab='v', zlab='phi'))
# m <- nrow(net$adjMat); n <- ncol(net$adjMat)
# U <- runif(m); V <- runif(n);
# Y <- matrix(rbinom(m*n, 1, t(sapply(U, function(u){sapply(V, function(v){phiGraphon(u, v, 2, 3)/phiGridMax})}))),
#             m, n)
image(Y[order(U), order(V)])

if(exportFig){png(paste0(dirFig, figName, '-LBM.png'))}
phiGrid <- t(sapply(uGrid, function(u){sapply(vGrid, function(v){phiLBM(u, v)})}))
plot(wireframe(phiGrid, shade = TRUE, light.source = c(3,3,3), aspect = c(61/87, 0.4), 
               xlab='u', ylab='v', zlab='phi'), xlim=c(0, 1), ylim=c(0, 1), zlim=c(0, 1))
if(exportFig){dev.off()}
# Y <- matrix(rbinom(m*n, 1, t(sapply(U, function(u){sapply(V, function(v){phiLBM(u, v)})}))), m, n)
# image(Y[order(U), order(V)])
