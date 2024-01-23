# Example from tree-fungus network

rm(list=ls());
library(sna); library(blockmodels)

exportPlot <- FALSE; mex <- .9
par(mfrow=c(1, 1), mex=mex)

# Data
dataDir <- '/home/robin/Bureau/RECHERCHE/RESEAUX/GOF-Network/vbemapp/gof-network/Data/FungiTree/'
dataName <- 'TreeFungi'
figDir <- '/home/robin/RECHERCHE/RESEAUX/EXPOSES/FIGURES/'

# # Make Rdata
# adj <- read.table(paste0(dataDir, 'SuppMatBin/csv/', dataName, 'Inter.csv'), sep=',', header=TRUE)
# adj <- adj[-nrow(adj) ,-ncol(adj)]
# fungiName <- adj[, 1]; adj <- adj[, -1]; 
# treeName <- colnames(adj)
# adj <- as.matrix(adj); colnames(adj) <- treeName; rownames(adj) <- fungiName
# save(treeName, fungiName, adj, file=paste0(dataDir, dataName, '.Rdata'))
load(paste0(dataDir, dataName, '.Rdata'))

# Network
m <- nrow(adj); n <- ncol(adj)
image(1:m, 1:n, adj, xlab='fungi', ylab='trees')

# Blockmodel
# LBM  <- BM_bernoulli('LBM', adj=adj)
# LBM$estimate()
# save(LBM, file=paste0(dataDir, dataName, '-LBM.Rdata'))
load(paste0(dataDir, dataName, '-LBM.Rdata'))

# Parms
bestNum <- which.max(LBM$ICL); 
tauRow <- LBM$memberships[[bestNum]]$Z1; K <- ncol(tauRow); 
tauRow <- tauRow[, order(colMeans(tauRow))]
tauCol <- LBM$memberships[[bestNum]]$Z2; L <- ncol(tauCol); 
tauCol <- tauCol[, order(colMeans(tauCol))]
orderRow <- order(tauRow%*%(1:K)); orderCol <- order(tauCol%*%(1:L))
Zrow <- apply(tauRow, 1, which.max); Zcol <- apply(tauCol, 1, which.max)
Nrow <- table(Zrow); Ncol <- table(Zcol)

# Classif
orderAdj <- adj[orderRow, orderCol]
image(1:m, 1:n, orderAdj, xlab='fungi', ylab='trees')
abline(h=.5+cumsum(Ncol), v=.5+cumsum(Nrow), col=4, lwd=2)

# Network
posRow <- cbind(seq(0, 1, length.out=m), rep(1, m))
posCol <- cbind(seq(0, 1, length.out=n), rep(0, n))
plot(0, xlim=c(0, 1), ylim=c(0, 1), col=0, xlab='', ylab='')
points(posRow, col=(Zrow[orderRow]), pch=20)
points(posCol, col=(Zcol[orderCol]), pch=20)
for(i in 1:m){for(j in 1:n){
   if(orderAdj[i, j]==1){
      lines(c(posRow[i, 1], posCol[j, 1]), c(posRow[i, 2], posCol[j, 2]), col=8)
   }
}}
           