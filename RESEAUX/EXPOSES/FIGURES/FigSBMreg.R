# Figures for SBMreg

rm(list=ls()); par(pch=20)
setwd('/home/robin/Bureau/RECHERCHE/RESEAUX/EXPOSES/FIGURES')
library(sna)

# Tree data
FigName = 'FigSBMreg-Tree'
load('/home/robin/Bureau/RECHERCHE/RESEAUX/GOF-Network/vbemapp/gof-network/pgm/exemples/Tree.Rdata')
n = nrow(Tree$Net); d = dim(Tree$EdgeCovar)[3]

# Edge covariates
par(mfrow=c(2, 2))
names(Tree$EdgeCovar[1, 1, ])
Log = c('T', 'F', 'F')
for (k in 1:d){
   pdf(paste(FigName, '-', names(Tree$EdgeCovar[1, 1, ])[k], '.pdf', sep=''))
   if(Log[k]){
      image(1:n, 1:n, -log(Tree$EdgeCovar[, , k]), xlab='', ylab='', col=grey((1:10)/10))
   }else{
      image(1:n, 1:n, -Tree$EdgeCovar[, , k], xlab='', ylab='', col=grey((1:10)/10))
   }
   dev.off()
}

# Binary network
pdf(paste(FigName, '-BinaryAdjacency.pdf', sep=''))
image(1:n, 1:n, 1-Tree$Net, xlab='', ylab='', col=grey(0:1))
dev.off()

pdf(paste(FigName, '-BinaryNetwork.pdf', sep=''))
gplot(Tree$Net, gmode="graph")
dev.off()

# Distribution of the covariates
par(mfrow=c(2, 2))
for (k in 1:d){
   hist(as.vector(Tree$EdgeCovar[, , k]), breaks=n, main=names(Tree$EdgeCovar[1, 1, ])[k])
}
plot(as.vector(Tree$EdgeCovar[, , 1]), as.vector(Tree$EdgeCovar[, , 3]), log='x')

# Interation vs covariates
par(mfrow=c(2, 2))
for (k in 1:d){
   plot(as.vector(Tree$EdgeCovar[, , k]), as.vector(Tree$Net), xlab=names(Tree$EdgeCovar[1, 1, ])[k])
}

