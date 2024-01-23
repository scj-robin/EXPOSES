rm(list=ls())

library(mixer)
library(lattice)
library(sna)
library(igraph)

source("jointPosterior.R")
source("usefulFunctions.R")

start.date = date()

# Dir
ResDir = '/home/robin/Bureau/1506-IMS-Singapore/Figs/'

# # PPI
# DataDir = '/home/robin/Bureau/RECHERCHE/RESEAUX/EXEMPLES/DIP/'; 
# DataName = 'Dmela20060402' # N = 7067
# DataName = 'Celeg20060402' # N = 2629
# DataName = 'Ecoli20060402' # N = 1838
# DataName = 'Hpylo20060402' # N = 705
# Net = read.table(paste(DataDir, DataName, '/', DataName, '.spm', sep=''))
# N = max(Net)
# X  <-  matrix(0, N, N)
# invisible(sapply(1:nrow(Net), function(i){X[Net[i, 1], Net[i, 2]] <<- 1}))
# X = X + t(X); X[which(X > 1)] = 1; X = X - diag(diag(X))

# # Dir
# FungiTree = load(paste(ResDir, 'FungiTreeNetworks.Rdata', sep=''))
# DataName = 'Fungi'; X = Fungi
# DataName = 'Tree'; X = Tree

# # Fungi Tree
# FungiTree = load(paste(ResDir, 'FungiTreeNetworks.Rdata', sep=''))
# DataName = 'Fungi'; X = Fungi; 
# # DataName = 'Tree'; X = Tree; 

# # Blog
# load('/home/robin/Bureau/BMA/Dropbox-vbemapp/code/mixer/data/blog.rda')
# DataName = 'Blog'
# X = as.matrix(blog$links); X = X+t(X); X[which(X>1)] = 1

# # Macaque
# load('/home/robin/Bureau/BMA/Dropbox-vbemapp/code/mixer/data/macaque.rda')
# DataName = 'Macaque'
# X = as.matrix(blog$links); X = X+t(X); X[which(X>1)] = 1

# Ecoli operon
X = read.table('/home/robin/Bureau/RECHERCHE/RESEAUX/VB-EM/S-Gazal/EColi oriente/data')
DataName = 'EcoliOperon'

N = nrow(X)

# Plot parms
L <- 50; U <- seq(0, 1, length=L); V <- seq(0, 1, length=L)

# Plot network
print( gden(X, mode="graph") )
png(paste(ResDir, DataName, '-graph.png', sep=''))
gplot(X, gmode="graph")
dev.off()

png(paste(ResDir, DataName, '-degree.png', sep=''))
hist(colSums(X), breaks = sqrt(N), main='', xlab='', ylab='')
dev.off()

# SBM + VBEM results
load(paste(ResDir, DataName, '-mixer.Rdata', sep=''))
a <- xout$output[[best_c]]$a
eta <- xout$output[[best_c]]$eta
zeta <- xout$output[[best_c]]$zeta
meanAlpha <- xout$output[[best_c]]$alphas
meanPi <- xout$output[[best_c]]$Pis
ord <- order(meanPi%*%meanAlpha)
#ord <- order(meanAlpha*(meanPi%*%meanAlpha))

meanPi <- meanPi[ord, ord]
meanAlpha <- meanAlpha[ord]
a <- a[ord]

# Plot estimated graphon (SBM)
# jfreq <- surfSBM(U%*%t(V), meanAlpha, meanPi)
# plot(wireframe(jfreq, shade = TRUE, light.source = c(3,3,3), aspect = c(61/87, 0.4)))

# Plot estimated graphon (bayesian SBM)
jpost <- jointPosterior(U, V, a)
meanphi <- matrix(0, L, L)
for(u in 1:L) {
  for(v in u:L){
    meanphi[u, v] <- sum(meanPi[upper.tri(meanPi, diag=TRUE)]*jpost[u, v, ] )
    meanphi[v, u] <- meanphi[u, v]  
  }
}

png(paste(ResDir, DataName, '-graphon.png', sep=''))
plot(wireframe(meanphi, shade = TRUE, light.source = c(3,3,3), aspect = c(61/87, 0.4)))
dev.off()

png(paste(ResDir, DataName, '-logitgraphon.png', sep=''))
plot(wireframe(qlogis(meanphi), shade = TRUE, light.source = c(3,3,3), aspect = c(61/87, 0.4)))
dev.off()

png(paste(ResDir, DataName, '-contour.png', sep=''))
contour(U, U, meanphi)
dev.off()

save(lout, xout, file=paste(ResDir, DataName, '-mixer.Rdata', sep=''))

print(xout$output[[best_c]]$criterion)

stop.date = date()
print(start.date)
print(stop.date)
