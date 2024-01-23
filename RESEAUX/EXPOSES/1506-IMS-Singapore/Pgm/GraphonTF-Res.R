rm(list=ls())

library(mixer)
library(lattice)
library(sna)
library(igraph)

source("jointPosterior.R")
source("usefulFunctions.R")

start.date = date()

# Dir
DataDir = '/home/robin/Bureau/RECHERCHE/RESEAUX/EXEMPLES/NetBiology/TF_networks/cell6312mmc4/'; 
ResDir = '/home/robin/Bureau/RECHERCHE/RESEAUX/EXEMPLES/NetBiology/TF_networks/cell6312mmc4/'

# Plot parms
L <- 50; U <- seq(0, 1, length=L); V <- seq(0, 1, length=L)

# Network list
load(paste(DataDir, 'TF_network_list.Rdata', sep=''))
for (DataName in DataList){
   cat(DataName, '')
   
   # Data & mixer
   load(paste(DataDir, DataName, '/', DataName, '.Rdata', sep=''))
   load(paste(ResDir, DataName, '/', DataName, '-mixer.Rdata', sep=''))
   X = Network; X = X + t(X); X[which(X>1)] = 1
   N = nrow(X)
   
   # Plot network
   print( gden(X, mode="graph") )
   #png(paste(ResDir, DataName, '-graph.#png', sep=''))
   gplot(X, gmode="graph")
   #dev.off()
   
   # Degrees
   #png(paste(ResDir, DataName, '-degree.#png', sep=''))
   hist(colSums(X), breaks = sqrt(N), main='', xlab='', ylab='')
   #dev.off()
   
   # Graphon
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
   
   # Plot estimated graphon (bayesian SBM)
   jpost <- jointPosterior(U, V, a)
   meanphi <- matrix(0, L, L)
   for(u in 1:L) {
      for(v in u:L){
         meanphi[u, v] <- sum(meanPi[upper.tri(meanPi, diag=TRUE)]*jpost[u, v, ] )
         meanphi[v, u] <- meanphi[u, v]  
      }
   }
   
   #png(paste(ResDir, DataName, '-graphon.#png', sep=''))
   plot(wireframe(meanphi, shade = TRUE, light.source = c(3,3,3), aspect = c(61/87, 0.4)))
   #dev.off()
   
   #png(paste(ResDir, DataName, '-logitgraphon.#png', sep=''))
   plot(wireframe(qlogis(meanphi), shade = TRUE, light.source = c(3,3,3), aspect = c(61/87, 0.4)))
   #dev.off()
   
   #png(paste(ResDir, DataName, '-contour.#png', sep=''))
   contour(U, U, meanphi)
   #dev.off()   
}
