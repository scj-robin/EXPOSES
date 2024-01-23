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

# Network list
load(paste(DataDir, 'TF_network_list.Rdata', sep=''))
DataNb = length(DataList)
# for (DataName in DataList){
for (i in (2:DataNb)){
   DataName = DataList[i]
   cat(DataName, '')
   
   # Data
#    DataName = DataList[1]
   load(paste(DataDir, DataName, '/', DataName, '.Rdata', sep=''))
   X = Network
   X = X + t(X); X[which(X>1)] = 1
   N = nrow(X)
   
   # SBM + VBEM
   Qmin <- 2
   Qmax <- 10
   nbrepeat <- 20
   lout <- list()
   tab_criteria <- matrix(0, nbrepeat, (Qmax-Qmin+1))
   for(start in 1:nbrepeat) {
         cat(start, '')   
         xout <- mixer(X, qmin=Qmin, qmax=Qmax, method="bayesian", directed=F)
         lout[[start]] <- xout
         tab_criteria[start, ] <- unlist(lapply(xout$output, function(x) x$criterion))
   }
   cat('\n')
   
   tab_criteria <- matrix(0, nbrepeat, (Qmax-Qmin+1))
   for(start in 1:nbrepeat) {
      cat(start, '')   
      xout <- lout[[start]]
      tab_criteria[start, ] <- unlist(lapply(xout$output, function(x) x$criterion))
   }
   cat('\n')
   
   best <- which.max(tab_criteria) 
   if(Qmax-Qmin>0) {
      best_l <- (best-1) %% nbrepeat + 1
      best_c <- (best-1) %/% nbrepeat + 1 
      xout <- lout[[best_l]]
      }else{
         xout <- lout[[best]]
         best_c <- 1
      }
   
   save(X, lout, xout, best_c, file=paste(ResDir, DataName, '/', DataName, '-mixer.Rdata', sep=''))
   print(xout$output[[best_c]]$criterion)
}
