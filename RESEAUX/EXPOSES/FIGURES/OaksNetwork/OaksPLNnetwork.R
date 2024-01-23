rm(list=ls())
dirFig <- '/home/robin/RECHERCHE/RESEAUX/EXPOSES/FIGURES/'
exportFig <- TRUE

library(PLNmodels)
dataName <- 'Oaks'
data(oaks)
attach(oaks)
p <- ncol(Abundance); n <- nrow(Abundance)
speciesSel <- 1:p

# # Data reduction
# p0 <- 21
# plot(sort(colMeans(Abundance/Offset)), log='y')
# abline(h=mean(Abundance[, ncol(Abundance)]/Offset[, ncol(Abundance)]))
# speciesRank <- rank(colMeans(Abundance/Offset))
# abline(h=mean(Abundance[, which(speciesRank==p-20)]/Offset[, which(speciesRank==p-20)]))
# 
# speciesSel <- which(speciesRank > p-p0)
# speciesCol <- rep(1, p); speciesCol[speciesSel] <- 2
# plot(colMeans(Abundance/Offset), log='y', col=speciesCol, pch=p0)
# abline(h=mean(Abundance[, which(speciesRank==p-p0)]/Offset[, which(speciesRank==p-p0)]))
# Abundance <- Abundance[, speciesSel]
# Offset <- Offset[, speciesSel]
# p <- ncol(Abundance)
# dataName <- paste0(dataName, '-p', p)

# PLN
plnFile <- paste0(dataName, '-pln.Rdata')
if(!file.exists(plnFile)){
  pln <- PLN(Abundance ~ tree + distTOground + orientation + offset(log(Offset)))
  save(pln, file=plnFile)
}else{load(plnFile)}

# Net
netFile <- paste0(dataName, '-net.Rdata')
if(!file.exists(netFile)){
  net <- PLNnetwork(Abundance ~ tree + distTOground + orientation + offset(log(Offset)))
  lambda <- exp(seq(log(max(net$penalties/2)), log(min(net$penalties/10)), length.out=50))
  net <- PLNnetwork(Abundance ~ tree + distTOground + orientation + offset(log(Offset)),
                    penalties=lambda)
  plot(net)
  # Stability selection
  subs <- replicate(20, sample.int(n, size = n/2), simplify = FALSE)
  stability_selection(net, subsamples = subs)
  save(net, file=paste0(dataName, '-net.Rdata'))
}else{load(netFile)}

# Best network
bestNet <- getBestModel(net, "StARS")
plot(net, 'stability')
if(exportFig==TRUE){png(paste0(dirFig, dataName, '-network.png'))}
plot(bestNet)
if(exportFig==TRUE){dev.off()}

# Tree by tree
treeList <- levels(tree); treeNb <- length(treeList)
treeNetList <- bestNetList <- list()
detach(oaks)
for(tt in 1:treeNb){
  data(oaks); attach(oaks)
  sampleSel <- which(tree==treeList[tt])
  # Selection
  speciesRank <- rank(colMeans(Abundance[sampleSel, ]/Offset[sampleSel, ]))
  speciesSel <- which(speciesRank > p-p0)
  # Network inference
  nTmp <- length(sampleSel)
  netTmp <- PLNnetwork(Abundance[sampleSel, ] ~ distTOground[sampleSel] + 
                         orientation[sampleSel] + offset(log(Offset[sampleSel, ])))
  lambda <- exp(seq(log(max(netTmp$penalties/2)), log(min(netTmp$penalties/20)), length.out=50))
  treeNetList[[tt]] <- PLNnetwork(Abundance[sampleSel, ] ~ distTOground[sampleSel] + 
                             orientation[sampleSel] + offset(log(Offset[sampleSel, ])),
                    penalties=lambda)
  plot(treeNetList[[tt]])
  # Selection
  subs <- replicate(20, sample.int(nTmp, size = nTmp/2), simplify = FALSE)
  stability_selection(treeNetList[[tt]], subsamples = subs)
  # Plot
  bestNetList[[tt]] <- getBestModel(treeNetList[[tt]], "StARS")
  plot(treeNetList[[tt]], 'stability')
  if(exportFig==TRUE){png(paste0(dirFig, dataName, '-network-', treeList[tt], '.png'))}
    plot(bestNetList[[tt]])
  if(exportFig==TRUE){dev.off()}
}

