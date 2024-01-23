library(PLNmodels)

# Data
dataName <- 'JC2BIM-OaksPLN'
fileName <- paste0(dataName, '.Rdata')
data("oaks")
exportFig <- TRUE
dirFig <- '/home/robin/RECHERCHE/RESEAUX/EXPOSES/FIGURES/'
dirFig <- '.'

# Data selection
sel <- which(oaks$tree != 'intermediate'); selName <- 'noIntermediate'
oaksTmp <- list()
oaksTmp$Abundance <- oaks$Abundance[sel, ]
oaksTmp$Offset <- oaks$Offset[sel, ]
oaksTmp$orientation <- oaks$orientation[sel]
oaksTmp$distTObase <- oaks$distTObase[sel]
oaksTmp$distTOtrunk <- oaks$distTOtrunk[sel]
oaksTmp$distTOground <- oaks$distTOground[sel]
oaksTmp$pmInfection <- oaks$pmInfection[sel]
dataName <- paste0(dataName, '-', selName); fileName <- paste0(dataName, '.Rdata')
oaks <- oaksTmp

# Loading previous results
load(fileName)

# # Fit PLN
# pln0 <- PLN(Abundance ~ offset(log(Offset)) + 1, data=oaks)
# pln1 <- PLN(Abundance ~ offset(log(Offset)) + orientation + distTObase + distTOtrunk + distTOground + pmInfection, data=oaks)
# save(pln0, pln1, file=fileName)

# # Fit PLNnetworks
# net0 <- PLNnetwork(Abundance ~ offset(log(Offset)) + 1, data=oaks, penalties=10^(seq(1, -.5, length.out=50)))
# net1 <- PLNnetwork(Abundance ~ offset(log(Offset)) + orientation + distTObase + distTOtrunk + distTOground+ pmInfection,
#                    data=oaks, penalties=10^(seq(1, -1.5, length.out=50)));
# save(pln0, pln1, net0, net1, file=fileName)

# Stability
# stab0 <- stability_selection(net0)
# stab1 <- stability_selection(net1)
# save(pln0, pln1, net0, net1, stab0, stab1, file=fileName)

# Plots
if(exportFig){png(paste0(dirFig, dataName, '-netSel-null.png'))}
plot(net0)
if(exportFig){dev.off()}
if(exportFig){png(paste0(dirFig, dataName, '-netSel-full.png'))}
plot(net1)
if(exportFig){dev.off()}
plot(net0, 'stability')
plot(net1, 'stability')

net0Best <- getBestModel(net0, "StARS")
plot(net0Best)

net1Best <- getBestModel(net1, "StARS")
plot(net1Best)

# # Tree / tree fit
# for(tree in levels(oaks$tree)){
#    oaksTmp <- oaks[which(oaks$tree==tree), ]
#    fileTmp <- fileName <- paste0(dataName, '-', tree, '.Rdata')
#    load(fileTmp)
#    # PLNA
#    pln0 <- PLN(Abundance ~ offset(log(Offset)) + 1, data=oaksTmp)
#    pln1 <- PLN(Abundance ~ offset(log(Offset)) + orientation + distTObase + distTOtrunk + distTOground + pmInfection, data=oaksTmp)
#    save(pln0, pln1, file=fileTmp)
#    # Network
#    net0 <- PLNnetwork(Abundance ~ offset(log(Offset)) + 1, data=oaksTmp, penalties=10^(seq(1, -.5, length.out=50)))
#    net1 <- PLNnetwork(Abundance ~ offset(log(Offset)) + orientation + distTObase + distTOtrunk + distTOground+ pmInfection,
#                       data=oaksTmp, penalties=10^(seq(1, -1.5, length.out=50)));
#    save(pln0, pln1, net0, net1, file=fileTmp)
#    # Stability
#    stab0 <- stability_selection(net0)
#    stab1 <- stability_selection(net1)
#    save(pln0, pln1, net0, net1, stab0, stab1, file=fileTmp)
# }

# Tree / tree results
tree <- "resistant"  # "susceptible"  "intermediate" "resistant"   
oaksTmp <- oaks[which(oaks$tree==tree), ]
fileTmp <- fileName <- paste0(dataName, '-', tree, '.Rdata')
load(fileTmp)
# Network
plot(net0)
plot(net1)
# Stability
plot(net0, 'stability')
plot(net1, 'stability')
# Plots
net0Best <- getBestModel(net0, "StARS"); plot(net0Best)
net1Best <- getBestModel(net1, "StARS"); plot(net1Best)
