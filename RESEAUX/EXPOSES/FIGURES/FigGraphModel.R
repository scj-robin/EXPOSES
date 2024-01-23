rm(list=ls())

##################################################################
# Functions
##################################################################
F_PlotNet <- function(NodePos, NodeColor, EdgeList, EdgeColor){
   p = nrow(NodePos); m = nrow(EdgeList); coef.lim = 1.05
   plot(NodePos[ ,1], NodePos[, 2], pch=20, xlab='', ylab='', axes=F,
        xlim=c(coef.lim*min(NodePos[, 1]), coef.lim*max(NodePos[, 1])), 
        ylim=c(coef.lim*min(NodePos[, 2]), coef.lim*max(NodePos[, 2])),
        col=0, cex=10)
   for (i in 1:m){lines(NodePos[EdgeList[i ,], 1], NodePos[EdgeList[i ,], 2], lwd=4, col=EdgeColor[i])}
   points(NodePos[ ,1], NodePos[, 2], pch=20, col=NodeColor, cex=12)
   text(NodePos[ ,1], NodePos[, 2], labels=(1:p), cex=3)
}

###################################################################
# Omitting a node
###################################################################
NodePos = rbind(c(-1, -1), c(-1, 1), c(-1/5, 1/5), c(1, 1), c(1, -1)); 
EdgeList0 = rbind(c(1, 2), c(2, 3), c(1, 3), c(3, 4), c(4, 5)); 
p = nrow(NodePos); m = nrow(EdgeList0)
NodeColor0 = rep(4, p); EdgeColor0 = rep(1, m)
pdf('FigGraphModel-5nodes.pdf')
F_PlotNet(NodePos, NodeColor0, EdgeList0, EdgeColor0)
dev.off()

EdgeListCond5 = rbind(c(1, 2), c(2, 3), c(1, 3), c(3, 4)); 
pdf('FigGraphModel-5nodesMarg5.pdf')
EdgeListMarg5 = EdgeListCond5 
F_PlotNet(NodePos, c(4, 4, 4, 4, 2), rbind(EdgeList0, EdgeListMarg5), c(rep(8, m), rep(1, m)))
dev.off()
pdf('FigGraphModel-5nodesCond5.pdf')
F_PlotNet(NodePos, c(4, 4, 4, 4, 2), rbind(EdgeList0, EdgeListCond5), c(rep(8, m), rep(1, m)))
dev.off()

pdf('FigGraphModel-5nodesMarg4.pdf')
EdgeListMarg4 = rbind(c(1, 2), c(2, 3), c(1, 3), c(3, 5)); 
F_PlotNet(NodePos, c(4, 4, 4, 2, 4), rbind(EdgeList0, EdgeListMarg4), c(rep(8, m), rep(1, m)))
dev.off()
pdf('FigGraphModel-5nodesCond4.pdf')
EdgeListCond4 = rbind(c(1, 2), c(2, 3), c(1, 3)); 
F_PlotNet(NodePos, c(4, 4, 4, 2, 4), rbind(EdgeList0, EdgeListCond4), c(rep(8, m), rep(1, m)))
dev.off()

pdf('FigGraphModel-5nodesMarg3.pdf')
EdgeListMarg3 = rbind(c(1, 2), c(1, 4), c(2, 4), c(4, 5)); 
F_PlotNet(NodePos, c(4, 4, 2, 4, 4), rbind(EdgeList0, EdgeListMarg3), c(rep(8, m), rep(1, m)))
dev.off()
pdf('FigGraphModel-5nodesCond3.pdf')
EdgeListCond3 = rbind(c(1, 2), c(4, 5)); 
F_PlotNet(NodePos, c(4, 4, 2, 4, 4), rbind(EdgeList0, EdgeListCond3), c(rep(8, m), rep(1, m)))
dev.off()
