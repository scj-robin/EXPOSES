# Fig of tree-fungi networks

library(igraph)

# Dir
DirData ='~/RECHERCHE/RESEAUX/VALUES/DATA/ArbreChampi-AoAS/'
DirFig = './'
DataName = 'TreeFungiInter'

# Data
Inter = as.matrix(read.table(paste(DirData, DataName, '.txt', sep='')))
NbFungi = dim(Inter)[1]
NbTree = dim(Inter)[2]

# Tree network
Tree = t(Inter)%*%Inter
Tree = Tree - diag(diag(Tree))
for (t in (1:max(Tree))){
   Tmp = matrix((Tree >= t), NbTree, NbTree)
   Sel = which(colSums(Tmp) > 0)
   Tmp = Tmp[Sel, Sel]
   G = graph.adjacency(Tmp, mode='undirected')
   png(paste(DirFig ,'TreeNetwork-t', t, '.png', sep=''))
   plot(G)
   dev.off() #readline()   
}

# Fungi network
Fungi = Inter%*%t(Inter)
Fungi = Fungi - diag(diag(Fungi))
for (t in (1:max(Fungi))){
   Tmp = matrix((Fungi >= t), NbFungi, NbFungi)
   Sel = which(colSums(Tmp) > 0)
   Tmp = Tmp[Sel, Sel]
   G = graph.adjacency(Tmp, mode='undirected')
   png(paste(DirFig ,'FungiNetwork-t', t, '.png', sep=''))
   plot(G)
   dev.off() #readline()   
}
