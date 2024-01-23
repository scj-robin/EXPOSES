rm(list=ls())
library(sna)
library(mvtnorm)

##################################################################
# Functions
##################################################################
F_Edge2Mat = function(E){
   p = max(E); A = matrix(0, p, p); for (i in 1:nrow(E)){A[E[i, 1], E[i, 2]]=1}
   return(A + t(A))
}
F_PlotEdge = function(E){p = max(E);
   gplot(F_Edge2Mat(E), gmode='graph', label=(1:p), label.pos=5, vertex.col=0, vertex.cex=2, label.cex=2)
}
F_Mat2Latex = function(M){
   for (i in 1:nrow(M)){cat(M[i, 1]); for (j in (2:ncol(M))){cat(' &', M[i, j])}; cat(' \\\\ \n')}
}
# F_Mat2Latex(M)

##################################################################
# Nb possible graphs
##################################################################
p = c(5, 10, 20, 50)
2^(p*(p-1)/2)
log10(2)*(p*(p-1)/2)

##################################################################
# Trees and non-trees
##################################################################
pdf('FigGGM-Chain.pdf'); F_PlotEdge(rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5))); dev.off()
pdf('FigGGM-Tree.pdf'); F_PlotEdge(rbind(c(1, 2), c(2, 3), c(2, 4), c(1, 5))); dev.off()
pdf('FigGGM-LargerTree.pdf');
F_PlotEdge(rbind(c(1, 2), c(1, 3), c(2, 4), c(2, 5), c(3, 6), c(3, 7), c(4, 8), c(4, 9), c(2, 10))); 
dev.off()
pdf('FigGGM-Star.pdf'); F_PlotEdge(rbind(c(1, 2), c(2, 3), c(2, 4), c(2, 5))); dev.off()
pdf('FigGGM-NonTree.pdf'); F_PlotEdge(rbind(c(1, 2), c(2, 3), c(2, 4), c(1, 5), c(2, 5))); dev.off()
pdf('FigGGM-TreeNonSpanning.pdf'); F_PlotEdge(rbind(c(1, 2), c(2, 3), c(2, 4), c(5, 5))); dev.off()
pdf('FigGGM-Loop.pdf'); F_PlotEdge(rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5, 1))); dev.off()
pdf('FigGGM-FullGraph.pdf'); F_PlotEdge(rbind(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(2, 3), c(2, 4), c(2, 5), c(3, 4), c(3, 5), c(4, 5))); dev.off()

##################################################################
# GGM : complete case example 
##################################################################
# Graph
G = F_Edge2Mat(rbind(c(1, 2), c(2, 3), c(3, 1), c(3, 4))) 
p = nrow(G)
pdf('FigGGM-4nodes.pdf')
gplot(G, gmode='graph', label=(1:p), label.pos=5, vertex.col=0, vertex.cex=2, label.cex=2)
dev.off()
# Covariance matrices
Omega = G + diag(rep(2, p))
diagOmega = diag(1/sqrt(diag(Omega)))
Omega = diagOmega%*%Omega%*%diagOmega
Sigma = solve(Omega)
diagSigma = diag(1/sqrt(diag(Sigma)))
# cholSigma = chol(Sigma)
F_Mat2Latex(G); F_Mat2Latex(round(diagSigma%*%Sigma%*%diagSigma, 2)); F_Mat2Latex(Omega)
# Data
n = 1e2
Y = rmvnorm(n, sigma=Sigma)
# Y = matrix(rnorm(n*p), n, p); Y = Y %*%  cholSigma
C = cor(Y)
Cinv = solve(C); Cinv = diag(1/sqrt(diag(Cinv)))%*%Cinv%*%diag(1/sqrt(diag(Cinv)))
F_Mat2Latex(G); F_Mat2Latex(round(C, 2)); F_Mat2Latex(round(Cinv, 2))

##################################################################
# GGM : incomplete case (ex. 1)
##################################################################
Ec = rbind(c(1, 6), c(2, 6), c(3, 6), c(4, 7), c(5, 7), c(6, 7)); Gc = F_Edge2Mat(Ec) 
Em = rbind(c(1, 2), c(2, 3), c(3, 1), c(4, 5)); Gm = F_Edge2Mat(rbind(Ec, Em)) 
p = 5; r = 2
pdf('FigGGM1-Complete.pdf')
gplot(Gc, gmode='graph', label=(1:(p+r)), label.pos=5, vertex.col=c(rep(0, p), rep(8, r)), vertex.cex=2, label.cex=2)
dev.off()
# Covariance matrices
Oc = Gc + diag(rep(2.5, (p+r))); diagOc = diag(1/sqrt(diag(Oc))); Oc = diagOc%*%Oc%*%diagOc
Sc = solve(Oc); diagSc = diag(1/sqrt(diag(Sc))); Sc = diagSc%*%Sc%*%diagSc
F_Mat2Latex(round(Oc, 1)); F_Mat2Latex(round(Sc, 1))
# 1 missing node
So1 = Sc[1:(p+1), 1:(p+1)]; 
Oo1 = solve(So1); diagOo1 = diag(1/sqrt(diag(Oo1))); Oo1 = diagOo1%*%Oo1%*%diagOo1
F_Mat2Latex(round(Oo1, 1)); F_Mat2Latex(round(So1, 1))
Go1 = 1*(abs(Oo1) > 1e-10)
pdf('FigGGM1-1missing.pdf')
gplot(Go1, gmode='graph', label=(1:(p+1)), label.pos=5, vertex.col=c(rep(0, p), 8), vertex.cex=2, label.cex=2)
dev.off()
# 2 missing nodes
So2 = Sc[1:p, 1:p]; 
Oo2 = solve(So2); diagOo2 = diag(1/sqrt(diag(Oo2))); Oo2 = diagOo2%*%Oo2%*%diagOo2
F_Mat2Latex(round(Oo2, 1)); F_Mat2Latex(round(So2, 1))
Go2 = 1*(abs(Oo2) > 1e-10)
pdf('FigGGM1-2missing.pdf')
gplot(Go2, gmode='graph', label=(1:p), label.pos=5, vertex.col=0, vertex.cex=2, label.cex=2)
dev.off()
# Data
n = 1e3
Y = rmvnorm(n, sigma=Sc); Yo = Y[, 1:p]
Co = cor(Yo)
pdf('FigGGM1-ObservedCor-2missing.pdf')
image(1:p, 1:p, 1-abs(Co), xlab='', ylab='')
dev.off()

##################################################################
# GGM : incomplete case (ex. 2)
##################################################################
Ec = rbind(c(1, 6), c(2, 6), c(3, 6), c(4, 7), c(5, 7), c(6, 7), c(1, 2)); Gc = F_Edge2Mat(Ec) 
# Em = rbind(c(1, 2), c(2, 3), c(3, 1), c(4, 5)); Gm = F_Edge2Mat(rbind(Ec, Em)) 
p = 5; r = 2
# pdf('FigGGM2-Complete.pdf')
gplot(Gc, gmode='graph', label=(1:(p+r)), label.pos=5, vertex.col=c(rep(0, p), rep(8, r)), vertex.cex=2, label.cex=2)
# dev.off()
# Covariance matrices
Oc = Gc + diag(rep(2.5, (p+r))); diagOc = diag(1/sqrt(diag(Oc))); Oc = diagOc%*%Oc%*%diagOc
Sc = solve(Oc); diagSc = diag(1/sqrt(diag(Sc))); Sc = diagSc%*%Sc%*%diagSc
F_Mat2Latex(round(Oc, 1)); F_Mat2Latex(round(Sc, 1))
# 1 missing node
So1 = Sc[1:(p+1), 1:(p+1)]; 
Oo1 = solve(So1); diagOo1 = diag(1/sqrt(diag(Oo1))); Oo1 = diagOo1%*%Oo1%*%diagOo1
F_Mat2Latex(round(Oo1, 1)); F_Mat2Latex(round(So1, 1))
Go1 = 1*(abs(Oo1) > 1e-10)
# pdf('FigGGM2-1missing.pdf')
gplot(Go1, gmode='graph', label=(1:(p+1)), label.pos=5, vertex.col=c(rep(0, p), 8), vertex.cex=2, label.cex=2)
# dev.off()
# 2 missing nodes
So2 = Sc[1:p, 1:p]; 
Oo2 = solve(So2); diagOo2 = diag(1/sqrt(diag(Oo2))); Oo2 = diagOo2%*%Oo2%*%diagOo2
F_Mat2Latex(round(Oo2, 1)); F_Mat2Latex(round(So2, 1))
Go2 = 1*(abs(Oo2) > 1e-10)
# pdf('FigGGM2-2missing.pdf')
gplot(Go2, gmode='graph', label=(1:p), label.pos=5, vertex.col=0, vertex.cex=2, label.cex=2)
# dev.off()
# Data
n = 1e3
Y = rmvnorm(n, sigma=Sc); Yo = Y[, 1:p]
Co = cor(Yo)
# pdf('FigGGM2-ObservedCor-2missing.pdf')
image(1:p, 1:p, 1-abs(Co), xlab='', ylab='')
# dev.off()

