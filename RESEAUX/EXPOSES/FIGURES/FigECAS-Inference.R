# Figures for ECAS / SFdS / NAIM : reminder on prob. dist

library(sna); library(mvtnorm); library(huge); library(LITree); library(ROCR); library(ape)
library(mixtools); library(latex2exp); library(EMtree)
source('FigECAS-Functions.R')
source('/home/robin/PgmR/General/FunctionsMatVec.R')
source('/home/robin/RECHERCHE/RESEAUX/S-Founas/NoisyNetworkInference/Pgm/Functions/SimulFunctions.R')
source('/home/robin/RECHERCHE/RESEAUX/S-Founas/NoisyNetworkInference/Pgm/Functions/ResultFunctions.R')

###############################################################################
# L0, L1, L2
alphaList = c(.8, .6, .4)
beta = c(1.5, 3); Sigma = matrix(c(1.1, 1.3, 1.3, 2), 2, 2)
plus.cex = 4; beta.cex = 3; lwd = 4

PlotBetaHat <- function(){
   plot(c(0, 0), xlim=c(-1.5, 4), ylim=c(-1.5, 5), col=0, axes=FALSE, cex.lab=2, xlab='', ylab='')
        # xlab=TeX('$\\beta_1'), ylab=TeX('$\\beta_2')); 
   abline(h=0, v=0)
   points(beta[1], beta[2], pch='+', cex=plus.cex); 
   text(beta[1], beta[2], pos=4, labels=TeX('$\\widehat{\\beta}'), cex=beta.cex)
   text(3.5, 0.15, pos=4, labels=TeX('$\\beta_1'), cex=beta.cex)
   text(0, 4.5, pos=4, labels=TeX('$\\beta_2'), cex=beta.cex)
   sapply(alphaList, function(alpha){
      ellipse(beta, Sigma, alpha=alpha, col=8, lwd=lwd, lty=2, type='l', xlab='', ylab='')
   })
}

pdf('FigECAS-RegularizationL0.pdf'); par(mex=.01)
PlotBetaHat()
abline(h=0, v=0, col=4, lwd=lwd)
ellipse(beta, Sigma, alpha=.37, col=8, lwd=lwd, type='l', xlab='', ylab='')
points(0, 1.2, col=2, pch='+', cex=plus.cex); 
dev.off()

pdf('FigECAS-RegularizationL1.pdf'); par(mex=.01)
PlotBetaHat()
l1Ball = matrix(c(1, 0, -1, 0, 1, 0, 1, 0, -1, 0), 5, 2)
lines(l1Ball, col=4, lwd=lwd)
ellipse(beta, Sigma, alpha=.34, col=8, lwd=lwd, type='l', xlab='', ylab='')
points(0, 1, col=2, pch='+', cex=plus.cex); 
dev.off()

pdf('FigECAS-RegularizationL2.pdf'); par(mex=.01)
PlotBetaHat()
alpha = 2*acos(-1)*seq(0, 1, by=.01); 
l2Ball = cbind(cos(alpha), sin(alpha))
lines(l2Ball, col=4, lwd=lwd)
ellipse(beta, Sigma, alpha=.37, col=8, lwd=lwd, type='l', xlab='', ylab='')
points(0.2, .95, col=2, pch='+', cex=plus.cex); 
dev.off()

pdf('FigECAS-TreeGM.pdf'); par(mex=.01)
p = 10; seed = 1; set.seed(seed)
Gtree = mst(as.matrix(dist(matrix(rnorm(2*p), p, 2), diag=TRUE)))
gplot(Gtree, gmode='graph', label=1:p, label.pos=5, vertex.col=0, vertex.cex=2, label.cex=2)
dev.off()

###############################################################################
# # Glasso performances for varying p and n
# # Simulation
# pList = round(10^c(1, 1.5, 2)); nList = pList; 
# seed = 3; set.seed(seed); scoreList = list(); sim = 0
# for (p in pList){
#    density = 1.5*log(p)/p
#    G = F_Vec2Sym(F_Sym2Vec(matrix(rbinom(p^2, 1, density), p, p)))
#    OmegaSigma = MakeOmegaSigma(G)
#    Y = rmvnorm(max(nList), sigma=OmegaSigma$Sigma)
#    for (n in nList){
#       sim = sim+1
#       cat('p =', p, 'n =', n, ': ')
#       scoreGlasso = fitHuge(Y[1:n, ], method='glasso')
#       print(dim(scoreGlasso))
#       scoreList[[sim]] = list(p=p, n=n, G=G, Glasso = scoreGlasso)
#       # scoreMB = fitHuge(Y[1:n, ], method='mb')
#       # scoreList[[sim]]$MB =  scoreMB
#       weightTree = -n/2*log(1 - cor(Y)^2)
#       scoreTree = EdgeProba(weightTree)
#       scoreList[[sim]]$Tree =  scoreTree
#       save(scoreList, file='FigECAS-scoreList.Rdata')
#    }
# }

# Plot
lwd = 2; 
load('FigECAS-scoreList.Rdata'); 
# pdf('FigECAS-NversusM.pdf')
pdf('FigECAS-NversusM-Glasso.pdf')
simNb = length(scoreList); 
par(mfrow=c(3, 3), mex=.5)
for (sim in 1:simNb){
   gLassoROC = performance(prediction(F_Sym2Vec(scoreList[[sim]]$Glasso), F_Sym2Vec(scoreList[[sim]]$G)), 
                           "tpr", "fpr")
   plot(gLassoROC, lwd=lwd, xlim=c(0, 1), ylim=c(0, 1), cex.lab=2, xlab='', ylab='', col=1, 
        main=paste0('p = ', scoreList[[sim]]$p, ',  n = ', scoreList[[sim]]$n), type='s'); abline(0, 1, lty=2)
   if(is.element("MB", names(scoreList[[sim]]))){
      MBROC = performance(prediction(F_Sym2Vec(scoreList[[sim]]$MB), F_Sym2Vec(scoreList[[sim]]$G)), 
                          "tpr", "fpr")
      plot(MBROC, lwd=lwd, col=2, , type='s', add=TRUE)
   }
   # if(is.element("Tree", names(scoreList[[sim]]))){
   #    TreeROC = performance(prediction(F_Sym2Vec(scoreList[[sim]]$Tree), F_Sym2Vec(scoreList[[sim]]$G)),
   #                          "tpr", "fpr")
   #    plot(TreeROC, lwd=lwd, col=4, , type='s', add=TRUE)
   # }
}
dev.off()
