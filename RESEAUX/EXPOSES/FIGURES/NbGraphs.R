# Nb spanninng trees for p nodes
# nb part. univer = 10^80 - 10^85

# Parms
rm(list=ls()); 
nbPart = 10^82.5
p.list = (1:100); p.max = max(p.list)
s.list = p.list*log(p.list)

# Functions
F_NbDAGs <- function(p){
   a = loga = rep(0, p+1); log2 = log(2)
   a[1:2] = 1; loga[1:2] = 0
   for (j in 2:p){
      for (k in 1:j){
         a[j+1] = a[j+1] + (-1)^(k-1) * exp(lchoose(j, k) + (k*(j-k))*log2 + loga[j-k+1])
      }
      loga[j+1] = log(a[j+1])
   }
   return(a)
}

nbDAGs = F_NbDAGs(p.max)[-1]
nbSpanTrees = (p.list)^(p.list-2)
nbUndirGraphs = 2^(p.list*(p.list-1)/2)
nbSparseGraphs = choose(p.list*(p.list-1)/2, s.list)

pdf(file='NbGraphs.pdf')
par(pch=20, lwd=3, mfrow=c(1, 1), cex=1.5)
plot(p.list, nbDAGs, col=2, log='y', type='l', xlab='number of nodes', ylab='')
lines(p.list, nbUndirGraphs, col=3)
lines(p.list, nbSpanTrees, col=4)
lines(p.list, nbSparseGraphs, col=5)
abline(h=nbPart, lty=1, col=8)
dev.off()

pdf(file='NbGraphs2.pdf')
par(pch=20, lwd=3, mfrow=c(1, 1), cex=1.5)
plot(p.list, nbDAGs, col=2, log='y', type='l', xlab='number of nodes', ylab='')
lines(p.list, nbUndirGraphs, col=3)
lines(p.list, nbSpanTrees, col=4)
abline(h=nbPart, lty=1, col=8)
dev.off()
