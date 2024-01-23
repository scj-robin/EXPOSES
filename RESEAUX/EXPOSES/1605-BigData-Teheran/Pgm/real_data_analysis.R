
library(mixer)
library(lattice)
library(sna)
library(igraph)

source("jointPosterior.R")
source("usefulFunctions.R")

L<-50
U<-seq(0, 1, length=L)
V<-seq(0, 1, length=L)


QMIN<-2
QMAX<-16
nbrepeat<-30
lout<-list()

data(blog)

id<-which(rowSums(as.matrix(blog$links))==0)

blog$links<-blog$links[-id, -id]
blog$politicalParty<-blog$politicalParty[-id]

X<-as.matrix(blog$links)

tab_criteria<-matrix(0, nbrepeat, (QMAX-QMIN+1))
for(start in 1:nbrepeat) {
  xout<-mixer(X, qmin=QMIN, qmax=QMAX, method="bayesian")
  lout[[start]]<-xout
  tab_criteria[start, ]<-unlist(lapply(xout$output, function(x) x$criterion))
}
best<-which.max(tab_criteria) 
if(QMAX-QMIN>0) {
  best_l<-(best-1) %% nbrepeat + 1
  best_c<-(best-1) %/% nbrepeat + 1 
  xout<-lout[[best_l]]
} else{
  xout<-lout[[best]]
  best_c<-1
}
a<-xout$output[[best_c]]$a
eta<-xout$output[[best_c]]$eta
zeta<-xout$output[[best_c]]$zeta

meanAlpha<-xout$output[[best_c]]$alphas
meanPi<-xout$output[[best_c]]$Pis

ord<-order(meanPi%*%meanAlpha)
#ord<-order(meanAlpha*(meanPi%*%meanAlpha))

meanPi<-meanPi[ord, ord]
meanAlpha<-meanAlpha[ord]
a<-a[ord]

# Plot estimated graphon (SBM)
jfreq<-surfSBM(U%*%t(V), meanAlpha, meanPi)

quartz()
plot(wireframe(jfreq, shade = TRUE, light.source = c(3,3,3), aspect = c(61/87, 0.4)))


# Plot estimated graphon (bayesian SBM)
jpost<-jointPosterior(U, V, a)

meanphi<-matrix(0, L, L)

for(u in 1:L) {
  for(v in u:L){
    meanphi[u, v]<-sum(meanPi[upper.tri(meanPi, diag=TRUE)]*jpost[u, v, ] )
    meanphi[v, u]<-meanphi[u, v]  
  }
  
}

contour(U, U, meanphi)
quartz()
plot(wireframe(meanphi, shade = TRUE, light.source = c(3,3,3), aspect = c(61/87, 0.4),  xlab="", ylab="", zlab=""))
# plot(wireframe((jfreq - phiTrue)^2,shade = TRUE, light.source = c(3,3,3), aspect = c(61/87, 0.4)))
# 
# plot(wireframe((meanphi - phiTrue)^2,shade = TRUE, light.source = c(3,3,3), aspect = c(61/87, 0.4)))
# 
# print(xout$output[[best_c]]$criterion)


Z<-round(t(xout$output[[best_c]]$Taus))
Z<-Z[, ord]

party<-blog$politicalParty
for(i in 1:(best_c+QMIN-1)) {
  print(party[which(Z[, i]==1)])
}

