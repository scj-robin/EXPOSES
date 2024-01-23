library(mixer)
library(lattice)
library(sna)
library(igraph)

source("jointPosterior.R")
source("usefulFunctions.R")

L<-50
U<-seq(0, 1, length=L)
V<-seq(0, 1, length=L)

# Parameters to be adjusted ...

N<-50 # number of nodes
lambda<-1

# rho sur l'enveloppe 
rho<- 1/((lambda+1)^2) 
beta<-1

# rho choisi arbitrairement. (Attention : rho doit rester sous l'enveloppe !!)
#rho<-0.1
#if(rho > 1/((lambda+1)^2)) {
#  stop("rho too large")
#}
#beta<- rho*(lambda+1)^2


# ... Parameters to be adjusted 

#beta<-rho*(lambda/log((1+exp(lambda/2))/(1+exp(-lambda/2))))^2

phiTrue<-surf(U%*%t(V), beta, lambda)


# generate graph ..
z<-runif(N)

X<-matrix(0, N, N)
prob<-surf(z%*%t(z), beta, lambda)
for(i in 1:N) {
  for(j in i:N) {
      if(i != j){
        X[i,j] <- rbinom(1, 1, prob[i, j])
        X[j,i] <- X[i,j]
      }
    }
}

# Compare theoric and estimated density ..
print( rho )
print( gden(X, mode="graph") )

gplot(X, gmode="graph")

# SBM + VBEM
QMIN<-2
QMAX<-10
nbrepeat<-20
lout<-list()
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
plot(wireframe(phiTrue, shade = TRUE, light.source = c(3,3,3), aspect = c(61/87, 0.4)))
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

#contour(U, U, meanphi)
plot(wireframe(meanphi, shade = TRUE, light.source = c(3,3,3), aspect = c(61/87, 0.4)))

plot(wireframe((jfreq - phiTrue)^2,shade = TRUE, light.source = c(3,3,3), aspect = c(61/87, 0.4)))

plot(wireframe((meanphi - phiTrue)^2,shade = TRUE, light.source = c(3,3,3), aspect = c(61/87, 0.4)))

print(xout$output[[best_c]]$criterion)

