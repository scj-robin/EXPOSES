# Figure Beta mixture identifiability
rm(list=ls())

#Dir = '/media/donnees/RECHERCHE/RUPTURES/Exposes/Figures/'

# Parms
n = 5e3
parms = t(matrix(c(1.2, 10, 5, 5, 10, 1.2), 2, 3))
K = dim(parms)[1]
prop = c(3, 2, 3); prop = prop/sum(prop)

# Densities
x = seq(0, 1, by=.001)
phi = matrix(0, length(x), K)
for (k in (1:K)){
   phi[, k] = dbeta(x, parms[k, 1], parms[k, 2])
}

# Simul
Z = (1:K)%*%rmultinom(n, 1, prop)
X = rbeta(n, parms[Z, 1], parms[Z, 2])
png('FigIndetif-BetaMixt0.png')
H = hist(X, breaks=sqrt(n), main='', xlab='', ylab='')
dev.off()
w = mean(diff(H$breaks))

# Plots
F_PlotBetaMixt <- function(X, x, p){
   D = dim(p)[2]
   f = phi %*% p
   hist(X, breaks=sqrt(n), main='', xlab='', ylab='', yaxt='n', , xaxt='n')
   for(d in (1:D)){
      lines(x, n*w*f[, d], col=(1+d), lwd=4)
   }   
   lines(x, n*w*rowSums(f), col=1, lwd=4, lty=2)
}

png('FigIndetif-BetaMixt1.png')
p = diag(prop)
F_PlotBetaMixt(X, x, p)
dev.off()

png('FigIndetif-BetaMixt2.png')
p = matrix(c(prop[1], prop[2], 0, 0, 0, prop[3]), 3, 2)
F_PlotBetaMixt(X, x, p)
dev.off()

png('FigIndetif-BetaMixt3.png')
p = matrix(c(prop[2]/2, prop[2], prop[2]/2, prop[1]-prop[2]/2, 0, prop[3]-prop[2]/2), 3, 2)
F_PlotBetaMixt(X, x, p)
dev.off()

