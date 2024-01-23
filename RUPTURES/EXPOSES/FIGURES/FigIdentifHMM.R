# Figure for talk on HMM identifiability

source('/home/robin/PgmR/FunctionsHMMGeneral.R')
library(lattice)

# Plot parms
Lwd = 3
par(mfrow=c(1, 1))

# Parms
n = 1e4
K = 3; Ku = 2
Pi = matrix(1, K, K) +diag(rep(100, K))
Pi = Pi / rowSums(Pi)
sigma = 1
mu = 3*(1:(Ku*K))

# Simulation
S = SimMC(n, rep(1/K, K), Pi) %*% (1:K)
U = ceiling(Ku*runif(n)) 
Z = S+(U-1)*K
Y = rnorm(n, mean=mu[Z], sd=sigma)

# Densities
Ysort = sort(Y)
phi = matrix(0, n, K*Ku)
for (d in (1:(K*Ku))){
  phi[, d] = dnorm(Ysort, mean=mu[d], sd=sigma)
}
g = phi %*% (rep(1/K/Ku, K*Ku))

# 1d plots
plot(Y, col=S, pch=20, cex=.5, main='', xlab='', ylab='')
# same color 
H = hist(Y, breaks=ceiling(sqrt(n)), main='', xlab='', ylab='')
w = mean(diff(H$breaks))
for (d in (1:(K*Ku))){
  lines(Ysort, w*n*phi[, d]/K/Ku, col=1, lwd=Lwd)
}
lines(Ysort, w*n*g, col=1, lwd=Lwd, lty=2)
# wrong colors 
H = hist(Y, breaks=ceiling(sqrt(n)), main='', xlab='', ylab='')
w = mean(diff(H$breaks))
for (d in (1:(K*Ku))){
  lines(Ysort, w*n*phi[, d]/K/Ku, col=ceiling(d/Ku), lwd=Lwd)
}
lines(Ysort, w*n*g, col=1, lwd=Lwd, lty=2)
# true colors 
H = hist(Y, breaks=ceiling(sqrt(n)), main='', xlab='', ylab='')
w = mean(diff(H$breaks))
for (d in (1:(K*Ku))){
  lines(Ysort, w*n*phi[, d]/K/Ku, col=1+(d%%K), lwd=Lwd)
}
lines(Ysort, w*n*g, col=1, lwd=Lwd, lty=2)

# 3d plots
# same color 
cloud(Y[3:n] ~ Y[1:(n-2)] * Y[2:(n-1)], col=1, pch=20, cex=.5, 
      default.scales=list(arrows=F), main='', xlab='Y(t-1)', ylab='Y(t)', zlab='Y(t+1)')
# wrong colors 
Y3 = cbind(Y[1:(n-2)], Y[2:(n-1)], Y[3:n])
W = ceiling(Z/Kz)
cloud(Y[3:n] ~ Y[1:(n-2)] * Y[2:(n-1)], col=W[2:(n-1)], pch=20, cex=.5, 
      default.scales=list(arrows=F), main='', xlab='Y(t-1)', ylab='Y(t)', zlab='Y(t+1)')
# true colors 
cloud(Y[3:n] ~ Y[1:(n-2)] * Y[2:(n-1)], col=S[2:(n-1)], pch=20, cex=.5, 
      default.scales=list(arrows=F), main='', xlab='Y(t-1)', ylab='Y(t)', zlab='Y(t+1)')
# wrong conditioning 
cloud(Y[3:n] ~ Y[1:(n-2)] * Y[2:(n-1)] | as.factor(W[2:(n-1)]), 
      col=1, pch=20, cex=.5, 
      default.scales=list(arrows=F), main='', xlab='Y(t-1)', ylab='Y(t)', zlab='Y(t+1)')
# true conditioning 
cloud(Y[3:n] ~ Y[1:(n-2)] * Y[2:(n-1)] | as.factor(S[2:(n-1)]), 
      col=1, pch=20, cex=.5, 
      default.scales=list(arrows=F), main='', xlab='Y(t-1)', ylab='Y(t)', zlab='Y(t+1)')

# Trnasition matrices
PiW = matrix(0, K, K); PiS = PiW
for (k in (1:K)){
  for (l in (1:K)){
    PiW[k, l] = sum((W[1:(n-1)]==k) & (W[2:n]==l))
    PiS[k, l] = sum((S[1:(n-1)]==k) & (S[2:n]==l))
  }
}
PiW = PiW / rowSums(PiW)
PiS = PiS / rowSums(PiS)
