# Choleski transform

rm(list=ls()); library(mvtnorm); library(ellipse); palette('R3')
seed <- 1; set.seed(seed)
exportFig <- TRUE
figName <- 'Cholevski'

# Parms
n <- 1e3
sigma <- c(2, 1); rho <- 0.75
Sigma <- diag(sigma)%*%matrix(c(1, rho, rho, 1), 2, 2)%*%diag(sigma)

# Choleski
eigSigma <- eigen(Sigma)
eigSigma$vectors
eigSigma$values
cholSigma <- (eigSigma$vectors)%*%diag(1/sqrt(eigSigma$values))%*%t(eigSigma$vectors)

# Data
x <- rmvnorm(n, sigma=Sigma)
cov(x)
xx <- x%*%cholSigma
cov(xx)

# Plots
if(exportFig){png(paste0(figName, '-raw.png'))}
plot(x, xlab='', ylab='', pch=20, xlim=c(-5, 5), ylim=c(-5, 5), axes=0, col=8)
abline(v=0, h=0, lty=2, lwd=2)
lines(ellipse(x=rho, scale=sigma), col=2, lwd=3)
if(exportFig){dev.off()}

# Plots
if(exportFig){png(paste0(figName, '-chol.png'))}
plot(xx, xlab='', ylab='', pch=20, xlim=c(-5, 5), ylim=c(-5, 5), axes=0, col=8)
abline(v=0, h=0, lty=2, lwd=2)
lines(ellipse(x=0, scale=c(1, 1)), col=2, lwd=3)
if(exportFig){dev.off()}
