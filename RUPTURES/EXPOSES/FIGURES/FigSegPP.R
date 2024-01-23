# Figures for Poisson point process segmentation

rm(list=ls())
source('/home/robin/RECHERCHE/RUPTURES/Hawkes/segHP/codes/PoissonCV-R/FunctionSegPP-CV.R')
source('/home/robin/RECHERCHE/RUPTURES/Hawkes/segHP/codes/Poisson-R/DynProg.R')
source('/home/robin/RECHERCHE/RUPTURES/Hawkes/segHP/codes/Poisson-R/Gsegmentation.R')
source('/home/robin/RECHERCHE/RUPTURES/Hawkes/segHP/codes/Poisson-R/G.contrast.R')

# Dirs
figName <- 'FigSegPP'
figDir <- './'
exportFig <- FALSE

################################################################################
# Examples
################################################################################
# Simul
seed <- 1; set.seed(seed)
K <- 5
tau <- c(0, sort(runif(K-1)), 1)
lambda <- rgamma(K, shape=.9, rate=.009)
Kmax <- 10; Kopt <- K
times <- c(); deltaTau <- diff(tau)
n <- length(times)
dataName <- paste0('simul-K', K, '-n', n, '-seed', seed)
Lambda <- 0
for(k in 1:K){Lambda <- c(Lambda, Lambda[k]+deltaTau[k]*lambda[k])}
for(k in 1:K){times <- c(times, tau[k]+deltaTau[k]*sort(runif(rpois(1, lambda[k]*deltaTau[k]))))}

# Data: chiropteres : s = 2295
dataDir <- '/home/robin/RECHERCHE/RUPTURES/Hawkes/DataSR/chiroptères/'
load(paste0(dataDir, 'sequences.Rdata'))
seqNum <- 2295; seq <- seqList[[seqNum]];
times <- seq$times
dataName <- paste0('Chiroptere-seq', seqNum, '-day', seq$day); Kmax <- 10; Kopt <- 5

# Data: vulcanos
dataDir <- '/home/robin/RECHERCHE/RUPTURES/Hawkes/segHP/data/vulcano/'
dataName <- 'Kilauea'; Kmax <- 10; Kopt <- 4
# dataName <- 'MaunaLoa'; Kmax <- 10; Kopt <- 6
times <- as.vector(read.table(paste0(dataDir, dataName, '-PP.csv'), sep=';')$x)
 
# Segmentation
n <- length(times); a <- 1; b <- 1/n
breaks <- TimesProc(times)
# cost <- CostPP(breaks$N.Potential.Times, breaks$Potential.Times, a=a, b=b) # Cost matrix
cost <- Gsegmentation(breaks$N.Potential.Times, breaks$Potential.Times, b=b) # Cost matrix
seg <- SegProc(times, Kmax, breaks, cost)

# Plots
if(exportFig){png(paste0(figDir, figName, '-', dataName, '-Path.png'))}
plot(c(0, times, 1), c((0:n), n), type='s', xlab='t', ylab='N(t)')
if(exportFig){dev.off()}

if(exportFig){png(paste0(figDir, figName, '-', dataName, '-Seg.png'))}
PlotSeg(times, seg[[Kopt]])
if(exportFig){dev.off()}

################################################################################
# Contrast
################################################################################
# Simul
seed <- 1; set.seed(seed)
# K <- 3; n <- 20; tau <- sort(runif(K-1)) 
# lambdak <- 5+2*(-1)^(1:K); nk <- rpois(K, lambdak); n <- sum(nk)
# times <- tau[1]*sort(runif(nk[1]))
# for(k in 2:K){times <- c(times, tau[k-1]+(tau[k]-tau[k-1])*sort(runif(nk[k])))}
# dataName <- paste0('simul-n', n, '-K', K, '-seed', seed)
n <- 10; times <- sort(runif(n)); tau <- lambdak <- nk <- NULL
dataName <- paste0('FigSegPPP-simul-n', n, '-K1-seed', seed)

# Simul
if(exportFig){png(paste0(figDir, dataName, '-Path.png'))}
plot(c(0, times, 1), c(0, 1:n, n), type='s', xlab='', ylab='N(t)', lwd=2)
abline(v=times, col=8, lwd=2)
abline(v=c(0, 1), col=1, lwd=3)
if(exportFig){dev.off()}

# Grid
nGrid <- 1e3; tauGrid <- seq(0, 1, length.out=nGrid)
Ngrid <- rep(0, nGrid)
for(i in 1:n){Ngrid[which(tauGrid>=times[i])] <- Ngrid[which(tauGrid>=times[i])]+1}
lambdaGrid <- rep(lambdak[1], nGrid)
for(k in 2:K){lambdaGrid[which((tauGrid>=tau[k-1]) & (tauGrid<=tau[k]))] <- lambdak[k]}

# Function
LogL <- function(dN, dT, b=1e-3, a=n*b){
  (a+dN)*(log(a+dN+(a+dN==0)) - log(b+dT + (b+dT==0)))
}

# Contrast
a <- b <- 0;
# b <- 1/n; a <- n*b
logL1 <- LogL(Ngrid, tauGrid, a=a, b=b)
logL3 <- LogL(n-Ngrid, 1-tauGrid, a=a, b=b)
logLtau <- logL1%o%rep(1, nGrid) + rep(1, nGrid)%o%logL3
logLtau[lower.tri(logLtau)] <- NA
for(i in 1:nGrid){
  for(j in i:nGrid){
    logLtau[i, j] <- logLtau[i, j] + LogL(Ngrid[j]-Ngrid[i], tauGrid[j]-tauGrid[i], a=a, b=b)
  }
}
library(fields)
if(exportFig){png(paste0(figDir, 'FigSegPP-', dataName, '-Contrast.png'))}
image.plot(tauGrid, tauGrid, -logLtau, xlab='tau1', ylab='tau2')
# image(tauGrid, tauGrid, -logLtau, xlab='tau1', ylab='tau2')
abline(v=times, h=times, col=8)
if(exportFig){dev.off()}

################################################################################
# Thining
################################################################################ 
seed <- 1; set.seed(seed)
n <- 15; v <- 2/3
times <- sort(runif(n))
train <- which(runif(n) < v); test <- (1:n)[-train]
if(exportFig){png(paste0(figDir, 'FigSegPP-ThiningOriginal.png'))}
plot(times, rep(0, n), pch=20, xlab='', ylab='', axes=0, cex=2)
abline(h=0)
if(exportFig){dev.off()}
if(exportFig){png(paste0(figDir, 'FigSegPP-ThiningSampling.png'))}
plot(times, rep(0, n), pch=20, xlab='', ylab='', axes=0, cex=2)
abline(h=0)
points(times[train], rep(0, length(train)), pch=20, cex=2, col='blue')
points(times[test], rep(0, length(test)), pch=20, cex=2, col='red')
if(exportFig){dev.off()}
