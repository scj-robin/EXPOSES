# Figures for Poisson point process segmentation

rm(list=ls()); palette('R3')
# source('/home/robin/RECHERCHE/RUPTURES/Hawkes/segHP/codes/PoissonCV-R/FunctionSegPP-CV.R')
# source('/home/robin/RECHERCHE/RUPTURES/Hawkes/segHP/codes/Poisson-R/DynProg.R')
# source('/home/robin/RECHERCHE/RUPTURES/Hawkes/segHP/codes/Poisson-R/Gsegmentation.R')
# source('/home/robin/RECHERCHE/RUPTURES/Hawkes/segHP/codes/Poisson-R/G.contrast.R')
source('/home/robin/RECHERCHE/RUPTURES/Hawkes/segHP/codes/Functions/FunctionCrossValidation.R')
source('/home/robin/RECHERCHE/RUPTURES/Hawkes/segHP/codes/Functions/FunctionsContrasts.R')
source('/home/robin/RECHERCHE/RUPTURES/Hawkes/segHP/codes/Functions/FunctionsSegmentation.R')
source('/home/robin/RECHERCHE/RUPTURES/Hawkes/segHP/codes/Functions/FunctionsSimulation.R')
source('/home/robin/RECHERCHE/RUPTURES/Hawkes/segHP/codes/Functions/FunctionValidationCriteria.R')

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
# breaks <- TimesProc(times)
# cost <- CostPP(breaks$N.Potential.Times, breaks$Potential.Times, a=a, b=b) # Cost matrix
breaks <- PotentialBreaks(times)
cost <- CostMatrix(breaks=breaks, ContrastSeg=ContrastSeg.PG.ab, parms=list(a=1, b=1/length(times), lga=0))
dp <- DynProg(cost, Kmax)
seg <- BestSeg.K(breaks, dp, Kopt)
seg$lambda <- seg$DN / seg$Dt

# Plots
if(exportFig){png(paste0(figDir, figName, '-', dataName, '-Path.png'))}
plot(c(0, times, 1), c((0:n), n), type='s', xlab='t', ylab='N(t)')
if(exportFig){dev.off()}

if(exportFig){png(paste0(figDir, figName, '-', dataName, '-Seg.png'))}
PlotSeg(times=times, seg=seg)
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
dataName <- paste0('simul-n', n, '-K1-seed', seed)

# Simul
if(exportFig){png(paste0(figDir, 'FigSegPP-', dataName, '-Path.png'))}
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
image.plot(tauGrid, tauGrid, -logLtau, xlab=expression(tau[1]), ylab=expression(tau[2]))
# image(tauGrid, tauGrid, -logLtau, xlab='tau1', ylab='tau2')
abline(v=times, h=times, col=8)
if(exportFig){dev.off()}

################################################################################
# Thining
################################################################################ 
seed <- 1; set.seed(seed); exportFig <- TRUE
tau <- c(0, .4, .55, .825, .925, 1)
lambda <- 10*c(4, 1, 3, 1.5, 5)
times <- SimulHPP(lambda=lambda, tau=tau)
n <- length(times)
v <- 4/5
train <- sort(sample(1:n, size=round(v*n))); test <- (1:n)[-train]

par(mfrow=c(2, 1))
if(exportFig){png(paste0(figDir, 'FigSegPP-ThiningOriginal.png'))}
plot(times, rep(0, n), pch=20, xlab='', ylab='', axes=0, cex=2)
abline(h=0)
if(exportFig){dev.off()}

if(exportFig){png(paste0(figDir, 'FigSegPP-ThiningOriginalLambda.png'))}
plot(0, 0, xlab='', ylab='', axes=0, col=0, xlim=c(0, 1), ylim=c(0, max(lambda)))
abline(h=0, v=0)
for(k in 1:length(lambda)){lines(tau[k:(k+1)], rep(lambda[k], 2), lwd=4)}
abline(v=tau, lty=2, col=8)
if(exportFig){dev.off()}

if(exportFig){png(paste0(figDir, 'FigSegPP-ThiningSampling.png'))}
plot(times, rep(0, n), pch=20, xlab='', ylab='', axes=0, cex=2)
abline(h=0)
points(times[train], rep(0, length(train)), pch=20, cex=2, col='blue')
points(times[test], rep(0, length(test)), pch=20, cex=2, col='red')
if(exportFig){dev.off()}

if(exportFig){png(paste0(figDir, 'FigSegPP-ThiningSamplingLambda.png'))}
plot(0, 0, xlab='', ylab='', axes=0, col=0, xlim=c(0, 1), ylim=c(0, max(lambda)))
abline(h=0, v=0)
for(k in 1:length(lambda)){lines(tau[k:(k+1)], rep(lambda[k], 2), lwd=4, lty=3)}
abline(v=tau, lty=2, col=8)
for(k in 1:length(lambda)){lines(tau[k:(k+1)], v*rep(lambda[k], 2), lwd=4, col=4)}
for(k in 1:length(lambda)){lines(tau[k:(k+1)], (1-v)*rep(lambda[k], 2), lwd=4, col=2)}
if(exportFig){dev.off()}

