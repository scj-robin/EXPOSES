# Figure pour U-stats & B-EDD

rm(list=ls())
library(gtools); library(sna); library(bmotif); library(polyaAeppli)

# Dirs
dirFig <- './'

# Figs
mex <- 1; cex <- 10; lwdEdge <- 4; lwdCurve <- 8; cex.main <- 8; mex <- .6
exportFig <- TRUE
figName <- paste0(dirFig, 'FigUstatsBEDD')

###############################################################################
# Modele B-EDD
###############################################################################
# Dims and Parms
m <- 30; n <- 50; lambda <- 10; mu.f <- 2.5; mu.g <- 4
seed <- 0; set.seed(seed)
F_Gmu <- function(u, mu){mu*u^(mu-1)}
F_Outline <- function(I, J, col=1){
   II <- c(min(I)-.5, max(I)+.5)
   JJ <- c(min(J)-.5, max(J)+.5)
   lines(rep(II[1], 2), JJ, col=col, lwd=2)
   lines(rep(II[2], 2), JJ, col=col, lwd=2)
   lines(II, rep(JJ[1], 2), col=col, lwd=2)
   lines(II, rep(JJ[2], 2), col=col, lwd=2)
}

# Data
U <- runif(m); V <- runif(n); Lambda <- lambda*F_Gmu(U, mu.f)%o%F_Gmu(V, mu.g)
Y <- matrix(rpois(m*n, Lambda), m, n)
if(exportFig){pdf(paste0(figName, '-mat-f', round(10*mu.f), '-g', round(10*mu.g), '.pdf')); par(mfrow=c(1, 1), mex=mex)}
image(1:m, 1:n, Y, xlab='', ylab='', asp=1, xlim=c(-2, m), ylim=c(-2, n), axes=0)
text(-1, n/2, 'plants', srt=90)
text(m/2, -1, 'insects')
if(exportFig){dev.off()}

# U-stat
if(exportFig){pdf(paste0(figName, '-ustat-f', round(10*mu.f), '-g', round(10*mu.g), '.pdf')); par(mfrow=c(1, 1), mex=mex)}
image(1:m, 1:n, Y, xlab='', ylab='', asp=1, xlim=c(-2, m), ylim=c(-2, n), axes=0)
text(-1, n/2, 'plants', srt=90)
text(m/2, -1, 'insects')
F_Outline(c(5, 6), c(15, 16), col=4)
F_Outline(c(25, 26), c(10, 10), col=3)
F_Outline(c(25, 26), c(17, 17), col=3)
F_Outline(c(10, 10), c(30, 30), col=1)
F_Outline(c(20, 20), c(30, 30), col=1)
F_Outline(c(10, 10), c(40, 40), col=1)
F_Outline(c(20, 20), c(40, 40), col=1)
if(exportFig){dev.off()}
