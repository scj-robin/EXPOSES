# Figures for continuous time SBM
rm(list=ls())
source('/home/robin/PgmR/Markov/FunctionsMarkov.R')
source('/home/robin/RECHERCHE/RESEAUX/DynamicNetwork/ctsbm/code/functions.r')

# Parms
# • the number of individuals – N .
# • the number of groups – K .
# • the expected number of interactions per pair per unit time – G.
# • the expected number of changes in group membership per node per unit time – R.
# • the scale between within- and across- group interaction rates – D
n = 5; K = 3; Tmax = 1; G = 4; R = .5; D = .75 # default in simulation study LMR18
n = 5; K = 3; Tmax = 1; G = 2.5; R = 3; D = .75 
seed = 1; set.seed(seed); lwd = 6;

# Functions
F_PlotHiddenPath <- function(times, nodepaths, lwd=4){
   # times = times; nodepaths = truenodepaths
   n = ncol(nodepaths); L = nrow(nodepaths); Tmax = max(times)
   plot(0, 0, col=0, main='', xlab='time', ylab='node', xlim=c(0, Tmax), ylim = c(.5, (n+.5)), 
        xaxt='n')   
   for (i in 1:n){
      path = nodepaths[, i]
      jumps = which(path[1:(L-1)]!=path[2:L])
      states = c(path[jumps], path[L])
      changes = c(0, times[which(path[1:(L-1)]!=path[2:L])], Tmax)
      invisible(sapply(1:(length(changes)-1), function(c){
         lines(c(changes[c], changes[c+1]), i*c(1, 1), col=1+states[c], lwd=lwd)
      }))
   }
}


# Simulation
simstudydata(n, K, Tmax, G, R, D, 'tmp.Rdata')
load('tmp.Rdata')
times = c(0, edges$edges$timestamp, Tmax)
M = length(times)-2

# Hidden path
pdf('FigCTSBM-pathZ.pdf')
par(mfrow=c(1, 1), mex=.75)
F_PlotHiddenPath(times, truenodepaths, lwd=lwd)
dev.off()

# Hidden + observed  path
pdf('FigCTSBM-pathZY.pdf')
par(mfrow=c(1, 1), mex=.75)
F_PlotHiddenPath(times, truenodepaths, lwd=lwd)
sapply(1:M, function(m){
   lines(edges$edges$timestamp[m]*c(1, 1), c(edges$edges$sender[m], edges$edges$receiver[m]), 
         lwd=2, type='b')
})
dev.off()

# Observed path
pdf('FigCTSBM-pathY.pdf')
par(mfrow=c(1, 1), mex=.75)
F_PlotHiddenPath(times, 7*(truenodepaths>0), lwd=lwd)
sapply(1:M, function(m){
   lines(edges$edges$timestamp[m]*c(1, 1), c(edges$edges$sender[m], edges$edges$receiver[m]), 
         lwd=2, type='b')
})
dev.off()
