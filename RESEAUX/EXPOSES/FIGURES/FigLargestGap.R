# Figures for the largest gap algorithm

rm(list=ls())
source('/home/robin/PgmR/General/FunctionsMatVec.R')
seed <- 1; set.seed(seed)
exportFig <- TRUE

# Dims
nList <- round(10^(seq(2, 4, length.out=3))); nNb <- length(nList)

# Parms
K <- 4
powerProp <- 3/2; probPower <- 1
prop <- (1:K)^powerProp; prop <- sort(prop/sum(prop))
dens <- .1
prob <- (1:K)^probPower%o%(1:K)^probPower
prob <- prob * (dens / (prop%*%prob%*%prop)[1, 1])
probBar <- as.vector(prob%*%prop)

# Sim
simFile <- paste0('FigLargestGap-Degree-seed', seed, '.Rdata')
if(file.exists(simFile)){
   load(simFile)
}else{
   degreeList <- list()
   for (i in 1:nNb){ # i <- 1
      n <- nList[i]
      Z <- t(rmultinom(n, 1, prop))
      Y <- F_Vec2Sym(F_Sym2Vec(matrix(rbinom(n^2, 1, Z%*%prob%*%t(Z)), n, n), diag=FALSE), diag=FALSE)
      degreeList[[i]] <- rowSums(Y)
   }
   save(degreeList, file=simFile)
}

# Plot
if(exportFig){pdf('FigLargestGap-Histograms.pdf')}
par(mfrow=c(nNb, 1), mex=.6)
for(i in 1:nNb){
   n <- nList[i]
   hist(degreeList[[i]]/(n-1), breaks=sqrt(max(nList)), xlab='', ylab='', main=paste('n =', nList[i]))
   abline(v = probBar, col=4, lwd=3)
}
if(exportFig){dev.off()}
