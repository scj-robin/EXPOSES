# Figures for the PLN model

# Parms
mu <- 1
sigma <- 2
B <- 1e8

# Simul
Z <- mu+ sigma*rnorm(B)
Y <- rpois(B, exp(Z))

# Conditional
pdf(paste0('FigPLN-pZcondY-mu', mu, '-sigma', sigma, '.pdf'))
yList <- c(3, 2, 1, 0); yNb <- length(yList); 
zGrid <- seq(min(Z), max(Z), length.out=1000)
par(mfrow=c(1, 1))
for(i in 1:yNb){
   Zy <- Z[which(Y==yList[i])]
   dens <- density(Zy)
   if(i==1){
      shift <- .05*max(dens$y)
      plot(dens, col=1+i, main='', xlab='', ylab='', lwd=3, ylim=c(0, 1.2*max(dens$y)), cex.axis=1)
   }else{
      dens <- density(Zy)
      lines(dens, col=1+i, lwd=3)
   }
   text(dens$x[which.max(dens$y)], max(dens$y)+shift, paste0('Y=', yList[i]), col=1+i, cex=2)
   lines(zGrid, dnorm(zGrid, mean=mean(Zy), sd=sd(Zy)), lwd=3, lty=3, col=1+i)
}
lines(zGrid, dnorm(zGrid, mean=mu, sd=sigma), lwd=2, col=1, lty=2)
dev.off()
