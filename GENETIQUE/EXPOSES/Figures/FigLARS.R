library(lars)
data(diabetes)
attach(diabetes)
Res = lars(x2,y,trace=TRUE,max.steps=80)
pmax = 30
gamma = -log(c(0, Res$lambda))

postscript("Fig-LARS-coef.ps")
plot(gamma[1:(1+pmax)], Res$beta[1:(1+pmax),1], col=1, type='l', lwd=2,
  ylim=c(min(Res$beta[1:(1+pmax),]), max(Res$beta[1:(1+pmax),])),
  xlab='-log(lambda)', ylab='theta', cex.lab=2)
for (j in (2:(dim(Res$beta)[2]-1)))
  {
  lines(gamma[1:(1+pmax)], Res$beta[1:(1+pmax),j], col=j, lwd=2)
  }
abline(h=0, lty=2, col=1)
dev.off()

postscript("Fig-LARS-CV.ps")
par(cex.lab=2)
cv.lars(x2,y,trace=TRUE,max.steps=80)
dev.off()

