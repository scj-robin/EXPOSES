a1 = 1; b1 = 1
a2 = 6; b2 = 4
n = 20; y = 4
p.vec = seq(0, 1, by=.001)

pdf('Fig-BetaBinom-Prior1.pdf')
plot(p.vec, dbeta(p.vec, a1, b1), lwd=4, lty=2, col=4, type='l', ylim=c(0, 5), main='', xlab='', ylab='')
dev.off()

pdf('Fig-BetaBinom-Prior1-Data.pdf')
plot(p.vec, dbeta(p.vec, a1, b1), lwd=4, lty=2, col=4, type='l', ylim=c(0, 5), main='', xlab='', ylab='')
abline(v=y/n, lwd=4, col=1)
dev.off()

pdf('Fig-BetaBinom-Prior1-Data-Posterior1.pdf')
plot(p.vec, dbeta(p.vec, a1, b1), lwd=4, lty=2, col=4, type='l', ylim=c(0, 5), main='', xlab='', ylab='')
abline(v=y/n, lwd=4, col=1)
lines(p.vec, dbeta(p.vec, a1+y, b1+n-y), lwd=4, col=4)
dev.off()

pdf('Fig-BetaBinom-Prior1-Data-ConfInter1.pdf')
plot(p.vec, dbeta(p.vec, a1, b1), lwd=4, lty=2, col=4, type='l', ylim=c(0, 5), main='', xlab='', ylab='')
abline(v=y/n, lwd=4, col=1)
lines(p.vec, dbeta(p.vec, a1+y, b1+n-y), lwd=4, col=4)
abline(v=qbeta(c(.025, .975), a1+y, b1+n-y), lwd=4, col=4, lty=3)
dev.off()

pdf('Fig-BetaBinom-Prior2.pdf')
plot(p.vec, dbeta(p.vec, a1, b1), lwd=4, lty=2, col=4, type='l', ylim=c(0, 5), main='', xlab='', ylab='')
abline(v=y/n, lwd=4, col=1)
lines(p.vec, dbeta(p.vec, a1+y, b1+n-y), lwd=4, col=4)
lines(p.vec, dbeta(p.vec, a2, b2), lwd=4, lty=2, col=2)
dev.off()

pdf('Fig-BetaBinom-Prior2-Data-Posterior2.pdf')
plot(p.vec, dbeta(p.vec, a1, b1), lwd=4, lty=2, col=4, type='l', ylim=c(0, 5), main='', xlab='', ylab='')
abline(v=y/n, lwd=4, col=1)
lines(p.vec, dbeta(p.vec, a1+y, b1+n-y), lwd=4, col=4)
lines(p.vec, dbeta(p.vec, a2, b2), lwd=4, lty=2, col=2)
lines(p.vec, dbeta(p.vec, a2+y, b2+n-y), lwd=4, col=2)
dev.off()

pdf('Fig-BetaBinom-Prior2-Data-ConfInter2.pdf')
plot(p.vec, dbeta(p.vec, a1, b1), lwd=4, lty=2, col=0, type='l', ylim=c(0, 5), main='', xlab='', ylab='')
abline(v=y/n, lwd=4, col=1)
lines(p.vec, dbeta(p.vec, a1+y, b1+n-y), lwd=4, col=4)
abline(v=qbeta(c(.025, .975), a1+y, b1+n-y), lwd=4, col=4, lty=3)
# lines(p.vec, dbeta(p.vec, a2, b2), lwd=4, lty=2, col=2)
lines(p.vec, dbeta(p.vec, a2+y, b2+n-y), lwd=4, col=2)
abline(v=qbeta(c(.025, .975), a2+y, b2+n-y), lwd=4, col=2, lty=3)
dev.off()

