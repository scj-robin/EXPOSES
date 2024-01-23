

x = seq(-4, 4, by=.01)
log.g = -log(1+exp(-x))
plot(x, log.g, type='l')

f = -log(exp(x/2) + exp(-x/2))
plot(x, f, type='l')
xi = .75
lower.f = - xi/2 - log(1+exp(-xi)) - 1/(4*xi)*tanh(xi/2)*(x^2 - xi^2)
lines(x, lower.f, col=4)

pdf('ApproxJaJ00.pdf')
plot(x, log.g, type='l', lwd=4, xlab='', ylab='')
lines(x, x/2+lower.f, col=4, lty=2, lwd=4)
abline(v=xi, lty=2, lwd=4)
dev.off()
