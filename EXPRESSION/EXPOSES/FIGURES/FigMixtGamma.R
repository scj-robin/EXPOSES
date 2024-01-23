# Plots of gamma mixture

# parms
a = 3^(1:4)
K = length(a)
p = rep(1/K, K)
x = seq(0, 3, by=.01)

# plot
png('FigMixtGamma.png')
plot(x, x, ylim=c(0, 1), col=0, xlab='', ylab='')
g = rep(0, length(x))
for (k in (1:K)){
  fk = dgamma(x, a[k], a[k])
  g = g + p[k] * fk
  lines(x, p[k]*fk, col=k, lwd=3)
}
#lines(x, g, col=1, lwd=2)
dev.off()