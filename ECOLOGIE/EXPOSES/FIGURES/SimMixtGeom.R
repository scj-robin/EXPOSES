# Param
K = 3;  n = 500
pi = (K:1);   pi = pi/sum(pi)
gamma = .6^(1:K)

# Simul
N = rmultinom(1, n, pi)
X = c()
for (k in (1:K))
  {X = c(X, rgeom(N[k], gamma[k]))}
H = hist(X, nclass=max(X)+1)

# Mixture
xx = H$breaks[-length(H$breaks)]
g = rep(0, length(xx))
for (k in (1:K))
  {g = g + pi[k]*dgeom(xx, gamma[k])}

# Plot
np = sum((X > 0))
#postscript("SimMixtGeom.ps")
plot(xx, np*g, col=0, main="", xlab="", ylab="")
points(xx[-1], H$counts[-1])
lines(xx[-1], np*g[-1], col=2, lwd=2)
lines(xx[1:2], np*g[1:2], col=4, lwd=2, lty=2)
#dev.off()
dev.copy2pdf(file="SimMixtGeom.pdf")
