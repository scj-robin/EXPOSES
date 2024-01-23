# Fig function phi for W-graph

# Parms
a = 3
b = 1
n = 50
x = seq(0, 1, by=1/(n-1))

# Phi
f = dbeta(x, a, b)
plot(x, f, type="l")
phi = matrix(0, n, n)
for (i in (1:n))
{
  for  (j in (1:n))
  {
    phi[i, j] = f[i]*f[j]
  }
}

# SBM
K = 5
alpha = (sqrt(a))^(-(1:K))
alpha = alpha / sum(alpha)
alcum = c(0, cumsum(alpha))
pi = matrix(0, K, K)
SBM = matrix(0, n, n)
for (k in (1:K))
{
  ik = which((x - alcum[k])*(x - alcum[k+1]) <= 0)
  for (l in (1:K))
  {
    jl = which((x - alcum[l])*(x - alcum[l+1]) <= 0)
    pi[k, l] = mean(phi[ik, jl])
    SBM[ik, jl] = pi[k, l]
  }
}
  
# Plot
theta = -5
#par(mfrow=c(2, 2))
postscript('FigWgraph-image.ps')
image(phi)
dev.off()
postscript('FigWgraph-persp.ps')
persp(phi, xlab='', ylab='', zlab='', border='black', theta=theta)
dev.off()
postscript(paste('FigWdiscret-K', K, '-image.ps', sep=''))
image(SBM)
dev.off()
postscript(paste('FigWdiscret-K', K, '-persp.ps', sep=''))
persp(SBM, xlab='', ylab='', zlab='', border='black', theta=theta)
dev.off()
