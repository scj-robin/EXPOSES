# Transformation for the negative binomial

# Parms
n = 1e4
p = .1*c(1, .01)
phi = 50*c(1, 2)

# Simulation
par(mfcol=c(3, 2))
X = rnbinom(n, phi[1], p[1])
hist(X, breaks = ceiling(sqrt(n)), main=paste('Counts: sd=', round(sd(X), 3)), ylab='', , xlab='')
L = log(1+X)
hist(A, breaks = ceiling(sqrt(n)), main=paste('log: sd=', round(sd(L), 3)), ylab='', , xlab='')
A = asinh(X/phi)
hist(L, breaks = ceiling(sqrt(n)), main=paste('arcsinh: sd=', round(sd(A), 3)), ylab='', , xlab='')

X = rnbinom(n, phi[2], p[1])
hist(X, breaks = ceiling(sqrt(n)), main=paste('Counts: sd=', round(sd(X), 3)), ylab='', , xlab='')
L = log(1+X)
hist(A, breaks = ceiling(sqrt(n)), main=paste('log: sd=', round(sd(L), 3)), ylab='', , xlab='')
A = asinh(X/phi)
hist(L, breaks = ceiling(sqrt(n)), main=paste('arcsinh: sd=', round(sd(A), 3)), ylab='', , xlab='')
