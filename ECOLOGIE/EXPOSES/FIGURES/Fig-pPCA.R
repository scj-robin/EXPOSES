rm(list=ls())
par(mfrow=c(1, 1), pch=20)

# Example of PCA
require(graphics)
PCA = prcomp(USArrests, scale = TRUE)
print(PCA)
plot(PCA)
summary(PCA)
p = ncol(USArrests)
pdf('Fig-pPCA-exPCA-inertia.pdf')
plot(0:p, c(0, cumsum(PCA$sdev^2))/p, type='b', main='', ylab='', xlab='', pch=20)
dev.off()
pdf('Fig-pPCA-exPCA-biplot.pdf')
biplot(PCA); abline(h=0, v=0, lty=2)
dev.off()

# Structure of Sigma
n = 1e2; p = 20; q = 5
W = matrix(rnorm(n*q), n, q)
B = matrix(rnorm(p*q), p, q)
Sigma = B%*%t(B) + diag(p)
lambda = eigen(Sigma)$val
maxLambda = max(lambda)
nbCol = 10*round(max(lambda))
pdf('Fig-pPCA-Sigma.pdf')
image(1:p, 1:p, maxLambda - Sigma, xlab='', ylab='', col=heat.colors(nbCol))
dev.off()
pdf('Fig-pPCA-rotSigma.pdf')
image(1:p, 1:p, maxLambda-diag(lambda), xlab='', ylab='', col=heat.colors(nbCol))
dev.off()
# image(1:p, 1:p, max(Sigma) - Sigma, xlab='', ylab='')
Y = W %*% t(B) + matrix(rnorm(n*p), n, p)
PCA = prcomp(Y)
plot(PCA)
P = PCA$rotation
pdf('Fig-pPCA-empCov.pdf')
image(1:p, 1:p, max(cov(Y)) - cov(Y), xlab='', ylab='')
dev.off()
pdf('Fig-pPCA-rotCov.pdf')
image(1:p, 1:p, max(cov(Y%*%P)) - cov(Y%*%P), xlab='', ylab='')
dev.off()

