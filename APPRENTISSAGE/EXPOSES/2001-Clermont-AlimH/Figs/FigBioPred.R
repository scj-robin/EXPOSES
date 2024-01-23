# Fig classif

# Parms
n = 200; r0 = 7; seed = 1; set.seed(seed); par(pch=20)

# Data
r = c(rnorm(n), r0+rnorm(n))
a = 2*pi*runif(2*n)
x = r*cos(a); y = r*sin(a)

# Classif
c1 = kmeans(cbind(x, y), centers=rbind(c(-r0, r0), c(0 ,0)))$cluster
c2 = kmeans(r, centers=c(0, r0))$cluster

# par(mfrow=c(2, 2), col.axis=0, fg=0, pch=20)
pdf('Classif0.pdf'); 
par(mfrow=c(1, 1), col.axis=0, fg=0, pch=20, cex=1.5, mex=.5)
plot(x, y, col=1, xlab='', ylab='')
dev.off()

pdf('Classif1.pdf'); 
par(mfrow=c(1, 1), col.axis=0, fg=0, pch=20, cex=1.5, mex=.5)
plot(x, y, col=2*c1, xlab='', ylab='')
dev.off()

pdf('Classif2.pdf'); 
par(mfrow=c(1, 1), col.axis=0, fg=0, pch=20, cex=1.5, mex=.5) 
plot(x, y, col=2*c2, xlab='', ylab='')
dev.off()

pdf('Classif3.pdf'); 
par(mfrow=c(1, 1), col.axis=0, fg=0, pch=20, cex=1.5, mex=.5)
plot(r, a, col=2*c2, xlab='', ylab='')
dev.off()
