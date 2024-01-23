# Figures for ECAS / SFdS / NAIM : reminder on prob. dist

rm(list=ls())
library(combinat); library(lattice); library(mvtnorm); library(mixtools); library(sna); library(latex2exp)
source('FigECAS-Functions.R')

# Joint discrete distribution
seed = 5; set.seed(seed)
x = hcube(c(2, 2, 2)) - 1; x = x[, (3:1)]
p = runif(nrow(x)); p = p / sum(p); p = round(p, 2)
Ftab(cbind(x, p))

# Joint continuous distribution
lattice.options(default.theme = standard.theme(color = FALSE))
mu = c(0, 0); Sigma = matrix(c(1, 1, 1, 2), 2, 2)
grid = seq(-3.5, 3.5, length.out=100)
x = mu[1] + sqrt(Sigma[1, 1])*grid; y = x #mu[2] + sqrt(Sigma[2, 2])*grid; 
z = matrix(0, length(x), length(y))
invisible(sapply(1:nrow(z), function(i){z[i, ] <<- dmvnorm(cbind(rep(x[i], length(y)), y), mu, Sigma)}))
pdf('FigECAS-BivariateNormal.pdf')
image(x, y, z, xlab=TeX('$x_1$'), ylab=TeX('$x_2$'), cex.lab=2)
ellipse(mu, Sigma, col=1, alpha=.1, lwd=3); 
ellipse(mu, Sigma, col=1, alpha=.5, lwd=3); 
ellipse(mu, Sigma, col=1, alpha=.9, lwd=3); 
dev.off()

# Conditional independence
n = 1e3; nodeLabels = c('X', 'Y', 'Z'); nodePos = matrix(c(-.77, 0, .77, 0, 0, -1), 3, 2, byrow=TRUE)
pdf('FigECAS-CondIndep-plot1.pdf'); par(font.lab=3)
mu1 = 1.5*c(-1, -1); mu2 = 1.5*c(1, 1); Sigma1 = Sigma2 = matrix(c(1.5, 0, 0, 1), 2, 2); 
XY1 = rmvnorm(n, mu1, Sigma1); XY2 = rmvnorm(n, mu2, Sigma2); 
Z = c(rep(2, n), rep(4, n))
# x = seq(min(XY1), max(XY1), length.out=1e2); y = seq(min(XY2), max(XY2), length.out=1e2)
# z = matrix(0, length(x), length(y))
# invisible(sapply(1:nrow(z), function(i){
#    z[i, ] <<- .5*dmvnorm(cbind(rep(x[i], length(y)), y), mu1, Sigma1) +
#       .5*dmvnorm(cbind(rep(x[i], length(y)), y), mu2, Sigma2)}))
# contour(x, y, z, xlab='x', ylab='y', cex.lab=3, levels=c(1e-3, 5e-3, 1e-2))
plot(rbind(XY1, XY2), col=Z, pch=20, cex=.2, xlab='x', ylab='y', cex.lab=3)
ellipse(mu1, Sigma1, col=2, alpha=.2, lwd=3); ellipse(mu2, Sigma2, col=4, alpha=.2, lwd=3); 
dev.off()
pdf('FigECAS-CondIndep-graph1.pdf'); par(mex=.1, font=3)
G1 = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3); 
gplot(G1, gmode='graph', coord=nodePos, label=nodeLabels, label.pos=5, 
      vertex.col=0, vertex.cex=4, label.cex=5)
dev.off()

pdf('FigECAS-CondIndep-plot2.pdf'); par(font.lab=3)
mu = c(0, 0); Sigma = matrix(c(1, .75, .75, 1), 2, 2);
XY = rmvnorm(2*n, mu, Sigma); Z = 2 + 2*(XY[, 1] > 0)
plot(XY, col=Z, xlab='x', ylab='y', pch=20, cex=.2, cex.lab=3, asp=1)
ellipse(mu, Sigma, col=1, alpha=.2, lwd=3)
dev.off()
pdf('FigECAS-CondIndep-graph2.pdf'); par(mex=.1, font=3)
G2 = matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), 3, 3); 
gplot(G2, gmode='graph', coord=nodePos, label=nodeLabels, label.pos=5, 
      vertex.col=0, vertex.cex=4, label.cex=5)
dev.off()


pdf('FigECAS-CondIndep-plot3.pdf'); par(font.lab=3)
mu1 = c(0, 0); mu2 = c(0, 0); 
Sigma1 = matrix(c(1, .75, .75, 1), 2, 2); Sigma2 = matrix(c(1, -.75, -.75, 1), 2, 2);
XY1 = rmvnorm(n, mu1, Sigma1); XY2 = rmvnorm(n, mu2, Sigma2); 
Z = c(rep(2, n), rep(4, n))
plot(rbind(XY1, XY2), col=Z, xlab='x', ylab='y', pch=20, cex=.2, cex.lab=3, asp=1)
ellipse(mu1, Sigma1, col=2, alpha=.2, lwd=3); ellipse(mu2, Sigma2, col=4, alpha=.2, lwd=3); 
dev.off()
pdf('FigECAS-CondIndep-graph3.pdf'); par(mex=.1, font=3)
G3 = matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3); 
gplot(G3, gmode='graph', coord=nodePos, label=nodeLabels, label.pos=5, 
      vertex.col=0, vertex.cex=4, label.cex=5)
dev.off()

pdf('FigECAS-CondIndep-plot4.pdf'); par(font.lab=3)
mu = c(0, 0); Sigma = matrix(c(1, 0, 0, 1), 2, 2)
XY = rmvnorm(2*n, mu, Sigma)
Z = 2+(1+sign(XY[, 1])*sign(XY[, 2]))
plot(XY, col=Z, xlab='x', ylab='y', pch=20, cex=.2, cex.lab=3, asp=1)
ellipse(mu, Sigma, col=1, alpha=.2, lwd=3)
dev.off()
pdf('FigECAS-CondIndep-graph4.pdf'); par(mex=.1, font=3)
G3 = matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3); 
gplot(G3, gmode='graph', coord=nodePos, label=nodeLabels, label.pos=5, 
      vertex.col=0, vertex.cex=4, label.cex=5)
dev.off()


