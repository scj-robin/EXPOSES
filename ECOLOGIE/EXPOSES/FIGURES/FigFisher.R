# FigFisher

library(SPECIES)

data(butterfly)
x = butterfly$j
Cx = butterfly$n_j
n = length(x)
x = x[-n]; Cx = Cx[-n]; n = n-1

# Data
plot(x, Sx, pch=20, cex=1.5, xlim=c(0, (max(x))), ylim = c(0, 150))

# Geometric fit
C = sum(Cx)
p = 1 - n / C
lines(x, n*dgeom((x-1), p), col='blue', lwd=2)

