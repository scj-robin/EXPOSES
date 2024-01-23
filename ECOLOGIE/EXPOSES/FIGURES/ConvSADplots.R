# plots for convex SAD estimation

rm(list=ls())

# Parms for geometric example
p = .2
xmax = 5

# Convex distribution
x = (0:xmax)
g = dgeom(x, p)
png(paste('GeomTri-p', round(10*p), '.png', sep=''))
plot(x, g, lwd=4, col='black', pch=20, cex=1.5, 
     ylim=c(0, max(g)), xlim=c(0, xmax), xlab='', ylab='')

# Triangular decomposition : Tj(i) = 2(j-i)/j/(j+1)
d2g = diff(diff(g))
j.list = (1:length(d2g))
prop = j.list*(j.list+1)*d2g/2
gt = rep(0, xmax)
for (j in j.list){
  lines(c(0, j), prop[j]*c(2/(j+1), 0), lty=2, lwd=2, col='blue')
  text(0, 2*prop[j]/(j+1)+0.02, labels='T', font=3, cex=2)
  text(0.12, 2*prop[j]/(j+1), labels=j, cex=1.5)
  gt[1:(j+1)] = gt[1:(j+1)] + prop[j] * 2/j/(j+1)*(j:0)
}
lines(x, g, lwd=3)
dev.off()

# # Convex SAD distribution
# g = dgeom(x, p)
# gSAD = g
# gSAD[1] = gSAD[1] - prop[1]
# 
# #png(paste('GeomTri-p', round(10*p), '.png', sep=''))
# plot(x, g, lwd=4, col='black', pch=20, cex=1.5, 
#      ylim=c(0, max(g)), xlim=c(0, xmax), xlab='', ylab='')
# 
# # Triangular decomposition : Tj(i) = 2(j-i)/j/(j+1)
# gt = rep(0, xmax)
# for (j in j.list){
#   lines(c(0, j), prop[j]*c(2/(j+1), 0), lty=2, lwd=2, col='blue')
#   text(0, 2*prop[j]/(j+1)+0.02, labels='T', font=3, cex=2)
#   text(0.12, 2*prop[j]/(j+1), labels=j, cex=1.5)
#   gt[1:(j+1)] = gt[1:(j+1)] + prop[j] * 2/j/(j+1)*(j:0)
# }
# lines(x, g, lwd=3)
# #dev.off()
