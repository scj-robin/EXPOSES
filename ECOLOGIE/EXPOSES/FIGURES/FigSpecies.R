# FigSpecies

rm(list=ls())

library(SPECIES)
data.list = c('butterfly', 'cottontail', 'insects', 'microbial', 'EST')

data(butterfly); 
n = dim(butterfly)[1]
png('butterfly.png')
plot(butterfly$j[-n], butterfly$n_j[-n], pch=20, cex=1.5, xlab='', ylab='')
dev.off()

data(cottontail); 
png('cottontail.png')
plot(cottontail$j, cottontail$n_j, pch=20, cex=1.5, xlab='', ylab='')
dev.off()

data(insects); 
n = dim(insects)[1]
png('insects.png')
plot(insects$j[-n], insects$n_j[-n], pch=20, cex=1.5, xlab='', ylab='')
dev.off()

data(microbial); 
png('microbial.png')
plot(microbial$X1, microbial$X2, pch=20, cex=1.5, xlab='', ylab='')
dev.off()

data(EST)
n = dim(EST)[1]
png('EST.png')
plot(EST[-n, 1], EST[-n, 2], pch=20, cex=1.5, xlab='', ylab='')
dev.off()