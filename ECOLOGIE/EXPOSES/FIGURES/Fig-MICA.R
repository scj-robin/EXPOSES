# Figure cours MICA

Dir = "D:/RECHERCHE/METAGENOME/EXPOSES/FIGURES/"

##---------------------------------------------------------------------------------------
## Ancova  parms
##---------------------------------------------------------------------------------------
#t = (0:4)
#T = length(t)
#a = c(1, 1)
#b = c(1.5, 2.5)
#c = c(1.5, 1.5)
#d = c(0, -.75)
#sigma = .75
#color = c(2, 3)
#char = c(24, 25)
#R = 3
#
#---------------------------------------------------------------------------------------
## Ancova Gaussien
#---------------------------------------------------------------------------------------
#data = c()
#plot(t, 6*t, col=0, xlab="dose", ylab="response", ylim=c(0, 35), cex.lab=1.5)
#for (i in (1:2))
#  {
#  for (j in (1:2))
#    {
#      for ( r in (1:R))
#        {
#        Ytmp = a[i] + b[i]*t + c[i]*t^2 + d[j] + sigma*rnorm(T)
#        points(t, Ytmp, col=color[i], lwd=2, pch=char[j], cex=1.5)
#        data = rbind(data, cbind(rep(i, T), rep(j, T), t, Ytmp))
#        }
#    }
#  }
#data
##Response = data[,4]
#Strain = as.factor(data[,1])
#Prep = as.factor(data[,2])
#Time = as.factor(data[,3])
#lm(Response ~ Strain+Time+Strain:Time+Prep)
#write.table(data, file="MICA-Ancova.txt")

#---------------------------------------------------------------------------------------
## Ancova Poisson
#---------------------------------------------------------------------------------------
#data = c();
#plot(t, 6*t, col=0, xlab="time", ylab="response", ylim=c(0, 35), cex.lab=1.5)
#for (i in (1:2))
#  {
#  for (j in (1:2))
#    {
#      for ( r in (1:R))
#        {
#        Ytmp = rpois(T, exp((a[i] + b[i]*t + c[i]*t^2 + d[j])/10))
#        points(t, Ytmp, col=color[i], lwd=2, pch=char[j], cex=1.5)
#        data = rbind(data, cbind(rep(i, T), rep(j, T), t, Ytmp))
#        }
#    }
#  }
#data
#Response = data[,4]
#Strain = as.factor(data[,1])
#Prep = as.factor(data[,2])
#Time = as.factor(data[,3])
#glm(Response ~ Strain+Time+Strain:Time+Prep, family=poisson(link = "log"))
#write.table(data, file="MICA-AncovaPois.txt")

#---------------------------------------------------------------------------------------
## Link functions
#---------------------------------------------------------------------------------------
#mu = seq(0, 1, by=0.01)
#plot(mu, log(mu/(1-mu)), type="l", lwd=4, col=2, xlab = "mu", ylab = "log(mu/1-mu)", cex.lab=1.5)
#lines(mu, rep(0, length(mu)), col=1, lwd=4)
#lines(rep(0, length(mu)), log(mu/(1-mu)), col=1, lwd=4)
#lines(rep(1, length(mu)), log(mu/(1-mu)), col=1, lwd=4, lty=2)
#
#mu = seq(0, 10, by=0.1)
#plot(mu, log(mu), type="l", lwd=4, col=2, xlab = "mu", ylab = "log(mu)", cex.lab=1.5)
#lines(mu, rep(0, length(mu)), col=1, lwd=4)
#lines(rep(0, length(mu)), log(mu), col=1, lwd=4)

##---------------------------------------------------------------------------------------
## Effet aléatoire
##---------------------------------------------------------------------------------------
#m = 7
#mu = m+c(-4, 1, 4)
#s = .3
#g = 1
#J = 2
#R = 3
#
#data = c()
#M = (1/g/sqrt(2*acos(-1)))
#plot(0, 0, col=0, xlab="response", ylab="", xlim=c(min(mu)-2, max(mu+2)), ylim = c(0, 1.5*M), cex.lab=1.5)
#for (i in (1:length(mu)))
#  {
#  lines(c(mu[i], mu[i]), c(0, 1.5*M), col=1, lwd=2, lty=1)
#  x = mu[i] + g*3*seq(-1, 1, by=.01)
#  lines(x, dnorm(x, mean=mu[i], sd=g), col=2, lwd=2, lty=2)
#  A = g*rt(J, 5)
#  for (j in (1:J))
#    {
#    lines(c(mu[i]+A[j], mu[i]+A[j]), c(0, M), col=2, lwd=2, lty=1)
#    x = mu[i] + A[j] + s*3*seq(-1, 1, by=.01)
#    lines(x, .4*s/g*dnorm(x, mean=mu[i]+A[j], sd=s), , col=4, lwd=2, lty=2)
#    Y = s*rnorm(R)
#    points(mu[i]+A[j]+Y, rep(0, R), pch='+', col=4, cex=2)
#    data = rbind(data, cbind(rep(i, R), rep(j, R), mu[i]+A[j]+Y))
#    }
#  }
#data
#Response = data[,3]
#Trt = data[,1]
#Indiv = data[,2]
#write.table(data, file="MICA-AnovaMixed.txt")

##---------------------------------------------------------------------------------------
## Effet aléatoire Poisson
##---------------------------------------------------------------------------------------
#m = 4.5
#mu = (m+c(-2, 0.5, 1))/1.5
#g = 10
#J = 2
#R = 3
#Ymax = 60
#
#data = c()
#M = 0.22
##plot(0, 0, col=0, xlab="response", ylab="", xlim=c(min(mu)-2, max(mu+2)), ylim = c(0, 1.5*M), cex.lab=1.5)
#plot(0, 0, col=0, xlab="response", ylab="", xlim=c(0, Ymax), ylim=c(0, 1.2*M), cex.lab=1.5)
#for (i in (1:length(mu)))
#  {
#  lines(exp(c(mu[i], mu[i])), c(0, 1.5*M), col=1, lwd=2, lty=1)
#  x = (0:Ymax)
#  #lines(x, dpois(x, exp(mu[i])), col=2, lwd=2, lty=2)
#  lines(x, 5*M*dgamma(x, shape=g, scale=exp(mu[i])/g), col=2, lwd=2, lty=2)
#  A = rgamma(J, shape=g, scale=1/g)
#  for (j in (1:J))
#    {
#    lines(exp(mu[i])*A[j]*c(1, 1), c(0, .75*M), col=2, lwd=2, lty=1)
#    #x = mu[i] + A[j] + s*3*seq(-1, 1, by=.01)
#    lines(x, M*dpois(x, exp(mu[i])*A[j]), col=4, lwd=2, lty=2)
#    Y = rpois(R, exp(mu[i])*A[j])
#    points(Y, rep(0, R), pch='+', col=4, cex=2)
#    data = rbind(data, cbind(rep(i, R), rep(j, R), Y))
#    }
#  }
#Response = data[,3]
#Trt = data[,1]
#Indiv = data[,2]
#write.table(data, file="MICA-AnovaMixedPois.txt")
#
##---------------------------------------------------------------------------------------
## Poisson / Binomiale negative
##---------------------------------------------------------------------------------------
#m = 4.5
#n = 10000
#g = 3
#X = rpois(n, m)
#G = rgamma(n, shape=g, scale=1/g)
#Y = rpois(n, m*G)
#par(mfrow=c(1, 2))
#hist(X, breaks=ceiling(sqrt(n)), col=4, xlim=c(0, 30), ylim=c(0, 1),
#  main="", xlab="", ylab="", freq=F)
#hist(Y, breaks=ceiling(sqrt(n)), col=4, xlim=c(0, 30), ylim=c(0, 1),
#  main="", xlab="", ylab="", freq=F)
#
##---------------------------------------------------------------------------------------
## FDR
##---------------------------------------------------------------------------------------
#m = 10000
#pi1 = 0.1
#lambda = 0.75
#rho = 0.95
#M = rep(0, m)-rbinom(m, 1, pi1)*rexp(m, lambda)
#X = M+rnorm(m)
#P = pnorm(X)
#par(mfrow=c(1, 2))
##hist(X, breaks=ceiling(sqrt(m)), main="", xlab="", ylab="", col=4)
#hist(P, breaks=ceiling(sqrt(m)), main="", xlab="", ylab="", col=4, freq=F)
#lines(c(0, 1), (1-pi1)*c(1, 1), lwd=3, col=2, lty=2)
#Y = rnorm(sqrt(m))
##Y = M + rho*rep(Y, sqrt(m)) + sqrt(1-rho^2)*rnorm(m)
#Q = pnorm(Y)
##hist(Y, breaks=ceiling(sqrt(m)), main="", xlab="", ylab="", col=4)
#hist(Q, breaks=ceiling(sqrt(m)), main="", xlab="", ylab="", col=4, freq=F)
#lines(c(0, 1), (1-pi1)*c(1, 1), lwd=3, col=2, lty=2)
#

S = 1e4
lambda = 5
sigma = c(0, .2, 5)
X = rnorm(S)
par(mfrow=c(1, length(sigma)))
for (s in sigma)
  {
  Y = rpois(S, lambda*exp(s*X))
  print(c(mean(Y), sd(Y), sd(Y)/mean(Y)))
  if (sd(Y)/mean(Y) <= 10)
    {
    hist(Y, breaks=ceiling(sqrt(S)), col=4, main="", xlab="", ylab="", freq=F)
    }
  if (sd(Y)/mean(Y) > 10)
    {
    hist(log10(Y), breaks=ceiling(sqrt(S)), col=4, main="", xlab="", ylab="", freq=F)
    }
  }
