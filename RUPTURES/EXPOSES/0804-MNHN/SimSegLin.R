rm(list=ls())
source("D:/RECHERCHE/RUPTURES/Exposes/0804-MNHN/Dynprog.R")

# Data simulation
l = c(5, 5, 5)
K = length(l)
n = sum(l)
x = (1:n)
y = vector()
for (k in (1:K))
  {
  xtmp = (1:l[k])
  a = rnorm(1)
  b = rnorm(1)/10
  ytmp = a + b*xtmp + rnorm(l[k])/5
  y = c(y, ytmp)
  }
plot(x, y)
rupt = cumsum(l)
for (k in (1:(K-1)))
  { 
  abline(v=rupt[k]+0.5, col=2, lwd=2)
  }


# Cost matrix
lmin = 2
C = matrix(Inf, n, n)
for (i in (1:(n-lmin)))
  {
  for (j in ((i+lmin):n))
    {
    reg = lm(y[i:j] ~ x[i:j])
    C[i, j] = sum(reg$residuals^2)
    }
  }
C

# Dynamic programmation
DP = DynProg(C, n)
DP$t.est
plot(DP$J.est)