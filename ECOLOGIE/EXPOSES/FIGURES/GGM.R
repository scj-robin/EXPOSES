# GGM -> Plots
library(sna);library(mvtnorm)
seed = 1; set.seed(seed)

# Function
F_Mat2Tex <- function(A){
   r = nrow(A); c = ncol(A)
   sapply(1:(r-1), function(i){
      sapply(1:(c-1), function(j){cat(A[i, j], ' & ')})
      cat(A[i, c], '\\\\ \n')
   })
   sapply(1:(c-1), function(j){cat(A[r, j], ' & ')})
   cat(A[r, c], '\n')
}

# Parms
Omega = matrix(0, 4, 4)
Omega[upper.tri(Omega)] = c(.5, 0, -.3, 0, .2, .6)
Omega = Omega+t(Omega)
diag(Omega) = rep(1, 4)
Sigma = solve(Omega)

F_Mat2Tex(round(Sigma, 1))
F_Mat2Tex(round(cov2cor(Sigma), 1))
F_Mat2Tex(Omega)

# png('GM4nodes.png')
gplot(1*(Omega!=0), gmode='graph', label=1:4, label.pos=5, vertex.col=0, label.cex=2, vertex.cex=1.5)
# dev.off()

# Data
n = 50
Y = rmvnorm(n, sigma=Sigma)
Sigma.hat = cov(Y)
Omega.hat = solve(Sigma.hat)
F_Mat2Tex(round(Sigma.hat, 2))
F_Mat2Tex(round(Omega.hat, 2))
