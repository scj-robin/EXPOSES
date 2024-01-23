# Figures for ECAS / SFdS / NAIM : function

# Latex display for a vector
Fvec <- function(a, linebreak=T){
   sapply(1:(length(a)-1), function(i){cat(a[i], '& ')})
   if(linebreak){cat(a[length(a)], '\\\\ \n')}else{cat(a[length(a)], '\n')}
}

# Latex display for a table
Ftab <- function(A){
   sapply(1:(nrow(A)-1), function(i){Fvec(A[i, ])})
   Fvec(A[nrow(A), ], linebreak=F)
}

# Create Sigma and Omega from G
MakeOmegaSigma <- function(G, sign=TRUE){
   p = ncol(G) 
   Gsign = F_Vec2Sym(F_Sym2Vec(G * matrix(2*rbinom(p^2, 1, .5)-1, p, p)))
   if(sign==FALSE){Gsign = G}
   lambda = 1.1; Omega = lambda*diag(rowSums(G)) - Gsign
   while(min(eigen(Omega)$values) < 1e-6){lambda = 1.1*lambda; Omega = lambda*diag(rowSums(G)) - G}
   sigma = 2*rgamma(p, 1, 1); 
   corSigma = cov2cor(solve(Omega)); Sigma = diag(sigma)%*%corSigma%*%diag(sigma)
   Omega = solve(Sigma); Omega = Omega * (abs(Omega) > 1e-10)
   return(list(Omega=Omega, Sigma=Sigma))
}

