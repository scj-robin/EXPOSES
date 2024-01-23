# Compute the joint posterior classification probability for each pair (u,v) of the grid and each pair (k,l) of classes

jointPosterior <- function(U, V, alpha) {

  K<-length(alpha) # number of classes
  N<-length(U) # size of the grid : N*N


  # Simulation
  #S = 1000000
  S=1000
  X = matrix(0, S, K)
  for (j in (1:K))
    {X[,j] = rgamma(S, shape=alpha[j])
   }
  X = X/rowSums(X)
  Xcum = t(apply(X, 1, cumsum))
  jpost<-array(0, dim=c(N, N, K*(K+1)/2))
  
  for(u in 1:N) {
    for(v in u:N){
      #Cuvkl<-matrix(0, K, K)
      Ruvkl<-matrix(0, K, K)
      for(k in 1:K) {
        kboundLeft = Xcum[, k-1]
        kboundRight = Xcum[, k]
        if(k==1) {
          kboundLeft = 0
        }
        else if(k==K) {
          kboundRight = 1
        }
        for(l in k:K) {
          #Cuvkl[k, l]<- F_P2Classif(U[u], V[v], k, l, alpha)
          lboundLeft = Xcum[, l-1]
          lboundRight = Xcum[, l]
          if(l==1){
            lboundLeft = 0
          }
          else if(l==K) {
            lboundRight = 1
          }
          
          Ruvkl[k, l] = sum((kboundLeft<=U[u])*(U[u]<=kboundRight)*(lboundLeft<=V[v])*(V[v]<=lboundRight))/S
         #   cat(U[u], V[v], k, l, Cuvkl[k, l],"\n")
         # cat(U[u], V[v], k, l, Ruvkl[k, l],"\n")
        }
      }
      #jpost[u, v, ]<-Cuvkl[upper.tri(Cuvkl, diag=TRUE)]
      jpost[u, v, ]<-Ruvkl[upper.tri(Ruvkl, diag=TRUE)]
     
      }
  }
  return(jpost)
}
