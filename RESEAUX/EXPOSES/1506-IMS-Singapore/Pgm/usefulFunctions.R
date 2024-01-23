# UsefulFunction.R

surfSBM<-function(grid, alpha, PI) {
  
  s<-dim(grid)[1]
  U<-seq(0, 1, length=s)
  V<-seq(0, 1, length=s)
  K<-length(alpha)
  Calpha<-c(0, cumsum(alpha[-K]))
  phi<-matrix(0, s, s)
  
  for(u in 1:s) {
    tU<-max(which(Calpha <= U[u]))
    for(v in u:s) {
      tV<-max(which(Calpha <= V[v]))
      phi[u, v]<-PI[tU, tV]
      phi[v, u]<-phi[u,v]
    }
  }
  return(phi)
}

surf<-function(grid, beta, lambda) {
  
  return( beta*(grid)^lambda )
  
}

surf2<-function(grid, beta, lambda) {
  s<-dim(grid)[1]
  U<-seq(0, 1, length=s)
  V<-seq(0, 1, length=s)
  phi<-matrix(0, s, s)
  for(u in 1:s) {
    for(v in u:s) {
      phi[u, v]<-beta/((1+exp(-lambda*(U[u]-1/2)))*(1+exp(-lambda*(V[v]-1/2))))
      phi[v, u]<-phi[u, v]
    } 
  }
  return(phi)
}

surf3<-function(grid, beta) {
  s<-dim(grid)[1]
  U<-seq(0, 1, length=s)
  V<-seq(0, 1, length=s)
  phi<-matrix(0, s, s)
  for(u in 1:s) {
    for(v in u:s) {
      phi[u, v]<-beta*sin(U[u])^2*sin(V[v])^2
      phi[v, u]<-phi[u, v]
    } 
  }
  return(phi)
  
}




