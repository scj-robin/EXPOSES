# Trace d'un nouvelle version du réseau
# des opérons d'E. coli

# Dir
DataDir = '/home/robin/Bureau/RECHERCHE/RESEAUX/VB-EM/S-Gazal/EColi oriente/'
Y = read.table(paste(DataDir, 'data2', sep=''))
Z = read.table(paste(DataDir, 'attributeVEM', sep=''), skip=1)$V3

# Groups
n = nrow(Y); K = max(Z)
N = sapply((1:K), function(k){sum(Z==k)})
Z = order(N, decreasing=T)[Z]
N = N[order(N, decreasing=T)]

# Number nodes within groups
C = rep(0, K)
Nik = rep(0, n)
invisible(sapply(1:n, function(i){
   C[Z[i]] <<- C[Z[i]]+1
   Nik[i] <<- C[Z[i]]
}))

# Coordinates
W = (N-1)^(1/5) #log(N)
e = 1 / max(W)/3 # rapport hauteur/largeur de l'ellipse
X = matrix(rnorm(2*n), n, 2)
invisible(sapply(1:n, function(i){
   theta = 2*pi*Nik[i]/N[Z[i]]
   X[i, ] <<- W[Z[i]]*c(cos(theta), e*sin(theta)) + c(0, Z[i])
}))

# Plot
png('im_EcoliVEM_2_hierachical.png')
par(xaxt='n', yaxt='n', bty='n')
plot(0, 0, xlim = c(min(X[, 1]), max(X[, 1])), 
     ylim = c(min(X[, 2]), max(X[, 2])), col=0, xlab='', ylab='')
invisible(sapply(1:n, function(i){
   sapply(1:n, function(j){
      if ((Z[i] < Z[j]) & Y[i, j]==1){
         lines(X[c(i, j), 1], X[c(i, j), 2], col=1, lwd=2)
         }})}))
invisible(sapply(1:n, function(i){
   sapply(1:n, function(j){
      if ((Z[i] == Z[j]) & Y[i, j]==1){
         lines(X[c(i, j), 1], X[c(i, j), 2], col=2, lwd=2)
         }})}))
invisible(sapply(1:n, function(i){
   sapply(1:n, function(j){
      if ((Z[i] > Z[j]) & Y[i, j]==1){
         lines(X[c(i, j), 1], X[c(i, j), 2], col=4, lwd=2)
         }})}))
points(X[, 1], X[, 2], pch=20, col=8, cex=.5)
dev.off()