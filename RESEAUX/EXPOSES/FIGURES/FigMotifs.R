# Motifs figures

library(igraph)
rm(list=ls())

E = list(); Motif = list(); 
Motif[[1]] = 'I'; E[[1]] = matrix(c(0, 1, 0, 0), 2, 2); 
Motif[[2]] = 'V'; E[[2]] = matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0), 3, 3); 
Motif[[3]] = 'Triangle'; E[[3]] = matrix(c(0, 1, 1, 0, 0, 1, 0, 0, 0), 3, 3); 
Motif[[4]] = 'Square'; E[[4]] = matrix(c(0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0), 4, 4); 
Motif[[5]] = 'Chain4'; E[[5]] = matrix(c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0), 4, 4); 
Motif[[6]] = 'Clique4'; E[[6]] = matrix(c(0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0), 4, 4); 
Motif[[7]] = 'Whisker'; E[[7]] = matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 4, 4); 
Motif[[8]] = 'Star3'; E[[8]] = matrix(c(0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 4, 4); 
Motif[[9]] = 'Star4'; E[[9]] = matrix(c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 5, 5); 
Motif[[10]] = 'Star5'; E[[10]] = matrix(c(0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 6, 6); 
Motif[[11]] = 'SquareDiag'; E[[11]] = matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0), 4, 4); 
M = length(E)

# Symetrisation
k = list()
lapply(1:M, function(m){
   E[[m]] <<- E[[m]] + t(E[[m]])
   k[[m]] <<- dim(E[[m]])[1]
   })
lapply(1:M, function(m){
   cat(Motif[[m]], k[[m]], '\n')
   print(E[[m]])
})

# Coordonnees
Coord = list()
lapply(1:M, function(m){
   if (substr(Motif[[m]], 1, 4) == 'Star'){
      Theta = (-0.5 + (1:(k[[m]]-1))) * 2 * acos(-1) / (k[[m]]-1)
      C = rbind(c(0, 0), cbind(cos(Theta), sin(Theta)));
   }else{
      Theta = (-0.5 + (1:k[[m]])) * 2 * acos(-1) / k[[m]]
      C = cbind(cos(Theta), sin(Theta));
   }
   Coord[[m]] <<- C
})

# Plots
lapply(1:M, function(m){
   FigFile = paste('FigMotif-', Motif[[m]], '.png', sep='')
   png(FigFile)
   par(plt=c(0, 1, 0, 1))
   plot(Coord[[m]][, 1], Coord[[m]][, 2], pch=20, cex=15, xlab='', ylab='', main='', axes=F, 
        xlim = c(-1.1, 1.1), ylim=c(-1.1, 1.1))
   for (i in (1:(k[[m]]-1))){
      for (j in ((i+1):k[[m]])){
         if (E[[m]][i, j] > 0){
            lines(c(Coord[[m]][i, 1], Coord[[m]][j, 1]), c(Coord[[m]][i, 2], Coord[[m]][j, 2]), 
                  lwd=15, col=1)
         }
      }
   }
   dev.off()
   plot(graph.adjacency(E[[m]]), main=Motif[[m]])
})
