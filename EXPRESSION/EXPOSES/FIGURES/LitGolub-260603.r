Golub.data<-read.table("e:\\cours\\bioinfo\\golub\\train-7070.prn",header=T)
MA <- apply(Golub.data[,2:28],1,mean)
MB <- apply(Golub.data[,29:39],1,mean)
sigma<