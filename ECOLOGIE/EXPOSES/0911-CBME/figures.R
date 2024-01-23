# Figure 1 : distribution

y = hist(rpois(300,2.9),breaks=seq(0,20)-0.5,plot=FALSE)
X11(width=6,height=4)
par(mar=c(2.5,2.5,1.5,0.5)+0.1,mgp=c(1.5,0.5,0),cex=1.,font.axis=2,font.lab=2)
plot(y$mids[-1],y$counts[-1], main="Frequency/Count data", xlab="rare species -------> abundant species", ylab="nb of observations")
dev.copy2pdf(file="figure1.pdf")

# A truncated mixture of f is equivalent to a mixture of truncated f (Bohning 2006)
em_kpoissont = function (data,k) { #{{{
    # EM fit mixture of k truncated Poisson 
    # starting values : classify range(data) into k groups
    hist0 = hist(data,breaks=seq(min(data)-0.5,max(data)+0.5,length.out=k+1),plot=FALSE)
    pi0 = hist0$counts / length(data)
    theta0 = hist0$mids
    Q0 = -Inf;
    data_hist = hist(data,breaks=seq(min(data)-0.5,max(data)+0.5,by=1),plot=FALSE)
    data_support = data_hist$mids[which(data_hist$counts > 0)];
    data_counts = data_hist$counts[which(data_hist$counts > 0)];
    mdata = max(data)

    # Helper function
    f = function(x,y) {return(x / (1-exp(-x))-y)}
    
    repeat {
    # E step calcul des tau
        tau = matrix(0,k,length(data_support)) # k rows, tau(x) in line
        for (j in seq.int(1,k)) {
            tau[j,] = pi0[j] * dpois(data_support,theta0[j]) }
        tau = tau / rep(colSums(tau),each=k)
    # M step
        pi1 = pmax(tau %*% matrix(data_counts) / sum(data_counts),0.001)
        # moyenne des tau avec counts
        temp_theta =  tau %*% matrix(data_counts*data_support) / tau %*% matrix(data_counts) 
        theta1 = temp_theta;
        # print(temp_theta);
        for (j in seq.int(1,k)) { # solve for theta
        theta1[j] = uniroot(f,c(0.0001,mdata),max(temp_theta[j],1.001))$root}
        
    # Update
        if (max(c(abs(pi1-pi0),abs(theta1-theta0))) < 0.000001) { break }
        pi0 = pi1; theta0 = theta1;
        #print(c(pi0,theta0))
    }
    return(c(pi0,theta0))
}#}}}
rkpoissont = function(n,pi0,lambda0) {#{{{
    # simulate mixture of k poisson with weights pi0
    k = length(pi0);
    classes = sample(seq(1,k),n,replace=TRUE,prob=pi0);
    temp = rpois(n,lambda0[classes])
    return(temp[which(temp>0)])
}#}}}

toto = rkpoissont(1000,c(0.25,0.25,0.25,0.25), c(5,10,20,25))
em_kpoissont(toto,4)

resb = matrix(0,nrow=100,ncol=8)
for (i in seq.int(1,100)) {
    toto = rkpoissont(4000,c(0.1,0.2,0.3,0.4), c(3,10,20,30))
    resb[i,] = em_kpoissont(toto,4)
    print(c(i,resb[i,]))
}
summary(resb)

resb = matrix(0,nrow=10000,ncol=8)
for (i in seq.int(1,10000)) {
    toto = rkpoissont(4000,c(0.22,0.24,0.26,0.28), c(3,10,20,30))
    resb[i,] = em_kpoissont(toto,4)
    print(c(i,resb[i,]))
}
summary(resb)




  
correctem = function (params) {
    # params = c(theta,lambda)
    n = length(params)/2
    P_0 = matrix(0,n); q_j = matrix(0,n)
    for (j in seq(1,n)) { P_0[j] = dpois(0,params[n+j]) }
    for (j in seq(1,n)) { q_j[j] = params[j] / (1-P_0[j]) }
    return(c(q_j/sum(q_j)))
}

toto = rkpoissont(1000,c(0.2,0.24,0.26,0.3), c(5,10,20,25))
params0 = em_kpoissont(toto,4)
params0
correctem(params0)

nb_samples = 2000
resb = matrix(0,nrow=nb_samples,ncol=15)
for (i in seq.int(1,nb_samples)) {
    toto = rkpoissont(2000,c(0.1,0.2,0.3,0.4), c(1,10,20,30))
    resb[i,seq.int(1,8)] = em_kpoissont(toto,4)
    resb[i,seq.int(1,4)+8] = correctem(resb[i,seq.int(1,8)])
    resb[i,13] = length(toto)
    resb[i,14] = length(toto) / (1-sum(resb[i,seq.int(1,4)] * dpois(0,resb[i,seq.int(5,8)])))
    resb[i,15] = length(toto) / (1-sum(resb[i,seq.int(9,12)] * dpois(0,resb[i,seq.int(5,8)])))
    print(i)
}
summary(resb)
save(resb,file="em_kpoissont_corrected.Rdata")

load("em_kpoissont_corrected.Rdata")

resb=resb[which(resb[,1]>0.002),]

dp1 = density(resb[,1])
dp2 = density(resb[,2])
dp3 = density(resb[,3])
dp4 = density(resb[,4])

X11(width=5,height=3)
par(mar=c(2.5,0.5,1.5,0.5)+0.1,mgp=c(1.5,0.5,0),cex=1.,font.axis=2,font.lab=2)
plot(NULL, xlim=c(0,0.5),ylim=c(0,70), yaxt="n", main="Proportions", xlab="", ylab="")

points(dp1$x,dp1$y,col=1)
points(dp2$x,dp2$y,col=2)
points(dp3$x,dp3$y,col=3)
points(dp4$x,dp4$y,col=4)
abline(v=c(0.1,0.2,0.3,0.4),col=seq(1,4))
text(c(0.1,0.2,0.3,0.4),rep.int(70,4),
   c(expression(pi[1]==0.1),expression(pi[2]==0.2),expression(pi[3]==0.3),expression(pi[4]==0.4)),col=seq(1,4))
dev.copy2pdf(file="figure2.pdf")


dp1 = density(resb[,5])
dp2 = density(resb[,6])
dp3 = density(resb[,7])
dp4 = density(resb[,8])

X11(width=5,height=3)
par(mar=c(2.5,0.5,1.5,0.5)+0.1,mgp=c(1.5,0.5,0),cex=1.,font.axis=2,font.lab=2)
plot(NULL, xlim=c(0,31),ylim=c(0,3), yaxt="n", main="Parameters", xlab="", ylab="")

points(dp1$x,dp1$y,col=1)
points(dp2$x,dp2$y,col=2)
points(dp3$x,dp3$y,col=3)
points(dp4$x,dp4$y,col=4)
abline(v=c(1,10,20,30),col=seq(1,4))
text(c(1,10,20,30),rep.int(3,4),
   c(expression(theta[1]==1),expression(theta[2]==10),expression(theta[3]==20),expression(theta[4]==30)),col=seq(1,4))
dev.copy2pdf(file="figure3.pdf")


dp1 = density(resb[,9])
dp2 = density(resb[,10])
dp3 = density(resb[,11])
dp4 = density(resb[,12])

X11(width=5,height=3)
par(mar=c(2.5,0.5,1.5,0.5)+0.1,mgp=c(1.5,0.5,0),cex=1.,font.axis=2,font.lab=2)
plot(NULL, xlim=c(0,0.5),ylim=c(0,70), yaxt="n", main="Corrected Proportions", xlab="", ylab="")

points(dp1$x,dp1$y,col=1)
points(dp2$x,dp2$y,col=2)
points(dp3$x,dp3$y,col=3)
points(dp4$x,dp4$y,col=4)
abline(v=c(0.1,0.2,0.3,0.4),col=seq(1,4))
text(c(0.1,0.2,0.3,0.4),rep.int(70,4),
   c(expression(pi[1]==0.1),expression(pi[2]==0.2),expression(pi[3]==0.3),expression(pi[4]==0.4)),col=seq(1,4))
dev.copy2pdf(file="figure4.pdf")


dp1 = density(resb[,13])
dp2 = density(resb[,14])
dp3 = density(resb[,15])

X11(width=5,height=3)
par(mar=c(2.5,0.5,1.5,0.5)+0.1,mgp=c(1.5,0.5,0),cex=1.,font.axis=2,font.lab=2)
plot(NULL, xlim=c(1900,2100),ylim=c(0,0.05), yaxt="n", main="Number of species", xlab="", ylab="")

points(dp1$x,dp1$y,col=1)
points(dp2$x,dp2$y,col=2)
points(dp3$x,dp3$y,col=3)
abline(v=2000,col=3)
text(c(1930,1980,2030),rep.int(0.05,3),
   c("Observed","Estimated","Corrected"),col=seq(1,3))
dev.copy2pdf(file="figure5.pdf")


