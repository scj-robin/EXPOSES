# Figures for significant recurrent alterations

rm(list=ls())
source('../../MinRegionCont/Pgm/F_SimBirthDeath.R')
par(mex=.3)

# #############################################################
# # Small L, small N: embedded MC
# #############################################################
# N = 20 # nb profiles
# M = ceiling(N/2+sqrt(N)) # threshold
# l = 3 # alteration length
# 
# # Transition matrix
# B = M # nb states below threshold
# A = N-M+1 # nb states below threshold
# P.size = B + l*A + 1
# P = matrix(0, P.size, P.size)
# Col = c(0, 8, 5, 3, 6, 2)
# Breaks = -.5+(0:length(Col))
# 
# png('Fig-MinReg-EmbMC-Porg.png')
# P[(1:(B+A)), (1:(B+A))] = 1
# image(P, x=(0:(P.size-1)), y=(0:(P.size-1)), xlab='', ylab='', col=Col, breaks=Breaks)
# dev.off()
# 
# png('Fig-MinReg-EmbMC-Thres.png')
# P[(1:B), ((B+1):(B+A))] = 2
# P[((B+1):(B+A)), ((B+1):(B+A))] = 3
# P[((B+1):(B+A)), (1:B)] = 4
# image(P, x=(0:(P.size-1)), y=(0:(P.size-1)), xlab='', ylab='', col=Col, breaks=Breaks)
# abline(v = B-.5, h = B-.5, col=1, lwd=3, lty=2); 
# dev.off()
# 
# png('Fig-MinReg-EmbMC-1st.png')
# P[((B+1):(B+A)), ((B+1):(B+A))] = 0
# P[((B+A+1):(B+2*A)), ((B+1):(B+A))] = 3
# image(P, x=(0:(P.size-1)), y=(0:(P.size-1)), xlab='', ylab='', col=Col, breaks=Breaks)
# abline(v = B-.5, h = B-.5, col=1, lwd=3, lty=2); 
# dev.off()
# 
# png('Fig-MinReg-EmbMC-2nd.png')
# P[(1:B), ((B+A+1):(B+2*A))] = 2
# P[((B+2*A+1):(B+3*A)), ((B+A+1):(B+2*A))] = 3
# image(P, x=(0:(P.size-1)), y=(0:(P.size-1)), xlab='', ylab='', col=Col, breaks=Breaks)
# abline(v = B-.5, h = B-.5, col=1, lwd=3, lty=2); 
# lines(-.5+c(0, B), -.5+c(B+A, B+A))
# dev.off()
# 
# png('Fig-MinReg-EmbMC-3rd.png')
# P[(1:B), ((B+2*A+1):(B+3*A))] = 2
# P[((B+3*A+1):(B+3*A+1)), ((B+2*A+1):(B+3*A))] = 3
# image(P, x=(0:(P.size-1)), y=(0:(P.size-1)), xlab='', ylab='', col=Col, breaks=Breaks)
# abline(v = B-.5, h = B-.5, col=1, lwd=3, lty=2); 
# lines(-.5+c(0, B), -.5+c(B+A, B+A))
# lines(-.5+c(0, B), -.5+c(B+2*A, B+2*A))
# lines(-.5+c(0, M), .5+c(N, N))
# dev.off()
# 
# png('Fig-MinReg-EmbMC-Abs.png')
# P[(B+3*A+1), (B+3*A+1)] = 5
# image(P, x=(0:(P.size-1)), y=(0:(P.size-1)), xlab='', ylab='', col=Col, breaks=Breaks)
# abline(v = B-.5, h = B-.5, col=1, lwd=3, lty=2); 
# lines(-.5+c(0, B), -.5+c(B+A, B+A))
# lines(-.5+c(0, B), -.5+c(B+2*A, B+2*A))
# lines(-.5+c(0, M), .5+c(N, N))
# dev.off()
# 
# 
# #############################################################
# # Example for min regions
# #############################################################
# # Data
# DataDir = '../../ETUDES/Curie/MinRegion/Data1/'
# ChromName = paste(DataDir, 'chromlimit_curie070606.txt', sep='')
# BoolName = paste(DataDir, 'curiebool0.8-0.5.txt', sep='')
# #BoolName = paste(DataDir, 'curiebool0.8-0.5-NaN-gainloss.txt', sep='')
# 
# # Chromosomes
# Chrom = read.table(ChromName)$V1
# L = length(Chrom)
# NbChrom = max(Chrom)
# Limits = c(-1, which(abs(diff(Chrom))>0))
# ChromCenter = .5*(Limits[2:(NbChrom+1)]+Limits[1:NbChrom])
# 
# # States
# Bool = as.matrix(sign(read.table(BoolName)))
# p = dim(Bool)[1]
# #image(t(Bool), col=c(2, 0, 3))
# Loss = matrix(0, dim(Bool)[1], dim(Bool)[2])
# Loss[which(Bool < 0)] = 1; Loss = as.matrix(Loss)
# png('Fig-MinReg-Data-Prof.png')
# image(t(Loss), col=c(0, 2), y=(0:p), x=(0:L), xlab='', ylab='')
# abline(v=.5+Limits, lty=2, lwd=2, col=4)
# abline(h=-.5+(1:p), lty=1, lwd=1, col=8)
# text(ChromCenter, p/2, (1:NbChrom))
# dev.off()
# 
# # Cumulated profiles
# CumProf = colSums(Loss)
# png('Fig-MinReg-Data-Cumul.png')
# plot(1:L, CumProf, type='s', lwd = 2, xlab='', ylab='', ylim=c(0, p))
# abline(v=.5+Limits, lty=2, lwd=2, col=4)
# text(ChromCenter, 0, (1:NbChrom))
# dev.off()
# 
# # Profiles for a specific chromosome
# ChrNum = 13
# Select = which(Chrom==ChrNum)
# png(paste('Fig-MinReg-Data-Prof', ChrNum, '.png', sep=''))
# image(t(Loss[, Select]), col=c(0, 2), y=(0:p), x=((Limits[ChrNum]):Limits[ChrNum+1]), xlab='', ylab='')
# abline(h=-.5+(1:p), lty=1, lwd=1, col=8)
# dev.off()
# 
# # Cumulated profiles for a specific chromosome
# #M = 31; In = 1189; Out = 1193
# png(paste('Fig-MinReg-Data-Cumul', ChrNum, '.png', sep=''))
# plot(Select, CumProf[Select], type='s', lwd = 2, xlab='', ylab='', ylim=c(0, p))
# dev.off()
# 
# # Cumulated profiles for a specific chromosome + region
# M = 32; In = 1709; Out = 1714
# png(paste('Fig-MinReg-Data-Reg1-', ChrNum, '.png', sep=''))
# plot(Select, CumProf[Select], type='s', lwd = 2, xlab='', ylab='', ylim=c(0, p))
# abline(h=M, col=4, lty=3, lwd=2)
# abline(v=c(In, Out), lty=2, lwd=2, col=8)
# lines(c(In, Out), c(M, M), lwd=6, col=2)
# text((min(Select)), (M+2), 'm', font=3, cex=2)
# text((In+Out)/2, (p-1), 'l', font=3, cex=2)
# dev.off()
# 
# M = 62; In = 1687; Out = 1688
# png(paste('Fig-MinReg-Data-Reg2-', ChrNum, '.png', sep=''))
# plot(Select, CumProf[Select], type='s', lwd = 2, xlab='', ylab='', ylim=c(0, p))
# abline(h=M, col=4, lty=3, lwd=2)
# abline(v=c(In, Out), lty=2, lwd=2, col=8)
# lines(c(In, Out), c(M, M), lwd=6, col=2)
# text((min(Select)), (M+2), 'm', font=3, cex=2)
# text((In+Out)/2, (p-1), 'l', font=3, cex=2)
# dev.off()

# #############################################################
# # Large L, large N: diffusion
# #############################################################
# p = 10000 # nb profiles
# N = 10000 # nb events in the simulation
# l = .05 # alteration length
# lambda = .5; mu = 2.5
# tau = lambda+mu; theta = lambda/tau; sigma = sqrt(theta*(1-theta))
# Lambda = (p:0)*lambda; Mu = (0:p)*mu
# set.seed(10)
# BD = F_SimBirthDeath(round(N*lambda/(lambda+mu)), Lambda, Mu, N)
# M = .3 * p*lambda/(lambda+mu) + .7*max(BD$Count) # threshold
# Mtilde = (M - p*theta)/sqrt(p)/sigma
# n = sum(BD$Time <= 1)
# BD$Time = BD$Time[1:n]
# BD$Count = BD$Count[1:n]
# Z = (BD$Count - p*theta)/sqrt(p)/sigma
# # plot(BD$Time, BD$Count, type='s', lwd = 2, xlab='', ylab='')
# # plot(BD$Time, Z, type='s', lwd = 2, xlab='', ylab='')
# 
# # Convergence toward an Ornstein-Uhlenbeck process : V2
# PsampList = 10^(1:4)
# Tmax = BD$Time[N]
# for (p.samp in PsampList){
#    png(paste('Fig-MinReg-Diffus-Cvg', p.samp, '.png', sep=''))
#    L.samp = (p.samp:0)*lambda; M.samp = (0:p.samp)*mu
#    set.seed(10)
#    BD = F_SimBirthDeath(round(p.samp*lambda/(lambda+mu)), L.samp, M.samp, 10000)
#    Z = (BD$Count - p.samp*theta)/sqrt(p.samp)/sigma
#    plot(BD$Time, Z, type='s', lwd = 2, xlab='', ylab='', 
#         xlim=c(0, 1), ylim=c(-2, 2))   
#    dev.off()
# }
# # Convergence toward an Ornstein-Uhlenbeck process : V1
# SimLengthList = 10^(1:4)
# Tmax = BD$Time[N]
# for (SimLength in SimLengthList){
#    png(paste('Fig-MinReg-Diffus-Cvg', SimLength, '.png', sep=''))
#    plot(BD$Time[1:SimLength], BD$Count[1:SimLength], type='s', lwd = 2, 
#         xlab='', ylab='', xlim=c(0, Tmax*SimLength/N))   
#    dev.off()
# }
# 
# # Cumulated profile
# png('Fig-MinReg-Diffus-Thres.png')
# plot(BD$Time, Z, type='s', lwd = 2, xlab='', ylab='')
# abline(h=Mtilde, lwd=3, lty=3, col=4)
# dev.off()
# 
# # Excursions
# In = c(); Out = c()
# for (i in (2:n)){
#    if ((Z[i] >= Mtilde) & (Z[i-1] < Mtilde)){In = c(In, i)}
#    if ((Z[i] < Mtilde) & (Z[i-1] >= Mtilde)){Out = c(Out, i)}
# }
# 
# # Longest excursion 
# Lg = BD$Time[Out] - BD$Time[In]
# MaxLg = which.max(Lg)
# png('Fig-MinReg-Diffus-Longest.png')
# plot(BD$Time, Z, type='s', lwd = 2, xlab='', ylab='', xlim=c(0, 1))
# abline(h=Mtilde, lwd=3, lty=3, col=4)
# abline(v=BD$Time[In[MaxLg]]+c(0, l), lty=2, lwd=2, col=8)
# if (Lg[MaxLg] < l){
#    lines(c(BD$Time[In[MaxLg]], BD$Time[Out[MaxLg]]), c(Mtilde, Mtilde), lwd=6, col=3)
# }else{
#    lines(c(BD$Time[In[MaxLg]], BD$Time[Out[MaxLg]]), c(Mtilde, Mtilde), lwd=6, col=2)
# }
# dev.off()
# 
# # Brownian & O-U paths
# set.seed(10)
# Nsim = 5e2
# NbSim = 5
# ZB = matrix(0, NbSim, Nsim); ZOU = ZB
# for (i in (1:NbSim)){
#    ZOU[i, ] = arima.sim(n=Nsim, list(ar=.9, ma=0))
#    ZB[i, ] = cumsum(rnorm(Nsim))
# }
# 
# png('Fig-MinReg-Diffus-Bpaths.png')
# plot((1:Nsim)/Nsim, ZB[1,], type='l', col=0, 
#      ylim = c(min(ZB), max(ZB)), xlab='', ylab='', main='')
# for (i in (1:NbSim)){
#    lines((1:Nsim)/Nsim, ZB[i, ], col=1)
# }
# dev.off()
# 
# png('Fig-MinReg-Diffus-B-OUpaths.png')
# plot((1:Nsim)/Nsim, ZB[1,], type='l', col=0, 
#      ylim = c(min(ZB), max(ZB)), xlab='', ylab='', main='')
# for (i in (1:NbSim)){
#    lines((1:Nsim)/Nsim, ZB[i, ], col=1)
#    lines((1:Nsim)/Nsim, ZOU[i, ], col=2)
# }
# dev.off()

# #######################################n######################
# # Continuity issue
# #############################################################
# 
# set.seed(10)
# Nsim = 1e3; xx = (1:Nsim)/Nsim; 
# Z = arima.sim(n=Nsim, list(ar=.9, ma=0))
# plot(xx, Z, type='l', lwd=2, col=1, xlab='', ylab='', main='')
# abline(h = 3.2, col=4, lwd=2, lty=2)
# x.lim = c(.8, .95)
# abline(v = x.lim, col=8, lwd=2)
# diff = .19*exp(-abs(xx-.88)^2*500)
# 
# png('Fig-MinReg-Diffus-Continuity.png')
# plot(xx, Z, type='l', lwd=2, col=1, xlim=x.lim, , xlab='', ylab='', main='')
# lines(xx, (1-diff)*Z, lwd=2, col=2)
# lines(xx, Z, lwd=2, col=1)
# abline(h = 3.1, col=4, lwd=2, lty=2)
# abline(v = c(.868, .90), col=1, lwd=2, lty=2)
# abline(v = c(.8685, .885), col=2, lwd=2, lty=2)
# dev.off()

# #############################################################
# # Large L, large N: importance sampling
# #############################################################
# # Densities
# x.max = 20
# x.scale = 10*x.max
# xx = seq(0, x.max, length.out=1e3)
# df = 3
# f = dnorm(xx, m=df, sd=sqrt(2*df)) + dnorm(-xx, m=df, sd=sqrt(2*df))
# g = dchisq(xx, df)
# 
# # Sampling + weights
# B = 1e3
# X = rchisq(B, df)
# W = (dnorm(X, m=df, sd=sqrt(2*df)) + dnorm(-X, m=df, sd=sqrt(2*df))) / dchisq(X, df)
# 
# # Histograms
# H = hist(X, breaks=sqrt(B), main='', xlab='', ylab='', lwd=2)
# h = mean(diff(H$breaks))
# Hg = rep(0, length(H$breaks))
# for (b in (1:length(H$breaks))){Hg[b] = sum(X < H$breaks[b])}
# Hg = Hg / B
# Hg = diff(c(Hg, 1))
# Hf = rep(0, length(H$breaks))
# for (b in (1:length(H$breaks))){Hf[b] = sum(W[X < H$breaks[b]])}
# Hf = Hf / sum(W)
# Hf = diff(c(Hf, 1))
# 
# # Plots
# png('Fig-MinReg-Diffus-ImpSamp0.png')
# plot(xx/x.scale, g*x.scale, col=1, lty=2, lwd=3, type='l', main='', xlab='', ylab='', xaxt='n', , yaxt='n')
# lines(xx/x.scale, f*x.scale, col=2, lty=2, lwd=3, main='', xlab='', ylab='')
# # lines(xx/x.scale, 10*f/g, col=4, lwd=3, main='', xlab='', ylab='')
# dev.off()
# png('Fig-MinReg-Diffus-ImpSamp1.png')
# plot(H$breaks/x.scale, Hg*x.scale, col=1, lwd=3, type='s', main='', xlab='', ylab='', xaxt='n', , yaxt='n')
# dev.off()
# png('Fig-MinReg-Diffus-ImpSamp2.png')
# plot(H$breaks/x.scale, Hg*x.scale, col=1, lwd=3, type='s', main='', xlab='', ylab='', xaxt='n', , yaxt='n')
# lines(xx/x.scale, h*g*x.scale, lwd=3, lty=2)
# dev.off()
# png('Fig-MinReg-Diffus-ImpSamp3.png')
# plot(H$breaks/x.scale, Hg*x.scale, col=1, lwd=3, type='s', main='', xlab='', ylab='', xaxt='n', , yaxt='n')
# lines(H$breaks/x.scale, Hf*x.scale, type='s', lwd=3, col=2)
# lines(xx/x.scale, 10*f/g, col=4, lwd=3, main='', xlab='', ylab='')
# dev.off()
# png('Fig-MinReg-Diffus-ImpSamp4.png')
# plot(H$breaks/x.scale, Hg*x.scale, col=0, lwd=3, type='s', main='', xlab='', ylab='', xaxt='n', , yaxt='n')
# lines(H$breaks/x.scale, Hf*x.scale, type='s', lwd=3, col=2)
# lines(xx/x.scale, h*f*x.scale, lwd=3, lty=2, col=2)
# abline(v = .04, lwd=3, col=3)
# dev.off()

#############################################################
# Large L, small N: birth & death
#############################################################
N = 20 # nb profiles
M = ceiling(sqrt(N)) # threshold
l = .05 # alteration length
lambda = .1; mu = .5
Lambda = (N:0)*lambda; Mu = (0:N)*mu
seed = 53 # 53, 39, 20, 18, 3, 13
# seed = seed+1
set.seed(seed) 
Plot = T

BD = F_SimBirthDeath(round(N*lambda/(lambda+mu)), Lambda, Mu, 100)
while(max(BD$Count) <= M)
{BD = F_SimBirthDeath(round(N*lambda/(lambda+mu)), Lambda, Mu, 100)}
n = length(BD$Count)
BD$Time = BD$Time / BD$Time[n]

# Truncation
BD$Time = 3*BD$Time
BD$Count = BD$Count[which(BD$Time < 1)]
BD$Time = BD$Time[which(BD$Time < 1)]
n = length(BD$Count)

# BD$Time = 1-BD$Time; #BD$Count = BD$Count[n:1]
FigName = 'Fig-MinReg-BirDea-Thres.png'
if(Plot){png(FigName)}
plot(BD$Time, BD$Count, type='s', lwd = 2, xlab='', ylab='')
abline(h=M, lwd=3, lty=3, col=4)
if(Plot){dev.off()}

# Excursions
In = c(); Out = c()
for (i in (2:n)){
   if ((BD$Count[i] >= M) & (BD$Count[i-1] < M)){In = c(In, i)}
   if ((BD$Count[i] < M) & (BD$Count[i-1] >= M)){Out = c(Out, i)}
}
NbExt = length(In)

# Test for valid excursions
e = 0; test = F
while (test==F){
   e = e+1
   FigName = paste('Fig-MinReg-BirDea-In', e, '.png', sep='')
   if(Plot){png(FigName)}
   plot(BD$Time, BD$Count, type='s', lwd = 2, xlab='', ylab='')
   abline(h=M, lwd=3, lty=3, col=4)
   abline(v=BD$Time[In[e]], lty=2, lwd=2, col=1)
   if(Plot){dev.off()}
   
   FigName = paste('Fig-MinReg-BirDea-Time', e, '.png', sep='')
   if(Plot){png(FigName)}
   plot(BD$Time, BD$Count, type='s', lwd = 2, xlab='', ylab='')
   abline(h=M, lwd=3, lty=3, col=4)
   abline(v=BD$Time[In[e]], lty=2, lwd=2, col=1)
   if (e > 1){
      lines(c(BD$Time[In[e-1]], BD$Time[Out[e-1]]), c(M, M), lwd=6, col=3)
      abline(v=c(BD$Time[In[e-1]], BD$Time[Out[e-1]]), lty=2, lwd=2, col=8)
   }
   if(Plot){dev.off()}
   
   FigName = paste('Fig-MinReg-BirDea-Excur', e, '.png', sep='')
   if(Plot){png(FigName)}
   plot(BD$Time, BD$Count, type='s', lwd = 2, xlab='', ylab='')
   abline(h=M, lwd=3, lty=3, col=4)
   for (i in (1:e)){
      lines(c(BD$Time[In[i]], BD$Time[Out[i]]), c(M, M), lwd=6, col=3)
      }
   abline(v=BD$Time[In[e]]+c(0, l), lty=2, lwd=2, col=8)
   if(Plot){dev.off()}
   if(BD$Time[Out[i]]-BD$Time[In[i]] > l){test=T}
}

FigName = 'Fig-MinReg-BirDea-OK.png'
if(Plot){png(FigName)}
plot(BD$Time, BD$Count, type='s', lwd = 2, xlab='', ylab='')
abline(h=M, lwd=3, lty=3, col=4)
for (i in (1:e)){
   lines(c(BD$Time[In[i]], BD$Time[Out[i]]), c(M, M), lwd=6, col=3)
}
lines(c(BD$Time[In[e]], BD$Time[Out[e]]), c(M, M), lwd=6, col=2)
abline(v=BD$Time[In[e]]+c(0, l), lty=2, lwd=2, col=2)
if(Plot){dev.off()}

