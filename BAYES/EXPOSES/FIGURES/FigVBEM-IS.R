# Pictures for VBEM-IS talk

rm(list=ls())
library(sna)
library(weights)
set.seed(1)
colStart <- 2; colTarget <- 4; 
# colStart <- 4; colTarget <- 2;
colCode <- paste0(colStart, 'to', colTarget)

# ######################################################################################
# # Tree data
source('/home/robin/PgmR/General/FunctionsMatVec.R')
DataDir = '/home/robin/Bureau/RECHERCHE/RESEAUX/GOF-Network/vbemapp/gof-network/Exemples/'
load(paste(DataDir, 'Tree.Rdata', sep=''))
n = nrow(Tree$Net)
# pdf('FigVBEM-IS-Tree-Network.pdf')
# gplot(Tree$Net, gmode='graph')
# # plot(graph_from_adjacency_matrix(Tree$Net, mode='undirected'))
# dev.off()
# pdf('FigVBEM-IS-Tree-GeneticDistance.pdf')
# image((1:n), (1:n), max(log10(Tree$EdgeCovar[, 1]))-F_Vec2Sym(log10(Tree$EdgeCovar[, 1])), xlab='', ylab='')
# dev.off()
# pdf('FigVBEM-IS-Tree-GeographicDistance.pdf')
# image((1:n), (1:n), max(Tree$EdgeCovar[, 2])-F_Vec2Sym(Tree$EdgeCovar[, 2]), xlab='', ylab='')
# dev.off()
# pdf('FigVBEM-IS-Tree-TaxonomicDistance.pdf')
# image((1:n), (1:n), max(Tree$EdgeCovar[, 3])-F_Vec2Sym(Tree$EdgeCovar[, 3]), xlab='', ylab='')
# dev.off()
 
# ######################################################################################
# # Importance sampling Beta-binomial case
n = 10; a.prior = b.prior = 1
# p = rbeta(1, a.prior, b.prior); Y = rbinom(1, n, p)
Y = n-1
a.post = a.prior+Y; b.post = b.prior+n-Y
a.wrong = b.post; b.wrong = a.post
p.list = seq(0, 1, by=.001)
# # plot(p.list, dbeta(p.list, a.post, b.post), main='', ylab='', xlab='', type='l', lwd=4, col=0)
# # lines(p.list, dbeta(p.list, a.prior, b.prior), lwd=4, col=2)
# # lines(p.list, dbeta(p.list, a.post, b.post), lwd=4, col=1)
# 
PlotIS <- function(a.IS, b.IS, title, a.target=a.post, b.target=b.post, M = 1e3){
   p.IS = rbeta(M, a.IS, b.IS);
   w.IS = dbeta(p.IS, a.target, b.target) / dbeta(p.IS, a.IS, b.IS)
   W.IS = w.IS / sum(w.IS)
   ESS.IS = mean(w.IS)^2/mean(w.IS^2)
   H = wtd.hist(p.IS, weight=w.IS, breaks=sqrt(M), xlim=c(0, 1), xlab='', ylab='',
                main=paste(title, ': ESS = ', signif(ESS.IS, 2), sep=''), cex.main=2.5, family='sans');
   w = mean(diff(H$breaks))
   lines(p.list, M*w*dbeta(p.list, a.target, b.target), lwd=4, col=colTarget)
   lines(p.list, M*w*dbeta(p.list, a.IS, b.IS), lwd=4, col=colStart, lty=2)
}
# 
# pdf(paste0('FigVBEM-IS-ISpost-', colCode, '.pdf'))
# PlotIS(a.post, b.post, 'posterior')
# dev.off()
# pdf(paste0('FigVBEM-IS-ISprior-', colCode, '.pdf'))
# PlotIS(a.prior, b.prior, 'prior')
# dev.off()
# pdf(paste0('FigVBEM-IS-ISwrong-', colCode, '.pdf'))
# PlotIS(b.post/2, a.post/2, 'bad')
# dev.off()
# pdf(paste0('FigVBEM-IS-ISgood-', colCode, '.pdf'))
# PlotIS(a.post/2, b.post/2, 'not too bad')
# dev.off()
# pdf(paste0('FigVBEM-IS-ISvb-', colCode, '.pdf'))
# PlotIS(10*a.post, 10*b.post, 'too peaked')
# dev.off()

######################################################################################
# Tempering
a.init = 2*b.post; b.init = 2*a.post
a.seq = exp(seq(log(a.init), log(a.post), length.out=5))
b.seq = exp(seq(log(b.init), log(b.post), length.out=5))
nstep = length(a.seq)-1

pdf(paste0('FigVBEM-IS-PropTarget-', colCode, '.pdf'))
plot(p.list, dbeta(p.list, a.init, b.init), main='', ylab='', xlab='', type='l', lwd=4, col=colStart)
lines(p.list, dbeta(p.list, a.post, b.post), lwd=4, col=colTarget)
dev.off();

pdf(paste0('FigVBEM-IS-Tempering-', colCode, '.pdf'))
plot(p.list, dbeta(p.list, a.init, b.init), main='', ylab='', xlab='', type='l', lwd=4, col=colStart)
lines(p.list, dbeta(p.list, a.post, b.post), lwd=4, col=colTarget)
for (s in (2:nstep)){
   lines(p.list, dbeta(p.list, a.seq[s], b.seq[s]), lwd=4, lty=2)
}
dev.off();

for (s in 1:nstep){
   pdf(paste0('FigVBEM-IS-Tempering-step', s, '-', colCode, '.pdf'))
   PlotIS(a.seq[s], b.seq[s], paste('step', s), a.seq[s+1], b.seq[s+1])
   dev.off()
}

######################################################################################
# Importance sampling Gaussian-Gaussian case
# mu ~ N(m.prior, s2.prior); Y_i iid ~ N(mu, s2)
# mu ~ N(m.post, s2.post); lambda = s2 / (s2/n + s2.prior); 
# m.post = lambda Ybar + (1 - lambda) m.prior; s2.post = 1 / (1/s2.prior + n /s2)
m.prior = 0; s.prior = 10; s2.prior = s.prior^2
n = 10; s = 1; s2 = s^2; Ybar = 2
lambda = s2.prior / (s2/n + s2.prior); 
m.post = lambda * Ybar + (1 - lambda) * m.prior
s2.post = 1 / (1/s2.prior + n /s2); s.post = sqrt(s2.post)
M.max = 5; M.list = seq(m.post-M.max, m.post+M.max, length.out=1000)

PlotISgauss <- function(m.IS, s.IS, title, m.target=m.post, s.target=s.post, M = 1e3){
   M.IS = rnorm(M, m.IS, s.IS); 
   w.IS = dnorm(M.IS, m.target, s.target) / dnorm(M.IS, m.IS, s.IS)
   W.IS = w.IS / sum(w.IS)
   ESS.IS = mean(w.IS)^2/mean(w.IS^2)
   H = wtd.hist(M.IS, weight=w.IS, breaks=sqrt(M), xlim=c(min(M.list), max(M.list)), plot=F); 
   w = mean(diff(H$breaks))
   f.target = dnorm(M.list, m.target, s.target)
   f.IS = dnorm(M.list, m.IS, s.IS)
   wtd.hist(M.IS, weight=w.IS, breaks=sqrt(M), xlim=c(min(M.list), max(M.list)), 
                ylim = M*w*c(0, max(max(f.target), max(f.IS))), xlab='', ylab='', 
                main=paste(title, ': ESS = ', signif(ESS.IS, 2), sep=''), cex.main=2.5, family='sans'); 
   lines(M.list, M*w*f.target, lwd=4, col=colTarget)
   lines(M.list, M*w*f.IS, lwd=4, col=colStart, lty=2)
}

# pdf(paste0('FigVBEM-IS-GaussMean-ISpost-', colCode, '.pdf'))
PlotISgauss(m.post, s.post, 'posterior')
# dev.off()
# pdf(paste0('FigVBEM-IS-GaussMean-ISprior-', colCode, '.pdf'))
PlotISgauss(m.prior, s.prior, 'prior')
# dev.off()
# pdf(paste0('FigVBEM-IS-GaussMean-ISwrong-', colCode, '.pdf'))
PlotISgauss(m.post/2, 5*s.post, 'bad')
# dev.off()
# pdf(paste0('FigVBEM-IS-GaussMean-ISgood-', colCode, '.pdf'))
PlotISgauss(1.1*m.post, 2*s.post, 'not too bad')
# dev.off()
# pdf(paste0('FigVBEM-IS-GaussMean-ISvb-', colCode, '.pdf'))
PlotISgauss(m.post, s.post/2, 'too peaked')
# dev.off()
