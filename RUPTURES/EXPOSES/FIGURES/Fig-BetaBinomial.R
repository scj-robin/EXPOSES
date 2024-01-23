# Figures -> beta-binomial
F_FigBetaBin <- function(a, b, y, n, pdf=F, col){
   # Plots for beta(a, b) prior and y success among n trials
   p.step=.001
   p = seq(p.step, 1-p.step, p.step)
   prior = dbeta(p, a, b)
   posterior = dbeta(p, a+y, b+n-y)
   ylim = c(0, 5)

   plot.name = gsub("\\.", "", paste('Fig-BetaBinomial-a', a, '-b', b, '-y', y, '-n', n, sep=''))
   if (pdf){pdf(paste(plot.name, '-prior.pdf', sep=''))}
   plot(p, prior, lwd=4, type='l', lty=2, main='', xlab='', ylab='', ylim=ylim, col=col)
   if(pdf){dev.off()}
   
   if (pdf){pdf(paste(plot.name, '-prior-data.pdf', sep=''))}
   plot(p, prior, lwd=4, type='l', lty=2, main='', xlab='', ylab='', ylim=ylim, col=col)
   abline(v=y/n, col=1, lty=2, lwd=4)
   if(pdf){dev.off()}

   if (pdf){pdf(paste(plot.name, '-all.pdf', sep=''))}
   plot(p, posterior, lwd=4, type='l', main='', xlab='', ylab='', ylim=ylim, col=col)
   abline(v=y/n, col=1, lty=2, lwd=4)
   lines(p, prior, lwd=4, lty=2, col=col)
   if(pdf){dev.off()}

   invisible(list(p=p, prior=prior, post=posterior))
}

y = 4; n = 20
BB = list()
BB[[1]] = F_FigBetaBin(1, 1, y, n, T, 2)
BB[[2]] = F_FigBetaBin(.5, .5, y, n, T, 3)
BB[[3]] = F_FigBetaBin(8, 3, y, n, T, 4)

pdf(paste('Fig-BetaBinomial-y', y, '-n', n, '-comp-prior.pdf', sep=''))
plot(BB1$p, BB1$post, main='', col=0, xlab='', ylab='', ylim=c(0, 5))
for (i in (1:3)){
   lines(BB[[i]]$p, BB[[i]]$prior, lwd=4, lty=2, col=1+i)
   lines(BB[[i]]$p, BB[[i]]$post, lwd=4, col=1+i)
}
abline(v=y/n, col=1, lty=2, lwd=4)
dev.off()

