# Check anova/adonis formulas

rm(list=ls());
# seed <- 1; set.seed(seed)

# Parms a dims
g <- 3; r <- 4; n <- g*r
mu <- 2*rnorm(g)
x <- rep(1:g, each=r)
y <- rnorm(n) + mu[x]
tab <- anova(lm(y ~ as.factor(x)))

# diffs
d <- as.matrix(dist(y, diag=TRUE))

# var
c(colSums(tab)[2],var(y)*(n-1), sum(d^2)/2/n)

# anova
c(tab[2, 2], sum(sapply(1:g, function(i){sum(d[which(x==i), which(x==i)]^2)/2/r})))
