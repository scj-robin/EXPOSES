# Posterior dist P(K|Y)

# Function
Prob <- function(LB){
   Prob = LB - mean(LB)
   Prob = exp(Prob)
   Prob = Prob / sum(Prob)
   return(Prob)
}

#Blog : 
#* sans covariable :  res$criterion_best
LB = c(-5062.166, -4668.137, -4430.175, -4174.351, -4046.210, -3934.743, -3777.563, -3708.774, -3661.853, -3625.978, -3603.810, -3599.215, -3583.430, -3596.835, -3618.848, -3621.867)
Prob(LB); plot(Prob(LB))
# * avec 1 covariable (même parti politique ou pas) :  
LB = c(-3920.037, -3797.661, -3531.063, -3529.182, -3435.559, -3479.923, -3526.242, -3541.893, -3578.929, -3546.849, -3505.651, -3510.716, -3533.796, -3532.282, -3526.279, -3546.953)
Prob(LB); plot(Prob(LB))

# Tree : 
# * sans covariable : 
LB = c(-757.6973, -573.7880, -562.1126, -500.2720, -463.4279, -460.3010, -459.0511, -460.3215, -456.9404, -464.1103, -472.0755, -471.6916)
Prob(LB); plot(Prob(LB))
# * avec les 3 covariables : 
LB = c(-808.1909, -534.1660, -454.3460, -491.0145, -472.9987, -477.4055, -486.0183, -488.8998, -480.7952, -497.4202, -503.6090, -506.7200)
Prob(LB); plot(Prob(LB))

Prob <- function(LB){
   Prob = LB - mean(LB)
   Prob = exp(Prob)
   Prob = Prob / sum(Prob)
   return(Prob)
}
plot(Prob(LB))
