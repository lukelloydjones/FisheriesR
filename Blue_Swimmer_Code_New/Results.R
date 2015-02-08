# RESULTS
# -------

# DATE: 07/02/2015

# females and males
# starting pars

num.inds <- length(months)                    # Number of individuals we have
pi.1     <- rep(1/3, num.months)              # Pi mixing prop group 1
pi.2     <- rep(1/3, num.months)              # Pi mixing prop group 2
pi.3     <- (1 - (pi.1 + pi.2))               # Pi group 3. Diff from 1
k0       <- 3                                 # K0 average K
linf     <- 200                               # Asym length
mu.yr.1  <- 40                                # First month's average length yr 1
mu.yr.2  <- 60                                # First month's average length yr 2
theta.1  <- -0.8                              # Seasonality parameter 1
theta.2  <- 0.2                               # Seasonality parameter 2
var.pars <- c(5, 1/100, 3, 1)                 # Variance fun parameter vector
pars     <- c(k0, theta.1, theta.2, linf,     # Parameters in a vector
              mu.yr.1, mu.yr.2, var.pars) 


[1] -6.752697e+04  7.041526e-06
 [1]   0.89809158   1.02346756   0.32301298 175.83800282  64.98820684  80.08088982  36.26829730   0.02588519
 [9] -41.56944464 356.30263720
 
 
 # Males
 
 num.inds <- length(months)                    # Number of individuals we have
pi.1     <- rep(1/3, num.months)              # Pi mixing prop group 1
pi.2     <- rep(1/3, num.months)              # Pi mixing prop group 2
pi.3     <- (1 - (pi.1 + pi.2))               # Pi group 3. Diff from 1
k0       <- -0.5                               # K0 average K
linf     <- 220                               # Asym length
mu.yr.1  <- 40                                # First month's average length yr 1
mu.yr.2  <- 60                                # First month's average length yr 2
theta.1  <- -0.1                              # Seasonality parameter 1
theta.2  <- 0.4                               # Seasonality parameter 2
var.pars <- c(5, 1/100, 3, 1)                 # Variance fun parameter vector
pars     <- c(k0, theta.1, theta.2, linf,     # Parameters in a vector
              mu.yr.1, mu.yr.2, var.pars) 