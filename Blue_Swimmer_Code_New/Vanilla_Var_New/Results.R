
# Initialise the parameters of the model 
# --------------------------------------
	
num.inds <- length(months)                    # Number of individuals we have
pi.1     <- rep(1/3, num.months)              # Pi mixing prop group 1
pi.2     <- rep(1/3, num.months)              # Pi mixing prop group 2
pi.3     <- (1 - (pi.1 + pi.2))               # Pi group 3. Diff from 1
k0       <- 1                                 # K0 average K
linf     <- 200                               # Asym length
mu.yr.1  <- 40                                # First month's average length yr 1
mu.yr.2  <- 60                                # First month's average length yr 2
theta.1  <- 1                                 # Seasonality parameter 1
theta.2  <- 0.2                               # Seasonality parameter 2
var.pars <- c(1, 0.01)                        # Variance fun parameter vector
pars     <- c(k0, theta.1, theta.2, linf,     # Parameters in a vector
              mu.yr.1, mu.yr.2, var.pars) 


[1] -6.754119e+04  9.684009e-06
[1]   0.88123275   1.05851108   0.31143653 176.25500525  64.48805451  79.94157463  33.34446030   0.02528065