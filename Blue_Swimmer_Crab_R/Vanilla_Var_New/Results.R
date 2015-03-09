
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


# Males 


[1] -3.828742e+04  8.969706e-06
[1]   1.03695032   0.98890520 172.74136718  70.06643329  81.63191704  82.07053143   0.03212765   0.31210488
[9]  55.12485354

# Initialise the parameters of the model 
# --------------------------------------
	
num.inds <- length(months)                    # Number of individuals we have
pi.1     <- rep(1/3, num.months)              # Pi mixing prop group 1
pi.2     <- rep(1/3, num.months)              # Pi mixing prop group 2
pi.3     <- (1 - (pi.1 + pi.2))               # Pi group 3. Diff from 1
k0       <- 0.8                               # K0 average K
linf     <- 180                               # Asym length
mu.yr.1  <- 65                                # First month's average length yr 1
mu.yr.2  <- 80                                # First month's average length yr 2
theta.1  <- 1.05                              # Seasonality parameter 1
theta.1.comb <- 1.02346756   
theta.2.comb <- 0.32301298
max.contr    <- (1 / (2 * pi)) * 
                acos(theta.1.comb /
                (sqrt(theta.2.comb ^ 2 + 
                theta.1.comb ^ 2)))           # Calculates max of seas curve
theta.2   <- (theta.1 * (sqrt(1 - cos(2 * 
              pi * max.contr)^2))) /
              cos(2 * pi * max.contr)         # Theta 2 constrained by max 
var.pars <- c(33, 0.025)                        # Variance fun parameter vector
pars     <- c(k0, theta.1, linf,     # Parameters in a vector
              mu.yr.1, mu.yr.2, var.pars) 
              
              
# Females

[1] -2.783855e+04  2.538416e-03
[1]   1.1502842   1.0969879 156.5559785  59.1764504  80.4743537  12.6213429   0.0162413   0.3462165 155.4198819

# Initialise the parameters of the model 
# --------------------------------------
	
num.inds <- length(months)                    # Number of individuals we have
pi.1     <- rep(1/3, num.months)              # Pi mixing prop group 1
pi.2     <- rep(1/3, num.months)              # Pi mixing prop group 2
pi.3     <- (1 - (pi.1 + pi.2))               # Pi group 3. Diff from 1
k0       <- 1.03                               # K0 average K
linf     <- 160                              # Asym length
mu.yr.1  <- 60                                # First month's average length yr 1
mu.yr.2  <- 80                                # First month's average length yr 2
theta.1  <- 0.9                              # Seasonality parameter 1
theta.1.comb <- 1.02346756   
theta.2.comb <- 0.32301298
max.contr    <- (1 / (2 * pi)) * 
                acos(theta.1.comb /
                (sqrt(theta.2.comb ^ 2 + 
                theta.1.comb ^ 2)))           # Calculates max of seas curve
theta.2   <- (theta.1 * (sqrt(1 - cos(2 * 
              pi * max.contr)^2))) /
              cos(2 * pi * max.contr)         # Theta 2 constrained by max 
var.pars <- c(10, 0.015)                        # Variance fun parameter vector
pars     <- c(k0, theta.1, linf,     # Parameters in a vector
              mu.yr.1, mu.yr.2, var.pars) 


