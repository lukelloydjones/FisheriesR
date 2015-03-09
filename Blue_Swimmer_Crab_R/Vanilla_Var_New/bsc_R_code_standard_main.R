###############################################################################
###############################################################################
###############################################################################
###                                                                         ###
###     SCRIPT TO IMPLEMENT ALGORITHM FOR SOLVING FOR GROWTH PARAMETERS     ###
###     VIA THE METHOD OUTLINED IN LLOYD-JONES ET AL. USES AN MM ALGORI     ###
###     THM AND THE OPTIM FUNCTION TO OPTIMISE OVER LENGTH FREQUENCY DA     ###
###     TA SETS.                                                            ###
###                                                                         ###	
###############################################################################
###############################################################################
###############################################################################

# Remove any objects to clear the slate

rm(list = ls( ))

# Source the function files needed
# --------------------------------

setwd("~/Dropbox/Git_Repos/Fisheries_R_Scripts/Blue_Swimmer_Code_New/Vanilla_Var_New")
source("bsc_seas_integral_func.R")
source("bsc_mean_length_func.R")
source("bsc_seas_root_func.R")
#source("bsc_variance_func.R")
source("bsc_variance_func_ricker.R")
source("bsc_log_likelihood_func.R")
source("bsc_pi_calc_func.R")
source("bsc_mean_var_optim_func.R")
source("bsc_plot_func.R")

# Data preliminaries
# ------------------

# Set the working directories

#setwd("~/Dropbox/Git_Repos/Fisheries_R_Scripts/Blue_Swimmer_Crab_Code_Sim/")
setwd("~/Dropbox/AAUni/APhD/Blueswimmer/RcodesData/Blue_Swimmer_Crab_Code_Sim/BSC_R_code_best/Diff_Variance_Function")


# Read in the data set on the asymptotic males that was gathered through pots

lfd.big.males.females <- read.table("LFD_bigMalesFem", header = T)


# Pull out the lengths and the dates from these data files

lfd.big.males.females.dates   <- as.Date(lfd.big.males.females$Date, "%d/%m/%y")
lfd.big.males.females.lengths <- lfd.big.males.females$Carapacewidth


# Pull out the year and month information from these dates

lfd.big.males.females.dates.year   <- format(lfd.big.males.females.dates, '%Y')
lfd.big.males.females.dates.months <- format(lfd.big.males.females.dates, '%m')


# Read in the trawl data on males that contains juvenile recruitment and adults

lfd.trawl.males.females <- read.table("LFD")


# Pull out the lengths and the dates from these data files

lfd.trawl.males.females.dates   <- as.Date(lfd.trawl.males.females $Date, 
									"%d/%m/%Y")
lfd.trawl.males.females.lengths <- lfd.trawl.males.females$Carapace.width


# Concatenate the necessary elements from each file into a common 
# dates and lengths array

lfd.dates   <- c(lfd.trawl.males.females.dates,   lfd.big.males.females.dates)
lfd.lengths <- c(lfd.trawl.males.females.lengths, lfd.big.males.females.lengths)


# Subset for males or females

# Males

# combined.sex <- c(lfd.trawl.males.females$Sex,   lfd.big.males.females$Sex)
# males        <- which(combined.sex == 1)
# lfd.dates    <- lfd.dates[males]
# lfd.lengths  <- lfd.lengths[males]

# Females

# combined.sex <- c(lfd.trawl.males.females$Sex,   lfd.big.males.females$Sex)
# females      <- which(combined.sex == 2)
# lfd.dates    <- lfd.dates[females]
# lfd.lengths  <- lfd.lengths[females]


# Pull out the year and month information from these dates

lfd.year   <- format(lfd.dates, '%Y')
lfd.months <- format(lfd.dates, '%m')


# Pull out the years and months that we are interested in 
# i.e., those that don't contain recruitment

str.month.yr.1 <- 2
end.month.yr.1 <- 8

str.month.yr.2 <- 2
end.month.yr.2 <- 5

lfd.85.feb.aug <- which((lfd.year == '1985') & 
				  (as.numeric(lfd.months) %in% (str.month.yr.1:end.month.yr.1)))
lfd.86.feb.may <- which((lfd.year == '1986') &
				  (as.numeric(lfd.months) %in% (str.month.yr.2:end.month.yr.2)))


# Initialise the data for the model 
# ---------------------------------

num.months       <- 11													
lengths          <- lfd.lengths[c(lfd.85.feb.aug, lfd.86.feb.may)]		
months.85        <- as.numeric(lfd.months[lfd.85.feb.aug]) - 1 # Jan = 0th month
months.86        <- as.numeric(lfd.months[lfd.86.feb.may]) + 11
months           <- c(months.85, months.86)
months.lst       <- as.numeric(names(table(months)))	
num.months.seq   <- seq(1, num.months)

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


# Initialise the likelihood and set tolerence
# -------------------------------------------

MeanVarOptim(pars)	
log.like.full <- -10e5
tol           <- 10e-6
log.like.old  <- -Inf	

	
# Run while loop over procedure until convergence
# -----------------------------------------------	

while (log.like.full - log.like.old > tol) {

  # Shift the current likelihood to the old likelihood
  
  log.like.old <- log.like.full 
  
  
  # Calculate the pi for each group in each month
  # ---------------------------------------------

  # Returns a vector of pi with each column representing a month
  # and wach row the groups. Row 1 the largest. Row 2 the yr olds
  # and row 3 the juveniles

  
  pi.all <- sapply(num.months.seq, PiCalc)
  pi.1   <- pi.all[1, ]
  pi.2   <- pi.all[2, ]
  pi.3   <- pi.all[3, ]


  # Optimise the parameters for the means
  # -------------------------------------
  
  # Initialise and optimise
  
  #pars            <- c(k0, theta.1, theta.2, linf, mu.yr.1, mu.yr.2, var.pars) 
  optim.means.var <- optim(pars, MeanVarOptim, control = list(maxit = 100000))	
  
  # Ask if optim converged
  
  print("Did optim converge?")
  print(optim.means.var$convergence)
  
  # Re-define the global parameters
  
  k0        <- optim.means.var$par[1]
  theta.1   <- optim.means.var$par[2]
  theta.2   <- optim.means.var$par[3]
  linf      <- optim.means.var$par[4]
  mu.yr.1   <- optim.means.var$par[5]
  mu.yr.2   <- optim.means.var$par[6]
  var.par.1 <- optim.means.var$par[7]
  var.par.2 <- optim.means.var$par[8]
  pars      <- optim.means.var$par
  
  # If male or female we keep thetas fixed so turn off thetas
  # above and turn those on below. Look in bsc_mean_var_func.R
  # for more details

  
  # Calculate the means again for the final likelihood update
  
  mean.2.yr <- sapply(months.lst, MeanLength, k0 = k0, theta.1 = theta.1, 
                     theta.2 = theta.2, linf = linf , mu.yr.1 = mu.yr.1, 
                     mu.yr.2 = mu.yr.2, yrs.old = 2, str.mnth = 1)       
  mean.1.yr <- sapply(months.lst, MeanLength, k0 = k0, theta.1 = theta.1, 
                     theta.2 = theta.2, linf = linf , mu.yr.1 = mu.yr.1, 
                     mu.yr.2 = mu.yr.2, yrs.old = 1, str.mnth = 1)
  mean.0.yr <- sapply(months.lst, MeanLength, k0 = k0, theta.1 = theta.1, 
                     theta.2 = theta.2, linf = linf , mu.yr.1 = mu.yr.1, 
                     mu.yr.2 = mu.yr.2, yrs.old = 0, str.mnth = 1)
                    
  # Calculate the variances again for the final likelihood update
             
  var.2.yr  <- sapply(mean.2.yr, BscVar, var.par.1 = var.par.1,
                     var.par.2 = var.par.2)
  var.1.yr  <- sapply(mean.1.yr, BscVar, var.par.1 = var.par.1,
                     var.par.2 = var.par.2)
  var.0.yr  <- sapply(mean.0.yr, BscVar, var.par.1 =  var.par.1,
                     var.par.2 = var.par.2)
                   
  
  # Evaluate the likelihood
  
  log.like.full <- sum(sapply(num.months.seq, LogLikelihood))
  
  
  # Give a plot of the current state of the model versus the data
  
  BscPlot(pars)


  # Print out the loglikelihood, tolerance, and parameters
  
  print(c(log.like.full, log.like.full - log.like.old))
  print(pars)
  
}
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	