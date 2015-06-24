# ==============================================================================
# Script to draw plots for all results in black and white
# ==============================================================================
setwd("~/Dropbox/Git_Repos/Fisheries_R_Scripts/Blue_Swimmer_Crab_R/Plot_Gen")
source("bsc_cm_plot_func.R")
source("bsc_simulation_plot_func.R")
# ------------------------------------------------------------------------------
# Generate the simulated data and then use the plot function to plot sim results 
# ------------------------------------------------------------------------------
sim.true.pars <- c(1, 2, 2, 190, 75.5, 40, 0.03)
pars          <- c(1.02, 1.84, 2.01, 189, 75.5, 40.2, 0.03)
LL            <- read.table("bsc_sim_lengths.txt")
LL            <- LL[, 1]
MM            <- read.table("bsc_sim_months.txt")
MM            <- MM[, 1]
BscPlotSim(pars)
# ------------------------------------------------------------------------------
# Generate the plots for the solutions to the real data
# ------------------------------------------------------------------------------
# Remove any objects to clear the slate
rm(list = ls( ))
# Source the function files needed
# --------------------------------
setwd("../Constrained_Max_Var_New")
source("bsc_cm_seas_integral_func.R")
source("bsc_cm_mean_length_func.R")
source("bsc_cm_seas_root_func.R")
source("bsc_cm_variance_func_ricker.R")
source("bsc_cm_log_likelihood_func.R")
source("bsc_cm_pi_calc_func.R")
source("bsc_cm_mean_var_optim_func.R")
# Data preliminaries
# ------------------
# Set the working directories
setwd("~/Dropbox/AAUni/APhD/Blueswimmer/RcodesData/Blue_Swimmer_Crab_Code_Sim/BSC_R_code_best/Diff_Variance_Function")
# ------------------------------------------------------------------------------
# Read in the data set on the asymptotic males that was gathered through pots
# Pull out the lengths and the dates from these data files
# Pull out the year and month information from these dates
# Read in the trawl data on males that contains juvenile recruitment and adults
# Pull out the lengths and the dates from these data files
# ------------------------------------------------------------------------------
lfd.big.males.females         <- read.table("LFD_bigMalesFem", header = T)
lfd.big.males.females.dates   <- as.Date(lfd.big.males.females$Date, "%d/%m/%y")
lfd.big.males.females.lengths <- lfd.big.males.females$Carapacewidth
lfd.big.males.females.dates.year   <- format(lfd.big.males.females.dates, '%Y')
lfd.big.males.females.dates.months <- format(lfd.big.males.females.dates, '%m')
lfd.trawl.males.females            <- read.table("LFD")
lfd.trawl.males.females.dates      <- as.Date(lfd.trawl.males.females $Date, 
									          "%d/%m/%Y")
lfd.trawl.males.females.lengths <- lfd.trawl.males.females$Carapace.width
# ------------------------------------------------------------------------------
# Concatenate the necessary elements from each file into a common 
# dates and lengths array
# ------------------------------------------------------------------------------
lfd.dates   <- c(lfd.trawl.males.females.dates,   lfd.big.males.females.dates)
lfd.lengths <- c(lfd.trawl.males.females.lengths, lfd.big.males.females.lengths)
lfd.year   <- format(lfd.dates, '%Y')
lfd.months <- format(lfd.dates, '%m')
# ------------------------------------------------------------------------------
# Pull out the years and months that we are interested in 
# i.e., those that don't contain recruitment
# ------------------------------------------------------------------------------
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
k0       <- 1.03                              # K0 average K
linf     <- 160                               # Asym length
mu.yr.1  <- 60                                # First month's average length yr 1
mu.yr.2  <- 80                                # First month's average length yr 2
theta.1  <- 0.9                               # Seasonality parameter 1
theta.1.comb <- 1.02346756   
theta.2.comb <- 0.32301298
max.contr    <- (1 / (2 * pi)) * 
                acos(theta.1.comb /
                (sqrt(theta.2.comb ^ 2 + 
                theta.1.comb ^ 2)))           # Calculates max of seas curve
theta.2   <- (theta.1 * (sqrt(1 - cos(2 * 
              pi * max.contr)^2))) /
              cos(2 * pi * max.contr)         # Theta 2 constrained by max 
# ------------------------------------------------------------------------------
# Combined
# ------------------------------------------------------------------------------
setwd("~/Dropbox/Git_Repos/Fisheries_R_Scripts/Blue_Swimmer_Crab_R/Plot_Gen")
source("bsc_plot_func.R")
pars           <- c(0.881, 1.06, 0.311, 176.3, 64.5, 79.9, 33.3, 0.0253, 0)
#pars          <- c(1.04, 0.989, 0.312, 172.7, 70.1, 81.6, 82.1, 0.0321, 0)
BscPlotNew(pars)
# ------------------------------------------------------------------------------
# Males
# ------------------------------------------------------------------------------
#Males
combined.sex <- c(lfd.trawl.males.females$Sex,   lfd.big.males.females$Sex)
males        <- which(combined.sex == 1)
lfd.dates    <- lfd.dates[males]
lfd.lengths  <- lfd.lengths[males]
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
k0       <- 1.03                              # K0 average K
linf     <- 160                               # Asym length
mu.yr.1  <- 60                                # First month's average length yr 1
mu.yr.2  <- 80                                # First month's average length yr 2
theta.1  <- 0.9                               # Seasonality parameter 1
theta.1.comb <- 1.02346756   
theta.2.comb <- 0.32301298
max.contr    <- (1 / (2 * pi)) * 
                acos(theta.1.comb /
                (sqrt(theta.2.comb ^ 2 + 
                theta.1.comb ^ 2)))           # Calculates max of seas curve
theta.2   <- (theta.1 * (sqrt(1 - cos(2 * 
              pi * max.contr)^2))) /
              cos(2 * pi * max.contr)         # Theta 2 constrained by max 
var.pars <- c(10, 0.015)                     # Variance fun parameter vector
pars     <- c(k0, theta.1, linf,             # Parameters in a vector
              mu.yr.1, mu.yr.2, var.pars) 
# Annnnddd plot
# -------------
pars          <- c(1.04, 0.989, 0.312, 172.7, 70.1, 81.6, 82.1, 0.0321, 1)
BscPlotNew(pars)
# ------------------------------------------------------------------------------
# Females
# ------------------------------------------------------------------------------
# Females
combined.sex <- c(lfd.trawl.males.females$Sex,   lfd.big.males.females$Sex)
females      <- which(combined.sex == 2)
lfd.dates    <- lfd.dates[females]
lfd.lengths  <- lfd.lengths[females]
#Pull out the year and month information from these dates
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
k0       <- 1.03                              # K0 average K
linf     <- 160                               # Asym length
mu.yr.1  <- 60                                # First month's average length yr 1
mu.yr.2  <- 80                                # First month's average length yr 2
theta.1  <- 0.9                               # Seasonality parameter 1
theta.1.comb <- 1.02346756   
theta.2.comb <- 0.32301298
max.contr    <- (1 / (2 * pi)) * 
                acos(theta.1.comb /
                (sqrt(theta.2.comb ^ 2 + 
                theta.1.comb ^ 2)))           # Calculates max of seas curve
theta.2   <- (theta.1 * (sqrt(1 - cos(2 * 
              pi * max.contr)^2))) /
              cos(2 * pi * max.contr)         # Theta 2 constrained by max 
#var.pars <- c(10, 0.015)                      # Variance fun parameter vector
#pars     <- c(k0, theta.1, linf,              # Parameters in a vector
#             mu.yr.1, mu.yr.2, var.pars) 
# Annnnddd plot
# -------------
pars          <- c(1.16, 1.10, 0.348, 156.2, 59.2, 80.5, 12.4, 0.0160, 1)
BscPlotNew(pars)