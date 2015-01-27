###############################################################################
###############################################################################
###############################################################################
###                                                                         ###
###     SCRIPT TO IMPLEMENT ALGORITHM FOR SOLVING FOR GROWTH PARAMETERS     ###
###	    VIA THE METHOD OUTLINED IN LLOYD-JONES ET AL. USES AN MM ALGORI     ###
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

setwd("~/Dropbox/Git_Repos/Fisheries_R_Scripts/Blue_Swimmer_Code_New")
source("bsc_seas_integral_func.R")
source("bsc_mean_length_func.R")
source("bsc_seas_root_func.R")
source("bsc_variance_func.R")


# Data preliminaries
# ------------------

# Set the working directories

setwd("~/Dropbox/Git_Repos/Fisheries_R_Scripts/Blue_Swimmer_Crab_Code_Sim/")
setwd("BSC_R_code_best/Diff_Variance_Function/")


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

num.months      <- 11													
lfd.lengths.sub <- lfd.lengths[c(lfd.85.feb.aug, lfd.86.feb.may)]		
months.85		<- as.numeric(lfd.months[lfd.85.feb.aug]) - 1 # Jan = 0th month
months.86 		<- as.numeric(lfd.months[lfd.86.feb.may]) + 11
months 			<- c(months.85, months.86)
months.lst      <- as.numeric(names(table(months)))	


# Initialise the parameters of the model 
# --------------------------------------
	
num.inds <- length(months)					# Number of individuals we have
pi.1     <- rep(1/3, num.months)			# Pi mixing prop group 1
pi.2     <- rep(1/3, num.months)			# Pi mixing prop group 2
k0       <- 1		 						# K0 average K
linf     <- 200								# Asym length
mu.yr.1  <- 40								# First month's average length yr 1
mu.yr.2  <- 40								# First month's average length yr 2
theta.1  <- 1.0	     						# Seasonality parameter 1
theta.2  <- 1.0								# Seasonality parameter 2
var.pars <- c(5, 1/100, 3, 1)				# Variance function parameter vector


# Initialise the taus 
# -------------------	

# There may be a better way to do this but does each tau get calculated for
# each individual in each month

tau.list <- list()
for (mm in 1:num.months)
	{
  		tau.list[[mm]] <- matrix(1/3, num.inds, 3)
	}

	
# Initialise the means for each of the three groups
# -------------------------------------------------

# Function paramters are month, k0, theta 1, theta 2, Linf, mean yr 1, 
# mean yr 2, years old, start month or first month included

mean.2.yr <- c()
mean.1.yr <- c()
mean.0.yr <- c()

for (mm in months.lst)
  {
    mean.2.yr[mm] <- MeanLength(mm, k0, theta.1, theta.2, linf, mu.yr.1, mu.yr.2, 2, 1)
    mean.1.yr[mm] <- MeanLength(mm, k0, theta.1, theta.2, linf, mu.yr.1, mu.yr.2, 1, 1)
    mean.0.yr[mm] <- MeanLength(mm, k0, theta.1, theta.2, linf, mu.yr.1, mu.yr.2, 0, 1)	
  }

	
	
	> MU2A_AUX
 [1] 178.3464 180.3816 182.0983 183.5115 184.6750 185.6661 186.5674 178.3464 180.3816 182.0983 183.5115
> MUA_AUX
 [1] 141.1393 146.6717 151.3381 155.1795 158.3424 161.0365 163.4865 141.1393 146.6717 151.3381 155.1795
> MUJ_AUX
 [1]  40.00000  55.03863  67.72336  78.16521  86.76302  94.08611 100.74605  40.00000  55.03863  67.72336
[11]  78.16521
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	