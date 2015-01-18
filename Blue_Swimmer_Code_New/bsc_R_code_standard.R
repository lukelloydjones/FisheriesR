###############################################################################
###############################################################################
###############################################################################
###													    	       ###
###     SCRIPT TO IMPLEMENT ALGORITHM FOR SOLVING FOR GROWTH PARAMETERS     ###
###	   VIA THE METHOD OUTLINED IN LLOYD-JONES ET AL. USES AN MM ALGORI     ###
###     THM AND THE OPTIM FUNCTION TO OPTIMISE OVER LENGTH FREQUENCY DA     ###
###     TA SETS.												  ###
###															  ###	
###############################################################################
###############################################################################
###############################################################################


# Remove any objects to clear the slate

rm(list = ls( ))


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
end.month.yr.1 <- 2

str.month.yr.2 <- 8
end.month.yr.2 <- 5

lfd.85.feb.aug <- which((lfd.year == '1985') & 
				  (as.numeric(lfd.months) %in% (str.month.yr.1:end.month.yr.1)))
lfd.86.feb.may <- which((lfd.year == '1986') &
				  (as.numeric(lfd.months) %in% (str.month.yr.2:end.month.yr.2)))


# Initialise the data for the model 
# ---------------------------------

num.months      <- 11													
lfd.lengths.sub <- lfd.lengths[c(lfd.85.feb.aug, lfd.86.feb.may)]		
months.85		<- as.numeric(lfd.months[lfd.85.feb.aug]) - 1 # Set January to be 0 th month
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
theta.1  <- 0.2								# Seasonality parameter 1
theta.2  <- 0.1								# Seasonality parameter 2
var.pars <- c(5, 1/100, 3, 1)				# Variance function parameter vector

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	