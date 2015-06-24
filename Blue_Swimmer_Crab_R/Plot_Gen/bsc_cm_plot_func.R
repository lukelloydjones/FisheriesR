# Function to visualise the data and the model within loop
# --------------------------------------------------------
	
BscPlot <- function(pars) {
  # Plots the data and the mixture fits as well as the 
  # the seasonal function, variance func, linf distribution
  # Args:
  #  pars: Vector containing the parameters to optimise
  #    k0:        Mean  k paramter for VB model
  #    theta.1:   Seasonality parameter 1
  #    theta.2:   Seasonality parameter 1
  #    linf:      Asymptotic length
  #    mu.yr.1:   First month's average length yr 1
  #    mu.yr.2:   First month's average length yr 2T0,T1,LINF,M01,M02
  #    var.pars.x The variance parameters
  # Returns:
  #  A plot of the data and the mixture fits as well as the 
  
  # Declare each of the parameters to names unique to inside the function
  
  k0.fun        <- pars[1]
  theta.1.fun   <- pars[2]
  linf.fun      <- pars[3]
  mu.yr.1.fun   <- pars[4]
  mu.yr.2.fun   <- pars[5]
  var.par.1.fun <- pars[6]
  var.par.2.fun <- pars[7]

  # If male or female we keep the maximum fixed so turn off thetas
  # above and turn those on below
  theta.2.fun   <- (theta.1.fun * (sqrt(1 - cos(2 * pi * max.contr)^2))) /
                   cos(2 * pi * max.contr) 
  
                   
  # Declare some plotting parameters
  
  par(mfrow = c(3, 5))
  x.vals <- 0 : 200
  x.lab.name <- c("Feb., 1985", "Mar., 1985", "Apr., 1985", "May, 1985",
                  "Jun., 1985", "Jul., 1985", "Aug., 1985", "Feb., 1986",
                  "Mar., 1986", "Apr., 1986", "May, 1986")
  
  # Plot of data and mixture fits
  
  for (mm in 1 : num.months) {
  	    
    lengths.sub    <- lengths[which(months == months.lst[mm])]
    hist(lengths.sub, breaks = 30, prob = T,
    xlim = c(20, 210), ylim = c(0, 0.055), xlab = x.lab.name[mm], 
    ylab = "", main = '')
    lines(x.vals, pi.1[mm] * dnorm(x.vals, mean.2.yr[mm], 
         sqrt(var.2.yr[mm])), col = 'red',   lwd = 2.5)
    lines(x.vals, pi.2[mm] * dnorm(x.vals, mean.1.yr[mm], 
         sqrt(var.1.yr[mm])), col = 'green', lwd = 2.5)
    lines(x.vals, pi.3[mm] * dnorm(x.vals, mean.0.yr[mm], 
         sqrt(var.0.yr[mm])), col = 'blue',  lwd = 2.5)
	         
  }


  # Plot the seasonal function  
  
  x.vals.2  <- seq(0, 1, length.out = 100)
  seas.func <- k0.fun + theta.1.fun * cos(2 * pi * x.vals.2) +
               theta.2.fun * sin(2 * pi * x.vals.2)
  plot(x.vals.2, ((seas.func + abs(seas.func)) / 2), type = 'l',
      xlab = "Fraction of a year since Jan 1", ylab = "k (per year)", 
      lwd = 2.5)
  
  
  # Plot the variance funtion
  
  var.linf.x <-  sapply(x.vals, BscVar, var.par.1 =  var.par.1,
                     var.par.2 = var.par.2)
  plot(x.vals, var.linf.x, type = 'l', xlab = "Mean length (mm)", 
      ylab = "Variance", lwd = 2.5)
  
  
  # Plot the distribution of linf

  var.linf <- BscVar(var.par.1.fun, var.par.2.fun, linf.fun)
  plot(density(rnorm(1000000, linf, sqrt(var.linf))), 
      xlab = "Length (mm)", main = "", lwd = 2.5)
  
}
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  