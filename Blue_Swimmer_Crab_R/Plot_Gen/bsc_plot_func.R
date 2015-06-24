# Function to visualise the data and the model within loop
# --------------------------------------------------------
	
BscPlotNew <- function(pars) {
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
  #    constrain: flag whether to use a constrained seasonal curve
  # Returns:
  #  A plot of the data and the mixture fits as well as the 
  
  # Declare each of the parameters to names unique to inside the function
  library(ggplot2)
  library(grid)
  k0        <- pars[1]
  theta.1   <- pars[2]
  theta.2   <- pars[3]
  linf      <- pars[4]
  mu.yr.1   <- pars[5]
  mu.yr.2   <- pars[6]
  var.par.1 <- pars[7]
  var.par.2 <- pars[8]
  constrain <- pars[9]
  if ( constrain == 1 ) {
  # If male or female we keep the maximum fixed so turn off thetas
  # above and turn those on below
  theta.2   <- (theta.1 * (sqrt(1 - cos(2 * pi * max.contr)^2))) /
                   cos(2 * pi * max.contr) 
  }
  #setwd("~/Dropbox/Git_Repos/Fisheries_R_Scripts/Blue_Swimmer_Crab_R/Constrained_Max_Var_New")
  #source("bsc_cm_seas_integral_func.R")
  #source("bsc_cm_mean_length_func.R")
  #source("bsc_cm_seas_root_func.R")
  #source("bsc_cm_variance_func_ricker.R")
  #source("bsc_cm_log_likelihood_func.R")
  #source("bsc_cm_pi_calc_func.R")
  #source("bsc_cm_mean_var_optim_func.R")

  # Declare some plotting parameters
  

  x.vals <- 0 : 200
  x.lab.name <- c("Feb., 1985", "Mar., 1985", "Apr., 1985", "May, 1985",
                  "Jun., 1985", "Jul., 1985", "Aug., 1985", "Feb., 1986",
                  "Mar., 1986", "Apr., 1986", "May, 1986")
  # Calculate the means again for the final likelihood update
  mean.2.yr <<- sapply(months.lst, MeanLength, k0 = k0, theta.1 = theta.1, 
                       theta.2 = theta.2, linf = linf , mu.yr.1 = mu.yr.1, 
                       mu.yr.2 = mu.yr.2, yrs.old = 2, str.mnth = 1)   
  print(mean.2.yr)    
  mean.1.yr <<- sapply(months.lst, MeanLength, k0 = k0, theta.1 = theta.1, 
                       theta.2 = theta.2, linf = linf , mu.yr.1 = mu.yr.1, 
                       mu.yr.2 = mu.yr.2, yrs.old = 1, str.mnth = 1)
  print(mean.1.yr)
  mean.0.yr <<- sapply(months.lst, MeanLength, k0 = k0, theta.1 = theta.1, 
                       theta.2 = theta.2, linf = linf , mu.yr.1 = mu.yr.1, 
                       mu.yr.2 = mu.yr.2, yrs.old = 0, str.mnth = 1)                  
  # Calculate the variances again for the final likelihood update
             
  var.2.yr  <<- sapply(mean.2.yr, BscVar, var.par.1 = var.par.1,
                       var.par.2 = var.par.2)
  var.1.yr  <<- sapply(mean.1.yr, BscVar, var.par.1 = var.par.1,
                       var.par.2 = var.par.2)
  var.0.yr  <<- sapply(mean.0.yr, BscVar, var.par.1 = var.par.1,
                       var.par.2 = var.par.2)
  print(mean.0.yr)
  pi.1     <- rep(1/3, num.months)              # Pi mixing prop group 1 
  pi.2     <- rep(1/3, num.months)
  pi.3     <- 1 - (pi.1 + pi.2)
  pi.all   <- sapply(num.months.seq, PiCalc)
  pi.1     <- pi.all[1, ]
  pi.2     <- pi.all[2, ]
  pi.3     <- pi.all[3, ]
  
  # Set up the plotting parameters
  dev.new(width = 15, height = 17/1.618, units = "cm")
  plot.nrow <- 3
  plot.ncol <- 5
  pushViewport(viewport(layout = grid.layout(plot.nrow, plot.ncol)))
  row.seq <- rep(seq(1, plot.nrow), each  = plot.ncol)
  col.seq <- rep(seq(1, plot.ncol), times = plot.nrow)
  print("I made it this far")
  # Plot of data and mixture fits

  for (mm in 1 : num.months) {
  	month.lengths <- lengths[which(months == months.lst[mm])]
    lengths.sub.df <- data.frame(month.lengths)
    m <- ggplot(lengths.sub.df, aes(x = month.lengths)) +
         geom_histogram(aes(y = ..density..), 
         fill = rgb(0.5, 0.5, 0.5, 0.5),     
         binwidth = 15) + ylim(c(0, 0.05))
    den.mu2 <- pi.1[mm] * dnorm(x.vals, mean.2.yr[mm], sqrt(var.2.yr[mm]))
    den.mu1 <- pi.2[mm] * dnorm(x.vals, mean.1.yr[mm], sqrt(var.2.yr[mm]))
    den.mu0 <- pi.3[mm] * dnorm(x.vals, mean.0.yr[mm], sqrt(var.0.yr[mm]))
    x.df.mu2 <- data.frame(x.vals, den.mu2)
    x.df.mu1 <- data.frame(x.vals, den.mu1)
    x.df.mu0 <- data.frame(x.vals, den.mu0)
	total.mixture <- den.mu2 + den.mu1 + den.mu0
	x.vals.tot    <- c(x.vals, x.vals, x.vals)
	tot.df        <- data.frame(x.vals.tot, total.mixture)
    plot.mm  <- m +  
          geom_line(data = x.df.mu2, aes(x = x.vals, y = den.mu2), 
          linetype = 1, size = 1) +
          geom_line(data = x.df.mu1, aes(x = x.vals, y = den.mu1), 
          linetype = 2, size = 1) +
          geom_line(data = x.df.mu0, aes(x = x.vals, y = den.mu0), 
          linetype = 3, size = 1) +
          geom_line(data = tot.df,   aes(x = x.vals.tot, y = total.mixture), 
          linetype = 1, size = 0.5) +
          xlab(x.lab.name[mm]) 
    print(c(row.seq[mm], col.seq[mm]))
    print(plot.mm, vp = viewport(layout.pos.row = row.seq[mm], 
          layout.pos.col = col.seq[mm]))
  }


  # Plot the seasonal function  
  
  x.vals.2  <- seq(0, 1, length.out = 100)
  yy <- k0 + theta.1 * cos(2 * pi * x.vals.2) +
               theta.2 * sin(2 * pi * x.vals.2)
  seas.func <- ((yy + abs(yy))/2)
  df.seas <- data.frame(x.vals.2, seas.func)
  seas.plot <- ggplot(df.seas, aes(x = x.vals.2 , y = seas.func)) + 
               xlab("Fraction of a year since Jan. 1") + ylab("k (per year)")
  print(seas.plot + geom_line(size = 1), vp = viewport(layout.pos.row = 3, 
        layout.pos.col = 2))
  
  
  # Plot the variance funtion
  
  var.linf.x <-  sapply(x.vals, BscVar, var.par.1 =  var.par.1,
                     var.par.2 = var.par.2)
  df.var   <- data.frame(x.vals, var.linf.x)
  var.plot <- ggplot(df.var, aes(x = x.vals, y = var.linf.x)) + 
               xlab("Mean length (mm)") + ylab("Variance")
  print(var.plot + geom_line(size = 1), vp = viewport(layout.pos.row = 3, 
        layout.pos.col = 3))
  
  
  # Plot the distribution of linf

  var.linf <- BscVar(var.par.1, var.par.2, linf)
  var.linf.vals <- rnorm(1000000, linf, sqrt(var.linf))
  var.linf.dist <- data.frame(var.linf.vals)
  m <- ggplot(var.linf.dist, aes(x = var.linf.vals)) +
       xlab("Length (mm)") + ylab("Density")
  print(m + geom_line(stat = "density", size = 1), 
        vp = viewport(layout.pos.row = 3, 
        layout.pos.col = 4))

  
}
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  