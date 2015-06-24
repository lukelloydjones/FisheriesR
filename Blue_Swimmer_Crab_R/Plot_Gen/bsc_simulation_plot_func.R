 # Function to visualise the results of the BSC simulation
# --------------------------------------------------------
	
BscPlotSim <- function(pars) {
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
  
  # Load the libraries
  
  library(ggplot2)
  library(grid)
  source('~/Dropbox/AAUni/APhD/Blueswimmer/Crab_Simulation/bsc_sim_funcs.R')
  
  # Declare each of the parameters to names unique to inside the function
  
  K0    <- pars[1]
  T0    <- pars[2]
  T1    <- pars[3]
  LINF  <- pars[4]
  M01   <- pars[5]
  VV.1  <- pars[6]
  VV.2  <- pars[7]
  VV    <- c(VV.1, VV.2)
                   
  # Declare some plotting parameters
  
  XX <- 0 : 200
  xlabNam <- c("Jan.", "Feb.","Mar.","Apr.","May","Jun.","Jul.","Aug.", 
                  "Sep.", "Oct.", "Nov.", "Dec.")
  MMLIST <- seq(0, 11)
  # Plot of data and mixture fits
  MON <- 12
  plots <- list()
  NN <- length(MM)
  # Calculate the pis
  MU2A_AUX <- c()
  MUA_AUX  <- c()
  MUJ_AUX  <- c()
  for (mm in 1:MON)
	{
 	MU2A_AUX[mm] <- MU2A(mm,K0,T0,T1,LINF,M01)
  	MUA_AUX[mm]  <-  MUA(mm,K0,T0,T1,LINF,M01)
 	MUJ_AUX[mm]  <-  MUJ(mm,K0,T0,T1,LINF,M01)
	}
  TAU <- list()
  for (mm in 1:(MON))
	{
  		TAU[[mm]] <- matrix(1/3, NN, 3)
	}
  PI1 <- c(0.50, 0.40, 0.30, 0.20, 0.10, 0.05, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02)
  PI2 <- c(0.1, 0.15, 0.1, 0.08, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09)
  for (mm in 1:MON)
  {
    	
    # Calculate the inclusion probablities for group 1 and 2 which forms the top of the mixture model algorithm

    Top1<-PI1[mm] * dnorm(LL, MU2A_AUX[mm], sqrt(VAR2A(mm,VV,K0,T0,T1,LINF,M01))) 
    Top2<-PI2[mm] * dnorm(LL, MUA_AUX[mm],  sqrt( VARA(mm,VV,K0,T0,T1,LINF,M01)))								
    	
    # Calculate the sum of each of the classes which forms the bottom of the mixture model algorithm

    Bot<- Top1 + Top2 + (1-(PI1[mm]+PI2[mm]))*dnorm(LL,MUJ_AUX[mm],sqrt(VARJ(mm,VV,K0,T0,T1,LINF,M01)))  
    	
    	
    # Calculate the tau scores for each of the groups, where the third group is what's left over from 1 after the   
    #first 2

    TAU[[mm]][,1] <- Top1/Bot																	
    TAU[[mm]][,2] <- Top2/Bot
	TAU[[mm]][,3] <- 1-(TAU[[mm]][,1]+TAU[[mm]][,2])	
		
	# Work out the PIs for groups 1 and 2
												  		    
    PI1[mm] <- sum((MM==MMLIST[mm])*TAU[[mm]][,1])/sum(MM==MMLIST[mm])					
	PI2[mm] <- sum((MM==MMLIST[mm])*TAU[[mm]][,2])/sum(MM==MMLIST[mm])										
  }
  
  dev.new(width = 15, height = 17/1.618, units = "cm")
  plot.nrow <- 3
  plot.ncol <- 5
  pushViewport(viewport(layout = grid.layout(plot.nrow, plot.ncol)))
  row.seq <- rep(seq(1, plot.nrow), each  = plot.ncol)
  col.seq <- rep(seq(1, plot.ncol), times = plot.nrow)
  for (mm in 1:MON)
  {
    lengths.months <- LL[which(MM==MMLIST[mm])]
    gg.df <- as.data.frame(lengths.months)
    m <- ggplot(gg.df, aes(x = lengths.months)) + geom_histogram(aes(y = ..density..), 
                fill = rgb(0.5, 0.5, 0.5, 0.5),     
                binwidth = 7) + ylim(c(0, 0.07))
    den.mu2 <- PI1[mm] * dnorm(XX, MU2A(mm, K0, T0, T1, LINF, M01), 
                               sqrt(VAR2A(mm,VV,K0,T0,T1,LINF,M01)))
    den.mu1 <- PI2[mm] * dnorm(XX, MUA(mm, K0, T0, T1, LINF, M01), 
                               sqrt(VARA(mm,VV,K0,T0,T1,LINF,M01)))
    den.mu0 <- (1-(PI1[mm]+PI2[mm])) * dnorm(XX, MUJ(mm, K0, T0, T1, LINF, M01), 
                               sqrt(VARJ(mm,VV,K0,T0,T1,LINF,M01)))
    x.df.mu2 <- data.frame(XX, den.mu2)
    x.df.mu1 <- data.frame(XX, den.mu1)
    x.df.mu0 <- data.frame(XX, den.mu0)
    total.mixture <- den.mu2 + den.mu1 + den.mu0
	x.vals.tot    <- c(XX, XX, XX)
	tot.df        <- data.frame(x.vals.tot, total.mixture)
    plot.mm  <- m +  
          geom_line(data = x.df.mu2, aes(x = XX, y = den.mu2), linetype = 1, 
          size = 1) +
          geom_line(data = x.df.mu1, aes(x = XX, y = den.mu1), linetype = 2, 
          size = 1) +
          geom_line(data = x.df.mu0, aes(x = XX, y = den.mu0), linetype = 3, 
          size = 1) +
          geom_line(data = tot.df,   aes(x = x.vals.tot, y = total.mixture), 
          linetype = 1, size = 0.5) +
          xlab(xlabNam[mm]) + theme(panel.background = element_blank())
    print(c(row.seq[mm], col.seq[mm]))
    print(plot.mm, vp = viewport(layout.pos.row = row.seq[mm], layout.pos.col = col.seq[mm]))

  }
  
  # Plot the seasonal function  
  XX2  <- seq(0, 1, length.out = 100)
  YY   <- K0 + T0 * cos(2 * pi * XX2) + T1 * sin(2 * pi *XX2)
  YY.2 <- ((YY + abs(YY))/2)
  df.seas <- data.frame(XX2, YY.2)
  seas.plot <- ggplot(df.seas, aes(x = XX2, y = YY.2)) + 
               xlab("Fraction of a year since Jan. 1") + ylab("k (per year)")
  print(seas.plot + geom_line(size = 1), vp = viewport(layout.pos.row = 3, layout.pos.col = 3))
 
  # Plot the variance funtion
  
  var.linf <- VV[1] * XX * exp(-VV[2] * XX)
  df.var   <- data.frame(XX, var.linf)
  var.plot <- ggplot(df.var, aes(x = XX, y = var.linf)) + 
               xlab("Mean length (mm)") + ylab("Variance")
  print(var.plot + geom_line(size = 1), vp = viewport(layout.pos.row = 3, layout.pos.col = 4))

  
  # Plot the distribution of linf
  var.linf.dist <- rnorm(1000000, LINF, sqrt(VV[1] * LINF * exp(-VV[2] * LINF)))
  var.linf.dist <- as.data.frame(var.linf.dist)
  m <- ggplot(var.linf.dist, aes(x = var.linf.dist)) +
       xlab("Length (mm)") + ylab("Density")
  print(m + geom_line(size = 1, stat='density'), vp = viewport(layout.pos.row = 3, layout.pos.col = 5))

}
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  
	  