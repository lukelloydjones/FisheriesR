# Function to optimise means using optim
# --------------------------------------
	
MeanVarOptim <- function(pars) {
  # Calculates the log likelihood 
  #
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
  #  The negative of the log likelihood for use in the OPTIM function
  
  # Declare each of the parameters to names unique to inside the fucntion
  
  k0.fun        <- pars[1]
  theta.1.fun   <- pars[2]
  theta.2.fun   <- pars[3]
  linf.fun      <- pars[4]
  mu.yr.1.fun   <- pars[5]
  mu.yr.2.fun   <- pars[6]
  var.par.1.fun <- pars[7]
  var.par.2.fun <- pars[8]

                   
  # It seems bad but we need to define the bottom bits as global variables
  # so that LogLikelihood can see them
  # Calculate the means given the current update of the parameters
  
  mean.2.yr <<- sapply(months.lst, MeanLength, k0 = k0.fun, theta.1 = theta.1.fun, 
                      theta.2 = theta.2.fun, linf = linf.fun , mu.yr.1 = mu.yr.1.fun, 
                      mu.yr.2 = mu.yr.2.fun, yrs.old = 2, str.mnth = 1)       
  mean.1.yr <<- sapply(months.lst, MeanLength, k0 = k0.fun, theta.1 = theta.1.fun, 
                      theta.2 = theta.2.fun, linf = linf.fun , mu.yr.1 = mu.yr.1.fun, 
                      mu.yr.2 = mu.yr.2.fun, yrs.old = 1, str.mnth = 1)
  mean.0.yr <<- sapply(months.lst, MeanLength, k0 = k0.fun, theta.1 = theta.1.fun, 
                      theta.2 = theta.2.fun, linf = linf.fun , mu.yr.1 = mu.yr.1.fun, 
                      mu.yr.2 = mu.yr.2.fun, yrs.old = 0, str.mnth = 1)
                    
  # Calculate the variances given the current update of the parameters
             
  var.2.yr  <<- sapply(mean.2.yr, BscVar, var.par.1 = var.par.1.fun,
                      var.par.2 = var.par.2.fun)
  var.1.yr  <<- sapply(mean.1.yr, BscVar, var.par.1 = var.par.1.fun,
                      var.par.2 = var.par.2.fun)
  var.0.yr  <<- sapply(mean.0.yr, BscVar, var.par.1 =  var.par.1.fun,
                      var.par.2 = var.par.2.fun)
                   
  # Calculate the log likelihood
               
  log.like.full <- sum(sapply(num.months.seq, LogLikelihood))

  # print(log.like.full)
  # Return the negative of the log likelihood
  
  return(-log.like.full)
}


