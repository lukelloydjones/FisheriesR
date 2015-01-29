# Variance function for calcuating the variance for the current month
# -------------------------------------------------------------------


BscVar <- function(month, var.par.1, var.par.2, var.par.3, var.par.4, mean.length) {
	
  # Computes the variance of the distribution for the current month
  # based on a bespoke varinace function that is positive in x and y.
  #
  # Args:
  #  month:    Current month to calulate mean at. Numbered from Jan=0
  #  k0:       Mean  k paramter for VB model
  #  theta.1:  Seasonality parameter 1
  #  theta.2:  Seasonality parameter 1
  #  linf:     Asymptotic length
  #  mu.yr.1:  First month's average length yr 1
  #  mu.yr.2:  First month's average length yr 2
  #  mean.length: Mean length for current month
  #
  # Returns:
  #  Variance for each month
  
  variance = max(var.par.1 * mean.length * exp(-var.par.2 * mean.length) + 
  exp(var.par.3 * (1-exp(-var.par.4 * mean.length))), 1)
  
  return(variance)
  }