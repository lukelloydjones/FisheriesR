# Function to calculate log likelihood
# ------------------------------------
	
LogLikelihood <- function(mm) {
  # Calculates the log likelihood 
  #
  # Args:
  #  k0:       Mean  k paramter for VB model
  #  theta.1:  Seasonality parameter 1
  #  theta.2:  Seasonality parameter 1
  #  time.var  time variable
  #
  # Returns:
  #  The log likelihood
  
  lengths.sub    <- lengths[which(months == months.lst[mm])]
  log.like.grp.1 <- pi.1[mm] * dnorm(LL2, MU2A_AUX[mm], sqrt(BscVar(mm, var.pars, K0, T0,T1,LINF,M01,M02)))
  log.like.grp.2 <- pi.2[mm] * dnorm(LL2,  MUA_AUX[mm], sqrt(BscVar(mm,  var.pars, K0, T0,T1,LINF,M01,M02)))
 
  log.like.grp.3 <- (1- (pi.1[mm] + pi.2[mm])) * dnorm(LL2,MUJ_AUX[mm],sqrt(VARJ(mm,VV,K0,T0,T1,LINF,M01,M02)))
  log.like <- sum(log(log.like.grp.1 + log.like.grp.2 + log.like.grp.3))
}


	  
	  for (mm in 1:MON) {
	
	  }