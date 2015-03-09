# Function to calculate log likelihood
# ------------------------------------
	
LogLikelihood <- function(mm) {
  # Calculates the log likelihood 
  #
  # Args:
  #  mm: Index for which month we are up to. Not the month itself
  #
  # Returns:
  #  The log likelihood for each month

  lengths.sub    <- lengths[which(months == months.lst[mm])]
  log.like.grp.1 <- pi.1[mm] * dnorm(lengths.sub, mean.2.yr[mm], sqrt(var.2.yr[mm]))
  log.like.grp.2 <- pi.2[mm] * dnorm(lengths.sub, mean.1.yr[mm], sqrt(var.1.yr[mm]))
  log.like.grp.3 <- pi.3[mm] * dnorm(lengths.sub, mean.0.yr[mm], sqrt(var.0.yr[mm]))
  log.like       <- sum(log(log.like.grp.1 + log.like.grp.2 + log.like.grp.3))
  
  
  return(log.like)

}

