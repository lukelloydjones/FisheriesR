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

# g1 <- list(pi = , mean = , var = )
# my_groups <- list(g1 = g1, g2 = g1, g3 )	  
# logLik <- function(group, mm) group$pi[mm] * dnorm(lengths.sub, group$mean[mm], sqrt(group$var[mm]))

# Map(logLik, my_group, list_of_mm)