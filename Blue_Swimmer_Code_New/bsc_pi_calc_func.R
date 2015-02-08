# Function to calculate pis for each month
# ----------------------------------------
	
PiCalc <- function(mm) {
  # Calculates the log likelihood 
  #
  # Args:
  #  mm: Index for which month we are up to. Not the month itself
  #
  # Returns:
  #  A set of pis for each group for the month
  
  # Subset the lengths
  lengths.sub    <- lengths[which(months == months.lst[mm])]
  
  # Calculate the inclusion probablities 
  
  top.1 <- pi.1[mm] * dnorm(lengths.sub, mean.2.yr[mm], sqrt(var.2.yr[mm]))			
  top.2 <- pi.2[mm] * dnorm(lengths.sub, mean.1.yr[mm], sqrt(var.1.yr[mm]))	
  bot   <- top.1 + top.2 + 
           pi.3[mm] * dnorm(lengths.sub, mean.0.yr[mm], sqrt(var.0.yr[mm])) 
  
  # Calculate the tau scores for each group
      
  taus.1 <- top.1 / bot	
  taus.2 <- top.2 / bot																													
  
  # Return the pis for each group	  	
	    
  pi.1.return <- sum(taus.1) / length(taus.1)										
  pi.2.return <- sum(taus.2) / length(taus.2)	
  pi.3.return <- 1 - (pi.1.return + pi.2.return)	
  
  # Return as a concatenated set
  
  return(c(pi.1.return,  pi.2.return, pi.3.return))
}