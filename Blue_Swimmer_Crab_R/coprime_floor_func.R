CoprimeFloor <- function(a, b)
{
	#install.packages("pracma")
	#library(pracma)
	
	stopifnot(gcd(a, b) == 1 )
	array.1 <- array(0, b-1) 
	for (i in seq(1, b-1))
	{
		array.1[i] <- floor(i*a/b)
	}
	
	
	array.2 <- array(0, a-1) 
	for (j in seq(1, a-1))
	{
		array.2[j] <- floor(j*b/a)
	}
	
	plot(array.1, cumsum(array.1), type = "l", col = "blue", 
	     xlim = c(0, max(c(array.1, array.2))), 
	     ylim = c(0, max(c(cumsum(array.1),cumsum(array.2)))) )
	lines(array.2, cumsum(array.2), type = "l", col = "red")
	ret <- list(array.1, array.2)
	
	return(ret)
}