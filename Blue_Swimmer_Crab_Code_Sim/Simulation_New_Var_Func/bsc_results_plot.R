##############################################################################
##############################################################################
##############################################################################
###																		   ###
###	     SCRIPT TO DRAW PLOTS FROM CRAB SIMULATION RESULTS				   ###
###																		   ###	
##############################################################################
##############################################################################
##############################################################################

# Need to generate some data for the visualisation to work below

# Plots to visualise the results from the bsc simulation study

# Set the par space and names for each of the windows
	  
par(mfrow = c(3, 5))
XX <- 0:210
xlabNam <- c("Jan.", "Feb.", "Mar.", "Apr.", "May", "Jun.", "Jul.", "Aug.",
		 "Sep.", "Oct.", "Nov.", "Dec.")
K0   <- 1.03
T0   <- 1.86
T1   <- 2.01
LINF <- 189
M01  <- 75.5
VV   <- c(42.1,  0.0301, 3.90,  987)
PI1  <- c(0.11, 0.17, 0.12, 0.12, 0.15, 0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.15)
PI2  <- c(0.39, 0.43, 0.58, 0.68, 0.75, 0.79, 0.82, 0.82, 0.82, 0.82, 0.82, 0.83)


# Plot the monthly data as a window with different colours for each of the 
# juvenile, 1 yr adults, and 2 yr adults

for (mm in 1:MON)
 {
 #mm=1
  hist(LL[which(MM==MMLIST[mm])], breaks = 30, prob = T, xlim = c(20, 210), 
      ylim = c(0, 0.055), xlab = xlabNam[mm], ylab = "", main = '')
  lines(XX, PI1[mm] * dnorm(XX, MU2A(mm, K0, T0, T1, LINF, M01), 
  	  sqrt(VAR2A(mm, VV, K0, T0, T1, LINF, M01))), col =  'red',  lwd = 2.5)
  lines(XX, PI2[mm] * dnorm(XX,  MUA(mm, K0, T0, T1, LINF, M01), 
      sqrt(VARA(mm, VV, K0, T0, T1, LINF, M01))), col = 'green',  lwd = 2.5)
  lines(XX,   (1 - (PI1[mm] + PI2[mm])) * dnorm(XX, MUJ(mm, K0, T0, T1, LINF, M01), 
      sqrt(VARJ(mm, VV, K0, T0, T1, LINF, M01))), col =' blue',   lwd = 2.5)
  }


# Plot of seasonal curve

XX2 <- seq(0, 1, length.out = 100)
YY  <- K0 + T0 * cos(2 * pi * XX2) + T1 * sin(2 * pi * XX2)
YY_TRUE <-1 + 2 * cos(2 * pi * XX2) + 2 * sin(2 * pi * XX2)
plot(XX2, ((YY + abs(YY)) / 2), type = 'l', xlab = "Fraction of a year since Jan 1", 
	ylab = "k (per year) ", lwd = 2.5)
lines(XX2, ((YY_TRUE + abs(YY_TRUE)) / 2), col = 'blue', lwd =2.5)


# Plot of variance function

BLAH  <- 0:210
PARA2 <- c(40.0,  0.0300,  4.00,  900)
var_linf <- PARA2[1] * BLAH * exp(-PARA2[2] * BLAH) + exp(PARA2[3] * 
			(1 - exp(-PARA2[4] * BLAH)))
PARA2_est <- c(42.1,  0.0301, 3.90,  987)
var_linf_true <- PARA2_est[1] * BLAH * exp(-PARA2_est[2] * BLAH) + exp(PARA2_est[3] 
				 * (1 - exp(-PARA2_est[4] * BLAH)))
PARA2_true[1] * 189.4 * exp(-PARA2_true[2] * 189.4) + exp(PARA2_true[3] * 
(1 - exp(-PARA2_true[4] * 189.4))) 
plot(BLAH, var_linf,type = 'l', xlab = "Mean length (mm)", ylab = "Variance", 
	lwd = 2.5)
lines(BLAH, var_linf_true, col = "blue", lwd = 2.5)


# Plot of asymptotic length density

plot(density(rnorm(1000000, LINF, sqrt(VV[1] * LINF * exp(-VV[2] * LINF) +
    exp(VV[3] * (1 - exp(-VV[4] * LINF)))))), xlab = "Length (mm)", main = "", 
    lwd = 2.5, col = "blue")
VV_tru <- c(40.0,  0.0300,  4.00,  900)
lines(density(rnorm(1000000, 190, sqrt(VV_tru[1] * 190 * exp(-VV_tru[2] * 190) +
    exp(VV_tru[3] * (1 - exp(-VV_tru[4] * 190)))))), lwd = 2.5)










