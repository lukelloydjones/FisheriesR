######################################################################################
######################################################################################
######################################################################################
###																				   ###
###			     SCRIPT TO DRAW PLOTS FROM CRAB SIMULATION RESULTS				   ###
###																				   ###	
######################################################################################
######################################################################################
######################################################################################


# Plots to visualise the results from the bsc simulation study

# Set the par space and names for each of the windows
	  
par(mfrow = c(3,5))
XX <- 0:200
xlabNam<-c("Jan.","Feb.", "Mar.","Apr.","May","Jun.","Jul.","Aug.","Sep.","Oct.","Nov.","Dec.")

# Plot the monthly data as a window with different colours for each of the juvenile, 1 yr adults, and 2 yr adults

for (mm in 1:MON)
 {
  hist(LL[which(MM==MMLIST[mm])],breaks=30,prob=T,xlim=c(20,210),ylim=c(0,0.055),xlab=xlabNam[mm],ylab="",main='')
  lines(XX,PI1[mm]*dnorm(XX,MU2A(mm,K0,T0,T1,LINF,M01),sqrt(VAR2A(mm,VV,K0,T0,T1,LINF,M01))),col='red',lwd=2.5)
  lines(XX,PI2[mm]*dnorm(XX,MUA(mm,K0,T0,T1,LINF,M01),sqrt(VARA(mm,VV,K0,T0,T1,LINF,M01))),col='green',lwd=2.5)
  lines(XX,(1-(PI1[mm]+PI2[mm]))*dnorm(XX,MUJ_AUX[mm],sqrt(VARJ(mm,VV,K0,T0,T1,LINF,M01))),col='blue',lwd=2.5)
  }
  
	  XX2 <- seq(0,1,length.out=100)
	  YY <- K0 + T0*cos(2*pi*XX2) + T1*sin(2*pi*XX2)
	  YY_TRUE<-1 + 2*cos(2*pi*XX2) + 2*sin(2*pi*XX2)
	  plot(XX2,((YY+abs(YY))/2),type='l',xlab="Fraction of a year since Jan 1",ylab="k (per year)",lwd=2.5)
	  lines(XX2, ((YY_TRUE +abs(YY_TRUE))/2), col='blue',lwd=2.5)
	  #VV=c(3,150,100)
	  BLAH <- 0:200
	  var_linf= PARA2[1]*BLAH*exp(-PARA2[2]*BLAH)+exp(PARA2[3]*(1-exp(-PARA2[4]*BLAH)))
	  PARA2_true<-c(30,0.03,-20,900)
	  var_linf_true = PARA2_true[1]*BLAH*exp(-PARA2_true[2]*BLAH)+exp(PARA2_true[3]*(1-exp(-PARA2_true[4]*BLAH)))
	  plot(BLAH,var_linf,type='l',xlab="Mean length (mm)", ylab="Variance",lwd=2.5)
	  lines(BLAH,var_linf_true,col="blue",lwd=2.5)
	  plot(density(rnorm(1000000,PARA[4],sqrt(VV[1]*PARA[4]*exp(-VV[2]*PARA[4])+exp(VV[3]*(1-exp(-VV[4]*PARA[4])))))),xlab="Length (mm)",main="",lwd=2.5)