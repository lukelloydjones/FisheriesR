# Mean function for calcuating the mean length for the current month
# ------------------------------------------------------------------


MeanLength <- function(mm, k0, theta.1, theta.2 , linf, mu.yr.1, mu.yr.2, 
			  yrs.old, str.mnth) {
	
  # Computes the mean length of the distribution for the current month
  # based on an integration over the seasonal curve from a von Bertal
  # anffy growth model.
  #
  # Args:
  #  mm: 	 	 Current month to calulate mean at
  #  k0:       Mean  k paramter for VB model
  #  theta.1:  Seasonality parameter 1
  #  theta.2:  Seasonality parameter 1
  #  linf:     Asymptotic length
  #  mu.yr.1:  First month's average length yr 1
  #  mu.yr.2:  First month's average length yr 2
  #  yrs.old:  How old the individuals in this group are
  #  str.mnth: The month that we start calculating means for 
  #
  # Returns:
  #  Mean length for the month in the current iteration
  
  
  # Set each of the months to be the middle of the month (1/24 to be used
  # with the seasonal function.
  
  mm.val  <- (mm %% 12) / 12 + 1 / 24 + yrs.old
  str.mid <- str.mnth / 12 + 1 / 24 + yrs.old
  end.mid <- str.mid + 1


  # Assess whether the parameters at this update cross the y=0 axis
  # This will allow us to assess whether we need to calculate roots 
  # for the seasonal function or not
  
  time      <- seq(0, 1, 0.01)
  seas.func <- k0 + theta.1 * cos( 2 * pi * time) + T1 * sin( 2 * pi * time)
  is.neg.1  <- min(seas.func)
  is.neg.2  <- max(seas.func)
	  
	  if (Isneg < 0 & Isneg2 > 0)
	  {
	    a=T0^2+T1^2
		b=2*K0*T0
		c=-(T1^2)+K0^2
		u1=(-b+sqrt(b^2-4*a*c))/(2*a)
		u2=(-b-sqrt(b^2-4*a*c))/(2*a)
		r11=1-acos(u1)/(2*pi)+yrsold
		r12=acos(u1)/(2*pi)+yrsold
		r21=1-acos(u2)/(2*pi)+yrsold
		r22=acos(u2)/(2*pi)+yrsold
		roots<-c(r11, r12, r21, r22)
		groot<-K0+T0*cos(2*pi* roots)+T1*sin(2*pi* roots)
		gmin<-round(groot)
		r1<-min(roots[which(gmin==0)])
		r2<-max(roots[which(gmin==0)])
	    
	    
	    
	    
	    #Regime
	    
		zMid1=K0*(r1-strmid) + (T0/(2*pi))*(sin(2*pi*r1)-sin(2*pi*strmid)) - (T1/(2*pi))*(cos(2*pi*r1)-cos(2*pi*strmid))
	    zEnd1=K0*(endmid-r2) + (T0/(2*pi))*(sin(2*pi*endmid)-sin(2*pi*r2)) - (T1/(2*pi))*(cos(2*pi*endmid)-cos(2*pi*r2))
	    KKyr<-yrsold*(zEnd1+zMid1)
	    
	    #Those less than r1
	    
	    if (mm2<r1) {KK = KKyr + K0*(mm2-strmid) +   (T0/(2*pi))*(sin(2*pi*mm2)-sin(2*pi*strmid))-(T1/(2*pi))*(cos(2*pi*mm2)-cos(2*pi*strmid))}
	    	
	    #Those between r1 and r2
	    
	    if (mm2>r1&mm2<r2) {KK= KKyr + K0*(r1-strmid) + (T0/(2*pi))*(sin(2*pi*r1)-sin(2*pi*strmid))-(T1/(2*pi))*(cos(2*pi*r1)-cos(2*pi*strmid))}
	    
	    #Those after r2
	    
	    if (mm2>r2)
	    {
	    zMid=K0*(r1-strmid) + (T0/(2*pi))*(sin(2*pi*r1)-sin(2*pi*strmid)) - (T1/(2*pi))*(cos(2*pi*r1)-cos(2*pi*strmid))
	    zEnd=K0*(mm2-r2)     + (T0/(2*pi))*(sin(2*pi*mm2)-sin(2*pi*r2))     - (T1/(2*pi))*(cos(2*pi*mm2)-cos(2*pi*r2))
	    KK<-zEnd+zMid+KKyr
	    }
	    
	  } else {KK=K0*(mm2-strmid)+(T0/(2*pi))*(sin(2*pi*mm2)-sin(2*pi*strmid))-(T1/(2*pi))*(cos(2*pi*mm2)-cos(2*pi*strmid)) + yrsold*(K0*(endmid-strmid)+(T0/(2*pi))*(sin(2*pi*endmid)-sin(2*pi*strmid))-(T1/(2*pi))*(cos(2*pi*endmid)-cos(2*pi*strmid)) )}
	
	 
	  if (MMLIST[mm]>max(MM1985)) {M02 + (LINF-M02)*(1-exp(-KK))} else (M01 + (LINF-M01)*(1-exp(-KK)))
	
	
	}
