######################################################################################
######################################################################################
######################################################################################
###																				   ###
###			     SCRIPT TO RUN SIMULATIN STUDY FOR BLUE SWIMMER CRABS			   ###
###																				   ###	
######################################################################################
######################################################################################
######################################################################################

# This script generates simulated data that mimics that of length frequency data gathered on the 
# blue swimmer crab and provides a testing ground for the algorithm written for estimating growth
# via a mixture model

# Remove any objects to clear the slate

rm(list=ls())

# First set the directory to write out the results

#setwd("/ibscratch/wrayvisscher/Luke/BSC_RcodesData/Simulation200")


# Generate three mean curves via von Bertalanffy
# Generate three mean curves via von Bertalanffy


# Set some initial paramters that are to be estimated in the next section. We set the seasonal function so that is passes below 
# zero. This allows us to test the algorithm for the ability to deal with this scenario.

days<-seq(0,1,1/366)
IL<-60
K0<-1
T0<-2
T1<-2
LINF<-190


# It is difficult to deal with the space between the mid section of each month and the start of the next year and thus
# we calculate the means over a fine mesh of values and take the average over the months as a proxy for the mean length
# for that month. This allows us to calculate the means over a whole year and then just add 1 and 2 to the integral
# values for the 1 year and 2 year old adults.

KK_store<-array(0,length(days)) # Initialise the array to store the mean values

# Cycle through each day and calculate the integral over the seasonal function

for (i in seq(1,length(days)))
	{
	
	# Calculate if the curve needs to be integrated or not and set the initial month which will be January for the simulation
	
  	mm2<-days[i]
  	yrsold<-0
  	strmnth<-0
	strmid<-0

	# Calculate the roots for the seasonal function. This allows us to integrate up to the root and set the integral between 
	# roots to 0. 
	
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
	groot<-K0+T0*cos(2*pi* roots)+T1*sin(2*pi*roots)
	gmin<-round(groot)
		
	# Set the roots
		
	r1<-min(roots[which(gmin==0)])
	r2<-max(roots[which(gmin==0)])
	    
    # The integral of those months less than root 1
	    
	if (mm2<r1) {KK = K0*(mm2-strmid) +   (T0/(2*pi))*(sin(2*pi*mm2)-sin(2*pi*strmid))-(T1/(2*pi))*(cos(2*pi*mm2)-								cos(2*pi*strmid))}
	    	
	# The integral of those months between root 1 and root 2. This will just be the integral up to root 1.
	    
	if (mm2>r1&mm2<r2) {KK=K0*(r1-strmid) + (T0/(2*pi))*(sin(2*pi*r1)-sin(2*pi*strmid))-(T1/(2*pi))*(cos(2*pi*r1)-										cos(2*pi*strmid))}
	    
    # The integral of those months after root 2
	    
	if (mm2>r2)
	   {
	   	zMid=K0*(r1-strmid) + (T0/(2*pi))*(sin(2*pi*r1)-sin(2*pi*strmid)) - (T1/(2*pi))*(cos(2*pi*r1)-cos(2*pi*strmid))
	   	zEnd=K0*(mm2-r2)    + (T0/(2*pi))*(sin(2*pi*mm2)-sin(2*pi*r2))     - (T1/(2*pi))*(cos(2*pi*mm2)-cos(2*pi*r2))
	   	KK<-zEnd+zMid
	   }
	
	# Print it out just to check
	
	print(KK)
	
	# Store each of the value so that we can calculate the means
	
	KK_store[i]<-KK
	
}


# Define the integrals over each year. Subsequent years are just multiples of the first year because the intergal is the same but 
# in multiples for each subsequent year.

KK_yr1=KK_store[1:365]
KK_yr2=KK_store[1:365]+1.44099
KK_yr3=KK_store[1:365]+2*1.44099


# Calculate the mean for each day over the year

MLY1<-(IL + (LINF-IL)*(1-exp(-KK_yr1)))
MLY2<-(IL + (LINF-IL)*(1-exp(-KK_yr2)))
MLY3<-(IL + (LINF-IL)*(1-exp(-KK_yr3)))


# Put each day into a set of months which is designated to be 30 days (only an approximation to true value). We then take the 
# mean over each of the columns to obtain the monthly mean

MLY1 <- colMeans(matrix(MLY1,ncol=12,nrow=30))
MLY2 <- colMeans(matrix(MLY2,ncol=12,nrow=30))
MLY3 <- colMeans(matrix(MLY3,ncol=12,nrow=30))


# Round these means for summary in the manuscript

paste(as.character(signif(MLY1,4)),collapse=" & ")
paste(as.character(signif(MLY2,4)),collapse=" & ")
paste(as.character(signif(MLY3,4)),collapse=" & ")


## Draw from the means to create the groups. Groups are drawn from a normal distribution to mimic our assumptions in real data
## Draw from the means to create the groups. Groups are drawn from a normal distribution to mimic our assumptions in real data


# Need to vary the number in each group in each month to better mimic reality
# Set up the initial parameters for the variance function


# Variance function test to see shape

VV<-c(40,0.03,4,900)
mu<-seq(0,200,0.5)
y=VV[1]*mu*exp(-VV[2]*mu)+exp(VV[3]*(1-exp(-VV[4]*mu)))
plot(mu,y)


# Set the proportions in each group varying them for each month

P1<-c(0.50, 0.40, 0.30, 0.20, 0.10, 0.05, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02)
P3<-c(0.1, 0.15, 0.1, 0.08, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09)
P2<-1-(P1+P3)


# Calculate the number of individuals in each group based on the proportions

N=2000 # Number of total individuals
N1<-as.integer(N*P1)
N2<-as.integer(N*P2)
N3<-N-(N1+N2)


# Make the empty arrays for the individuals to go in

G<-matrix(0,nrow=12,ncol=N)

# For each of the twelve months create the three distributions

for (i in seq(1,12))
{

# Calculate the variances

var1<-VV[1]*MLY1[i]*exp(-VV[2]*MLY1[i])+exp(VV[3]*(1-exp(-VV[4]*MLY1[i])))
var2<-VV[1]*MLY2[i]*exp(-VV[2]*MLY2[i])+exp(VV[3]*(1-exp(-VV[4]*MLY2[i])))
var3<-VV[1]*MLY3[i]*exp(-VV[2]*MLY3[i])+exp(VV[3]*(1-exp(-VV[4]*MLY3[i])))

# Draw the individuals from each of the distributions

G1<-rnorm(N1[i],MLY1[i], sqrt(var1))
G2<-rnorm(N2[i],MLY2[i], sqrt(var2))
G3<-rnorm(N3[i],MLY3[i], sqrt(var3))

G[i,]<-c(G1,G2,G3)
}



### RUN THE CODE AND GET THE ESTIMATES
### RUN THE CODE AND GET THE ESTIMATES
### RUN THE CODE AND GET THE ESTIMATES

# Given that the data are generated we attempt to estimate the true parameters using our model
# Given that the data are generated we attempt to estimate the true parameters using our model


MON <- 12
MONTHS<-rep(seq(0,11),each=N)
MM<-cbind(matrix(MONTHS,ncol=12,nrow=N))
MM<-t(MM)
MM <-as.numeric(MM)
LL<-as.numeric(G)

# Initialise the model parameters

NN <- length(MM)
PI1 <- P3
PI2 <- P2
K0 <- 1
LINF <- 185
M01<-40
T0 <- 2
T1 <- 2
MMLIST <- seq(0,11)

#Initialise the starting values for the variance quadratic update

VV <- c(5,1/100,3,1)										

## FUNCTIONS
## FUNCTIONS

#ASYMPTOTIC MEAN FUNCTION

MU2A <- function(mm,K0,T0,T1,LINF,M01)
{
  #mm=12
  mmval <- MMLIST[mm]%%12
  yrsold<-2
  strmnth<-0
  strmid<-strmnth/12+1/24+yrsold
  endmid<-strmid+1
  mm2=(mmval/12)+1/24+yrsold
  t<-seq(0,1,0.01)
  g<-K0+T0*cos(2*pi*t)+T1*sin(2*pi*t)
  #plot(t,g)
  Isneg<-min(g)
  Isneg2<-max(g)
  if (Isneg<0&Isneg2>0)
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
    
    if (mm2<r1) {KK = KKyr + K0*(mm2-strmid) +   (T0/(2*pi))*(sin(2*pi*mm2)-sin(2*pi*strmid))-(T1/(2*pi))*(cos(2*pi*mm2)-						cos(2*pi*strmid))}
    	
    #Those between r1 and r2
    
    if (mm2>r1&mm2<r2) {KK= KKyr + K0*(r1-strmid) + (T0/(2*pi))*(sin(2*pi*r1)-sin(2*pi*strmid))-(T1/(2*pi))*(cos(2*pi*r1)-								cos(2*pi*strmid))}
    
    #Those after r2
    
    if (mm2>r2)
    {
    zMid=K0*(r1-strmid) + (T0/(2*pi))*(sin(2*pi*r1)-sin(2*pi*strmid)) - (T1/(2*pi))*(cos(2*pi*r1)-cos(2*pi*strmid))
    zEnd=K0*(mm2-r2)     + (T0/(2*pi))*(sin(2*pi*mm2)-sin(2*pi*r2))     - (T1/(2*pi))*(cos(2*pi*mm2)-cos(2*pi*r2))
    KK<-zEnd+zMid+KKyr
    }
    
  } else 
  {KK=K0*(mm2-strmid)+(T0/(2*pi))*(sin(2*pi*mm2)-sin(2*pi*strmid))-(T1/(2*pi))*(cos(2*pi*mm2)-cos(2*pi*strmid)) + 				                           	yrsold*(K0*(endmid-strmid)+(T0/(2*pi))*(sin(2*pi*endmid)-sin(2*pi*strmid))-(T1/(2*pi))*(cos(2*pi*endmid)-cos(2*pi*strmid)) )}

 
  M01 + (LINF-M01)*(1-exp(-KK))


}


#1 YEAR OLD ADULT's MEAN FUNCTION


MUA <- function(mm,K0,T0,T1,LINF,M01)
{
  mmval <- MMLIST[mm]%%12
  yrsold<-1
  strmnth<-0
  strmid<-strmnth/12+1/24+yrsold
  endmid<-strmid+1
  mm2=(mmval/12)+1/24+yrsold
  t<-seq(0,1,0.01)
  g<-K0+T0*cos(2*pi*t)+T1*sin(2*pi*t)
  Isneg<-min(g)
  Isneg2<-max(g)
  if (Isneg<0&Isneg2>0)
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
    
    	if (mm2<r1) {KK = KKyr + K0*(mm2-strmid) +   (T0/(2*pi))*(sin(2*pi*mm2)-sin(2*pi*strmid))-(T1/(2*pi))*(cos(2*pi*mm2) -						cos(2*pi*strmid))}
    	
    	#Those between r1 and r2
    
    	if (mm2>r1&mm2<r2) {KK= KKyr + K0*(r1-strmid) + (T0/(2*pi))*(sin(2*pi*r1)-sin(2*pi*strmid))-(T1/(2*pi))*(cos(2*pi*r1) -								cos(2*pi*strmid))}
    
    	#Those after r2
    
    	if (mm2>r2) {
    		zMid <- K0*(r1-strmid) + (T0/(2*pi))*(sin(2*pi*r1)-sin(2*pi*strmid)) - (T1/(2*pi))*(cos(2*pi*r1)-cos(2*pi*strmid))
    		zEnd <- K0*(mm2-r2)     + (T0/(2*pi))*(sin(2*pi*mm2)-sin(2*pi*r2))     - (T1/(2*pi))*(cos(2*pi*mm2)-cos(2*pi*r2))
    		KK   <- zEnd+zMid+KKyr
    		} }else {
    		KK  <- K0*(mm2-strmid)+(T0/(2*pi))*(sin(2*pi*mm2)-sin(2*pi*strmid))-(T1/(2*pi))*(cos(2*pi*mm2)-cos(2*pi*strmid)) + 					yrsold*(K0*(endmid-strmid)+(T0/(2*pi))*(sin(2*pi*endmid)-sin(2*pi*strmid))-(T1/(2*pi))*(cos(2*pi*endmid) -						cos(2*pi*strmid)) )}

 	
 	
  (M01 + (LINF-M01)*(1-exp(-KK)))


}


#THE JUVENUILES MEAN FUNCTION

MUJ <- function(mm,K0,T0,T1,LINF,M01)
{
	#mm=2
	mmval <- MMLIST[mm]%%12
  	yrsold<-0
  	strmnth<-0
  	strmid<-strmnth/12+1/24+yrsold
  	endmid<-strmid+1
  	mm2=(mmval/12)+1/24+yrsold
  	t<-seq(0,1,0.01)
	g<-K0+T0*cos(2*pi*t)+T1*sin(2*pi*t)
	Isneg<-min(g)
	Isneg2<-max(g)
 	if (Isneg<0&Isneg2>0)
  		{
    	a=T0^2+T1^2
		b=2*K0*T0
		c=-(T1^2)+K0^2
		u1=(-b+sqrt(b^2-4*a*c))/(2*a)
		u2=(-b-sqrt(b^2-4*a*c))/(2*a)
		r11=1-acos(u1)/(2*pi)
		r12=acos(u1)/(2*pi)
		r21=1-acos(u2)/(2*pi)
		r22=acos(u2)/(2*pi)
		roots<-c(r11, r12, r21, r22)
		groot<-K0+T0*cos(2*pi* roots)+T1*sin(2*pi* roots)
		gmin<-round(groot)
		r1<-min(roots[which(gmin==0)])
		r2<-max(roots[which(gmin==0)])
		
		
    	#Regime

    	#Those less than r1
    
    		if (mm2<r1) {KK = K0*(mm2-strmid) + (T0/(2*pi))*(sin(2*pi*mm2)-sin(2*pi*strmid))-(T1/(2*pi))*(cos(2*pi*mm2) -								 cos(2*pi*strmid))}
    	
    		#Those between r1 and r2
    
    		if (mm2>r1&mm2<r2) {KK=K0*(r1-strmid) + (T0/(2*pi))*(sin(2*pi*r1)-sin(2*pi*strmid))-(T1/(2*pi))*(cos(2*pi*r1) -										cos(2*pi*strmid))}
    
    		#Those after r2
    
    		if (mm2>r2)
    			{
    			zMid=K0*(r1-strmid) + (T0/(2*pi))*(sin(2*pi*r1)-sin(2*pi*strmid)) - (T1/(2*pi))*(cos(2*pi*r1)  -								cos(2*pi*strmid))
    			zEnd=K0*(mm2-r2)    + (T0/(2*pi))*(sin(2*pi*mm2)-sin(2*pi*r2))    - (T1/(2*pi))*(cos(2*pi*mm2) - 								cos(2*pi*r2))
    			KK<-zEnd+zMid
    			}} else {
    	KK=K0*(mm2-strmid)+(T0/(2*pi))*(sin(2*pi*mm2)-sin(2*pi*strmid))-(T1/(2*pi))*(cos(2*pi*mm2)-cos(2*pi*strmid))}	


 M01 + (LINF-M01)*(1-exp(-KK))
}


#INITIALISE THE TAU's which are the inclusion probablities to all be a third


TAU <- list()
for (mm in 1:(MON))
	{
  		TAU[[mm]] <- matrix(1/3,NN,3)
	}


#INITIALISE THE AUXILLARY STORAGE ARRAYS FOR INSIDE LOOP ASSIGNMENT

MU2A_AUX <- c()
MUA_AUX  <- c()
MUJ_AUX  <- c()
for (mm in 1:MON)
	{
 	MU2A_AUX[mm] <- MU2A(mm,K0,T0,T1,LINF,M01)
  	MUA_AUX[mm]  <- MUA(mm,K0,T0,T1,LINF,M01)
 	MUJ_AUX[mm]  <- MUJ(mm,K0,T0,T1,LINF,M01)
	}


#INITIALISE THE VARIANCE FUNCTIONS. THEY ARE ALL QUADRATIC FUNCTIONS OF THE MEAN FOR EACH COHORT

VARJ <- function(mm,VV,K0,T0,T1,LINF,M01)
	{
		
		#Try a heavy tailed function like the Weibull for more flexibility
		
		mu=MUJ(mm,K0,T0,T1,LINF,M01)
		max(VV[1]*mu*exp(-VV[2]*mu)+exp(VV[3]*(1-exp(-VV[4]*mu))),1)
	}

VARA <- function(mm,VV,K0,T0,T1,LINF,M01)
	{
		
		#Try a heavy tailed function like the Weibull for more flexibility
		
		mu=MUA(mm,K0,T0,T1,LINF,M01)
		max(VV[1]*mu*exp(-VV[2]*mu)+exp(VV[3]*(1-exp(-VV[4]*mu))),1)
	}

VAR2A <- function(mm,VV,K0,T0,T1,LINF,M01)
	{
		
		#Try a heavy tailed function like the Weibull for more flexibility
		
		mu <- MU2A(mm,K0,T0,T1,LINF,M01)
		max(VV[1]*mu*exp(-VV[2]*mu)+exp(VV[3]*(1-exp(-VV[4]*mu))),1)
	}



#INITIALISE THE LIKELIHOOD FUNCTION

LOGLIKE <- 0
for (mm in 1:MON)
  {
	LL2 <- LL[which(MM==MMLIST[mm])]
	LikGrp1 <- PI1[mm]*dnorm(LL2,MU2A_AUX[mm],sqrt(VAR2A(mm,VV,K0,T0,T1,LINF,M01)))
	LikGrp2 <- PI2[mm]*dnorm(LL2,MUA_AUX[mm],sqrt(VARA(mm,VV,K0,T0,T1,LINF,M01)))
	LikGrp3 <- (1 - (PI1[mm] + PI2[mm]))*dnorm(LL2,MUJ_AUX[mm],sqrt(VARJ(mm,VV,K0,T0,T1,LINF,M01)))
	LOGLIKE <- LOGLIKE + sum(log(LikGrp1 + LikGrp2 + LikGrp3))
  }



	
	
	
######WHILE LOOP BEGINS######WHILE LOOP BEGINS######WHILE LOOP BEGINS
######WHILE LOOP BEGINS######WHILE LOOP BEGINS######WHILE LOOP BEGINS
######WHILE LOOP BEGINS######WHILE LOOP BEGINS######WHILE LOOP BEGINS
######WHILE LOOP BEGINS######WHILE LOOP BEGINS######WHILE LOOP BEGINS
	
	
#THIS LOOP DOES ALL THE WORK AND INCLUDES TWO NELDER MEAD STEPS TO OPTIMISE THE NONLINEAR MEANS AND VARIANCES


TOL <- 10^-6
LOGOLD <- -Inf

while (LOGLIKE - LOGOLD > TOL)
{
  LOGOLD <- LOGLIKE      #Assign the current likelihood value to the an old value so we can evaluate the update

 
	# UPDATE THE TAU SCORES AND THEN CALCULATE THE PIs GIVEN THESE TAU SCORES

	for (mm in 1:MON)
  	{
    	
    	# Calculate the inclusion probablities for group 1 and 2 which forms the top of the mixture model algorithm

    	Top1<-PI1[mm]*dnorm(LL,MU2A_AUX[mm],sqrt(VAR2A(mm,VV,K0,T0,T1,LINF,M01))) 
    	Top2<-PI2[mm]*dnorm(LL,MUA_AUX[mm],sqrt(VARA(mm,VV,K0,T0,T1,LINF,M01)))								
    	
    	# Calculate the sum of each of the classes which forms the bottom of the mixture model algorithm

    	Bot<- Top1 + Top2 + (1-(PI1[mm]+PI2[mm]))*dnorm(LL,MUJ_AUX[mm],sqrt(VARJ(mm,VV,K0,T0,T1,LINF,M01)))  
    	
    	
    	# Calculate the tau scores for each of the groups, where the third group is what's left over from 1 after the first 2

    	TAU[[mm]][,1] <- Top1/Bot																	
    	TAU[[mm]][,2] <- Top2/Bot
		TAU[[mm]][,3] <- 1-(TAU[[mm]][,1]+TAU[[mm]][,2])
		
		
		# Work out the PIs for groups 1 and 2
												  		    
    	PI1[mm] <- sum((MM==MMLIST[mm])*TAU[[mm]][,1])/sum(MM==MMLIST[mm])					
		PI2[mm] <- sum((MM==MMLIST[mm])*TAU[[mm]][,2])/sum(MM==MMLIST[mm])										
	}


#DEFINE THE LIKELIHOOD FUNCTION OUTSIDE THE LOOP TO UPDATE THE MEANS

  PARA <- c(K0,T0,T1,LINF,M01)
  OPTIFUN <- function(PARA)
  {
  	LOGLIKE <- 0
  		for (mm in 1:MON)
  			{
			LL2 <-LL[which(MM==MMLIST[mm])]
			LikGrp1 <- PI1[mm]*dnorm(LL2,MU2A(mm,PARA[1],PARA[2],PARA[3],PARA[4],PARA[5]),
			sqrt(VAR2A(mm,VV,PARA[1],PARA[2],PARA[3],PARA[4],PARA[5])))
			
			LikGrp2 <-PI2[mm]*dnorm(LL2,MUA(mm,PARA[1],PARA[2],PARA[3],PARA[4],PARA[5]),
			sqrt(VARA(mm,VV,PARA[1],PARA[2],PARA[3],PARA[4],PARA[5]))) 
			
			LikGrp3 <- (1 - (PI1[mm] + PI2[mm]))*dnorm(LL2,MUJ(mm,PARA[1],PARA[2],PARA[3],PARA[4],PARA[5]),
			sqrt(VARJ(mm,VV,PARA[1],PARA[2],PARA[3],PARA[4],PARA[5])))
			LOGLIKE <- LOGLIKE + sum(log(LikGrp1 + LikGrp2 + LikGrp3))
  			}
	-LOGLIKE
  }
  
  
  # Run OPTIM to find the best variances given the taus and Pis
  
  OPTIM <- optim(PARA,OPTIFUN, control = list(maxit = 10000))		
  print("Did optim 1 converge?")
  print(OPTIM$convergence)
  
  
  #Assign the estimates from optim to the parameters
  
  K0 <- OPTIM$par[1]
  T0 <- OPTIM$par[2]
  T1 <- OPTIM$par[3]
  LINF <- OPTIM$par[4]
  M01 <- OPTIM$par[5]
  

  #Calculate the means for each group based on these new estimates
  
  MU2A_AUX <- c()
  MUA_AUX <- c()
  MUJ_AUX <- c()
  for (mm in 1:MON)
  {
  	#Compute the means for each month 
    MU2A_AUX[mm] <- MU2A(mm,K0,T0,T1,LINF,M01)
    MUA_AUX[mm] <- MUA(mm,K0,T0,T1,LINF,M01)
    MUJ_AUX[mm] <- MUJ(mm,K0,T0,T1,LINF,M01)
  }
  

  #OPTIMISE THE VARIANCES GIVEN THE NEW MEANS
  
  #DEFINE THE LIKELIHOOD FUNCTION OUTSIDE THE LOOP TO UPDATE THE VARIANCES 
  
  PARA2 <- VV
  VAROPTIFUN <- function(PARA2)
  {
  LOGLIKE <- 0
  		for (mm in 1:MON)
  			{
			LL2<-LL[which(MM==MMLIST[mm])]
			LikGrp1<- PI1[mm]*dnorm(LL2,MU2A_AUX[mm],sqrt(VAR2A(mm,PARA2,K0,T0,T1,LINF,M01)))
			LikGrp2<- PI2[mm]*dnorm(LL2,MUA_AUX[mm],sqrt(VARA(mm,PARA2,K0,T0,T1,LINF,M01)))
			LikGrp3<-(1-(PI1[mm]+PI2[mm]))*dnorm(LL2,MUJ_AUX[mm],sqrt(VARJ(mm,PARA2,K0,T0,T1,LINF,M01)))
			LOGLIKE <- LOGLIKE + sum(log(LikGrp1 + LikGrp2 + LikGrp3))
  			}
	-LOGLIKE
  }
	
		
  # Run OPTIM to find the best variances given the taus and Pis

  OPTIM2 <- optim(PARA2,VAROPTIFUN,control = list(maxit = 10000))													
  VV <- OPTIM2$par
  print("Did optim 2 converge?")
  print(OPTIM2$convergence)  
  
  
  #COMPUTE THE LIKELIHOOD GIVEN ALL THESE NICE NEW UPDATES
      
  LOGLIKE <- 0
  for (mm in 1:MON)
  {
	LL2<-LL[which(MM==MMLIST[mm])]
	LikGrp1<-PI1[mm]*dnorm(LL2,MU2A_AUX[mm],sqrt(VAR2A(mm,VV,K0,T0,T1,LINF,M01)))
	LikGrp2<-PI2[mm]*dnorm(LL2,MUA_AUX[mm],sqrt(VARA(mm,VV,K0,T0,T1,LINF,M01)))
	LikGrp3<-(1-(PI1[mm]+PI2[mm]))*dnorm(LL2,MUJ_AUX[mm],sqrt(VARJ(mm,VV,K0,T0,T1,LINF,M01)))
	LOGLIKE <- LOGLIKE + sum(log(LikGrp1 + LikGrp2 + LikGrp3))
  }

  
  #PRINT OUT OUR NEW UPDATES FOR THE PARAMETERS AND THE LIKELIHOOD
  
  print(c(LOGLIKE,LOGLIKE-LOGOLD))
  print(c(PARA,PARA2))
 
}

# Write all results to a text file for later reference

name<-paste("Sim_Res",sample(seq(1,100000),1),sep="")
write.table(matrix(c(PARA,PARA2),nrow=1,ncol=length(c(PARA,PARA2))),name,append=T,row.names=F, col.names=F,quote=F)
write.table(matrix(PI1,nrow=1,ncol=length(PI1)),name,append=T,row.names=F,col.names=F,quote=F)
write.table(matrix(PI2,nrow=1,ncol=length(PI2)),name,append=T,row.names=F,col.names=F,quote=F)
write.table(matrix(1-(PI1+PI2),nrow=1,ncol=length(PI2)),name,append=T,row.names=F,col.names=F,quote=F)
write.table(matrix(MU2A_AUX,nrow=1,ncol=length(MU2A_AUX)),name,append=T,row.names=F,col.names=F,quote=F)
write.table(matrix(MUA_AUX,nrow=1,ncol=length(MUA_AUX)),name,append=T,row.names=F,col.names=F,quote=F)
write.table(matrix(MUJ_AUX,nrow=1,ncol=length(MUJ_AUX)),name,append=T,row.names=F,col.names=F,quote=F)














