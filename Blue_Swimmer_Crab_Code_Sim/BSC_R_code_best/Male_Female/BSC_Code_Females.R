###################CODE TO ESTIMATE GROWTH PARAMETERS FOR FEMALE BLUE SWIMMER CRAB FROM LENGTH FREQUENCY DATA###########################
###################CODE TO ESTIMATE GROWTH PARAMETERS FOR FEMALE BLUE SWIMMER CRAB FROM LENGTH FREQUENCY DATA###########################
###################CODE TO ESTIMATE GROWTH PARAMETERS FOR FEMALE BLUE SWIMMER CRAB FROM LENGTH FREQUENCY DATA###########################
###################CODE TO ESTIMATE GROWTH PARAMETERS FOR FEMALE BLUE SWIMMER CRAB FROM LENGTH FREQUENCY DATA###########################


#Set the necessary working directories

#setwd("/Users/uqllloyd/Dropbox/AAUni/APhD/Blueswimmer/CrabStuff2withMMAlg/RcodesData")
setwd("~/Dropbox/AAUni/APhD/Blueswimmer/RcodeData_03:10:2014/All_Months")


#Remove any old objects

rm(list = ls())


#DATA PRELIMINARIES

#Read in the data set on the asymptotic males that was gathered through pots

LFD_bigMalesFemales<-read.table("LFD_bigMalesFem",header=T)


#Pull out the lengths and the dates from these data files

LFD_bigFem<-LFD_bigMalesFemales[which(LFD_bigMalesFemales$Sex==2),]
Dates1<-as.Date(LFD_bigFem$Date,"%d/%m/%y")               
Lengths1<-LFD_bigFem$Carapacewidth


#Read in the trawl data on males that contains juvenile recruitment and adults

LFD<-read.table("LFD")
LFD<-LFD[which(LFD$Sex==2),]


#Pull out the lengths and the dates from these data files

Dates2<-as.Date(LFD$Date,"%d/%m/%Y")
Lengths2<-LFD$Carapace.width


#Concatenate the necessary elements from each file into a common dates and lengths array
 
Dates <- c(Dates2,Dates1)
Lengths <- c(Lengths2,Lengths1)


#Pull out the year and month information from these dates

YEAR <- format(Dates,'%Y')
MONTHS <- format(Dates,'%m')


#Pull out the years and months that we a re interested in i.e., those that don't contain 	recruitment

WHICH85JO <- which((YEAR == '1985') & (as.numeric(MONTHS) %in% (2:8)))
WHICH86FM <- which((YEAR == '1986') & (as.numeric(MONTHS) %in% (2:5)))



#INITIALISE


MON <- 11													#Assign how many months we would like to model
LL <- Lengths[c(WHICH85JO,WHICH86FM)]						#Assign the number of individuals
MM1985 <- as.numeric(MONTHS[WHICH85JO])-1					#Assign January to be the 0th months
MM1986 <- as.numeric(MONTHS[WHICH86FM])+11					#Thus January of next year will be the 12 month. Plus we do this so that we don't have to estimate yr
MM<-c(MM1985,MM1986)										#Concatenate these months values together
MMLIST <- as.numeric(names(table(MM)))						#Make a list of the months from January first i.e., January 1986 is gets assigned a 12



#Parameter Initialise

NN <- length(MM)											#Initialise the number of individuals we have
PI1 <- rep(1/3,MON)											#Initialise the PIs
PI2 <- rep(1/3,MON)											#Initialise the PIs
K0 <- 3														#Initialise K0 average K
LINF <- 200													#Initialise asym length
M01<-40														#Initialise first month's average length
M02<-40														#Initialise second month's average length
T0 <- 1														#Initialise first seasonality parameter
T1 <- 1														#Initialise second seasonality parameter
VV <- c(100,0,0)										   #Initialise the starting values for the variance quadratic update


#Setting values for drawing plots 

# K0=0.84062808;  LINF=159.65387413; M01=58.86155838; M02=81.29860313; T0=-0.07925965; T1=0.57794716; VV=c( -4.77584147,   6.19436755  ,-0.03489903)
# #0.84062808  -0.07925965   0.57794716 159.65387413  58.86155838  81.29860313  -4.77584147   6.19436755  -0.03489903

#FUNCTIONS


#ASYMPTOTIC MEAN FUNCTION

MU2A <- function(mm,K0,T0,T1,LINF,M01,M02)
{
  #mm=1
  mmval <- MMLIST[mm]%%12
  yrsold<-2
  strmnth<-min(MM1985)
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


#1 YEAR OLD ADULT's MEAN FUNCTION


MUA <- function(mm,K0,T0,T1,LINF,M01,M02)
{
  mmval <- MMLIST[mm]%%12
  yrsold<-1
  strmnth<-min(MM1985)
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


#THE JUVENUILES MEAN FUNCTION

MUJ <- function(mm,K0,T0,T1,LINF,M01,M02)
{
	#mm=2
	mmval <- MMLIST[mm]%%12
  	yrsold<-0
  	strmnth<-min(MM1985)
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
    
    		if (mm2<r1) {KK = K0*(mm2-strmid) +   (T0/(2*pi))*(sin(2*pi*mm2)-sin(2*pi*strmid))-(T1/(2*pi))*(cos(2*pi*mm2)-cos(2*pi*strmid))}
    	
    		#Those between r1 and r2
    
    		if (mm2>r1&mm2<r2) {KK=K0*(r1-strmid) + (T0/(2*pi))*(sin(2*pi*r1)-sin(2*pi*strmid))-(T1/(2*pi))*(cos(2*pi*r1)-cos(2*pi*strmid))}
    
    		#Those after r2
    
    		if (mm2>r2)
    			{
    			zMid=K0*(r1-strmid) + (T0/(2*pi))*(sin(2*pi*r1)-sin(2*pi*strmid)) - (T1/(2*pi))*(cos(2*pi*r1)-cos(2*pi*strmid))
    			zEnd=K0*(mm2-r2)     + (T0/(2*pi))*(sin(2*pi*mm2)-sin(2*pi*r2))     - (T1/(2*pi))*(cos(2*pi*mm2)-cos(2*pi*r2))
    			KK<-zEnd+zMid
    			}
    
  	} else {KK=K0*(mm2-strmid)+(T0/(2*pi))*(sin(2*pi*mm2)-sin(2*pi*strmid))-(T1/(2*pi))*(cos(2*pi*mm2)-cos(2*pi*strmid))}	


  if (MMLIST[mm]==min(MM1985))  {M01} else if (MMLIST[mm]==min(MM1986)) {M02} else if (MMLIST[mm]>min(MM1985)&MMLIST[mm]<min(MM1986)) {M01 + (LINF-M01)*(1-exp(-KK))} else {M02 + (LINF-M02)*(1-exp(-KK))}
}


#INITIALISE THE TAU's which are the inclusion probablities to all be a third


TAU <- list()
for (mm in 1:(MON))
	{
  		TAU[[mm]] <- matrix(1/3,NN,3)
	}


#INITIALISE THE AUXILLARY STORAGE ARRAYS FOR INSIDE LOOP ASSIGNMENT

MU2A_AUX <- c()
MUA_AUX <- c()
MUJ_AUX <- c()
for (mm in 1:MON)
	{
 	MU2A_AUX[mm] <- MU2A(mm,K0,T0,T1,LINF,M01,M02)
  	MUA_AUX[mm] <- MUA(mm,K0,T0,T1,LINF,M01,M02)
 	MUJ_AUX[mm] <- MUJ(mm,K0,T0,T1,LINF,M01,M02)
	}


#INITIALISE THE VARIANCE FUNCTIONS. THEY ARE ALL QUADRATIC FUNCTIONS OF THE MEAN FOR EACH COHORT

VARJ <- function(mm,VV,K0,T0,T1,LINF,M01,M02)
	{
		max(exp(VV[1])+VV[2]*MUJ(mm,K0,T0,T1,LINF,M01,M02)+VV[3]*MUJ(mm,K0,T0,T1,LINF,M01,M02)^2,1)
	}

VARA <- function(mm,VV,K0,T0,T1,LINF,M01,M02)
	{
		max(exp(VV[1])+VV[2]*MUA(mm,K0,T0,T1,LINF,M01,M02)+VV[3]*MUA(mm,K0,T0,T1,LINF,M01,M02)^2,1)
	}

VAR2A <- function(mm,VV,K0,T0,T1,LINF,M01,M02)
	{
		max(exp(VV[1])+VV[2]*MU2A(mm,K0,T0,T1,LINF,M01,M02)+VV[3]*MU2A(mm,K0,T0,T1,LINF,M01,M02)^2,1)
	}



#INITIALISE THE LIKELIHOOD FUNCTION

LOGLIKE <- 0
for (mm in 1:MON)
  {
	LL2<-LL[which(MM==MMLIST[mm])]
	LikGrp1<-PI1[mm]*dnorm(LL2,MU2A_AUX[mm],sqrt(VAR2A(mm,VV,K0,T0,T1,LINF,M01,M02)))
	LikGrp2<-PI2[mm]*dnorm(LL2,MUA_AUX[mm],sqrt(VARA(mm,VV,K0,T0,T1,LINF,M01,M02)))
	LikGrp3<-(1-(PI1[mm]+PI2[mm]))*dnorm(LL2,MUJ_AUX[mm],sqrt(VARJ(mm,VV,K0,T0,T1,LINF,M01,M02)))
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
  LOGOLD <- LOGLIKE                     #Assign the current likelihood value to the an old value so we can evaluate the update

 
	#UPDATE THE TAU SCORES AND THEN CALCULATE THE PIs GIVEN THESE TAU SCORES
	
	for (mm in 1:MON)
  	{
    
    
    	Top1<-PI1[mm]*dnorm(LL,MU2A_AUX[mm],sqrt(VAR2A(mm,VV,K0,T0,T1,LINF,M01,M02)))							#Calculate the inclusion probablities for group 1
    	Top2<-PI2[mm]*dnorm(LL,MUA_AUX[mm],sqrt(VARA(mm,VV,K0,T0,T1,LINF,M01,M02)))								#Calculate the inclusion probablities for group 2
    	Bot<- Top1 + Top2 +(1-(PI1[mm]+PI2[mm]))*dnorm(LL,MUJ_AUX[mm],sqrt(VARJ(mm,VV,K0,T0,T1,LINF,M01,M02)))  #Calculate the sum of each of the classes
    
    	TAU[[mm]][,1] <- Top1/Bot																				#Calculate the first group's tau scores
    	TAU[[mm]][,2] <- Top2/Bot																				#Calculate the second group's tau scores
   		TAU[[mm]][,3] <- 1-(TAU[[mm]][,1]+TAU[[mm]][,2])													  	#Calculate the third group's tau scores
    
    	PI1[mm] <- sum((MM==MMLIST[mm])*TAU[[mm]][,1])/sum(MM==MMLIST[mm])										#Work out the PIs for group 1
    	PI2[mm] <- sum((MM==MMLIST[mm])*TAU[[mm]][,2])/sum(MM==MMLIST[mm])										#Work out the PIs for group 2
  	}


  #DEFINE THE LIKELIHOOD FUNCTION OUTSIDE THE LOOP TO UPDATE THE MEANS

  PARA <- c(K0,T0,T1,LINF,M01,M02)
  OPTIFUN <- function(PARA)
  {
  	PARA[2]=0.490
	PARA[3]=0.40
  	LOGLIKE <- 0
  		for (mm in 1:MON)
  			{
			LL2<-LL[which(MM==MMLIST[mm])]
			LikGrp1<-PI1[mm]*dnorm(LL2,MU2A(mm,PARA[1],PARA[2],PARA[3],PARA[4],PARA[5],PARA[6]),sqrt(VAR2A(mm,VV,PARA[1],PARA[2],PARA[3],PARA[4],PARA[5],PARA[6])))
			LikGrp2<-PI2[mm]*dnorm(LL2,MUA(mm,PARA[1],PARA[2],PARA[3],PARA[4],PARA[5],PARA[6]),sqrt(VARA(mm,VV,PARA[1],PARA[2],PARA[3],PARA[4],PARA[5],PARA[6]))) 
			LikGrp3<-(1-(PI1[mm]+PI2[mm]))*dnorm(LL2,MUJ(mm,PARA[1],PARA[2],PARA[3],PARA[4],PARA[5],PARA[6]),sqrt(VARJ(mm,VV,PARA[1],PARA[2],PARA[3],PARA[4],PARA[5],PARA[6])))
			LOGLIKE <- LOGLIKE + sum(log(LikGrp1 + LikGrp2 + LikGrp3))
  			}
	-LOGLIKE
  }
	  
  OPTIM <- optim(PARA,OPTIFUN, control = list(maxit = 10000))														#RUN OPTIM to find the best parameters given the taus and Pis
  print("Did optim 1 converge?")
  print(OPTIM$convergence)
  
  #Assign the estimates from optim to the parameters
  
  K0 <- OPTIM$par[1]
  #T0 <- OPTIM$par[2]
  #T1 <- OPTIM$par[3]
  # Fixed thetas
  T0=0.490
 T1=0.40
  LINF <- OPTIM$par[4]
  M01 <- OPTIM$par[5]
  M02<-OPTIM$par[6]
  

  #Calculate the means for each group based on these new estimates
  
  MU2A_AUX <- c()
  MUA_AUX <- c()
  MUJ_AUX <- c()
  for (mm in 1:MON)
  {
  	#Compute the means for each month 
    MU2A_AUX[mm] <- MU2A(mm,K0,T0,T1,LINF,M01,M02)
    MUA_AUX[mm] <- MUA(mm,K0,T0,T1,LINF,M01,M02)
    MUJ_AUX[mm] <- MUJ(mm,K0,T0,T1,LINF,M01,M02)
  }
  
  
  #OPTIMISE THE VARIANCES GIVEN THE NEW MEANS
  
  #DEFINE THE LIKELIHOOD FUNCTION OUTSIDE THE LOOP TO UPDATE THE VARIANCES 
  PARA2 <- VV
  VAROPTIFUN <- function(PARA2)
  {
  LOGLIKE <- 0
  		for (mm in 1:MON)
  			{
  			#mm=1
			LL2<-LL[which(MM==MMLIST[mm])]
			LikGrp1<- PI1[mm]*dnorm(LL2,MU2A_AUX[mm],sqrt(VAR2A(mm,PARA2,K0,T0,T1,LINF,M01,M02)))
			LikGrp2<- PI2[mm]*dnorm(LL2,MUA_AUX[mm],sqrt(VARA(mm,PARA2,K0,T0,T1,LINF,M01,M02)))
			LikGrp3<-(1-(PI1[mm]+PI2[mm]))*dnorm(LL2,MUJ_AUX[mm],sqrt(VARJ(mm,PARA2,K0,T0,T1,LINF,M01,M02)))
			LOGLIKE <- LOGLIKE + sum(log(LikGrp1 + LikGrp2 + LikGrp3))
  			}
	-LOGLIKE
  }

  OPTIM2 <- optim(PARA2,VAROPTIFUN,control = list(maxit = 10000))													#RUN OPTIM to find the best variances given the taus and Pis
  VV <- OPTIM2$par
  print("Did optim 2 converge?")
  print(OPTIM2$convergence)  
    
  #COMPUTE THE LIKELIHOOD GIVEN ALL THESE NICE NEW UPDATES
      
  LOGLIKE <- 0
  for (mm in 1:MON)
  {
	LL2<-LL[which(MM==MMLIST[mm])]
	LikGrp1<-PI1[mm]*dnorm(LL2,MU2A_AUX[mm],sqrt(VAR2A(mm,VV,K0,T0,T1,LINF,M01,M02)))
	LikGrp2<-PI2[mm]*dnorm(LL2,MUA_AUX[mm],sqrt(VARA(mm,VV,K0,T0,T1,LINF,M01,M02)))
	LikGrp3<-(1-(PI1[mm]+PI2[mm]))*dnorm(LL2,MUJ_AUX[mm],sqrt(VARJ(mm,VV,K0,T0,T1,LINF,M01,M02)))
	LOGLIKE <- LOGLIKE + sum(log(LikGrp1 + LikGrp2 + LikGrp3))
  }

  
  #PRINT OUT OUR NEW UPDATES FOR THE PARAMETERS AND THE LIKELIHOOD
  
  print(c(LOGLIKE,LOGLIKE-LOGOLD))
  print(c(PARA,PARA2))
 
  
  #DRAW SOME NICE REAL TIME PLOTS SO WE CAN MONITOR CONVERGENCE
  
  par(mfrow = c(3,5))
  XX <- 0:200
  xlabNam<-c("Feb., 1985","Mar., 1985","Apr., 1985","May, 1985","Jun., 1985","Jul, 1985","Aug., 1985","Feb., 1986","Mar., 1986","Apr., 1986","May, 1986")
  for (mm in 1:MON)
  {
    hist(LL[which(MM==MMLIST[mm])],breaks=30,prob=T,xlim=c(20,210),ylim=c(0,0.055),xlab=xlabNam[mm],ylab="",main='')
    lines(XX,PI1[mm]*dnorm(XX,MU2A(mm,K0,T0,T1,LINF,M01,M02),sqrt(VAR2A(mm,VV,K0,T0,T1,LINF,M01,M02))),col='red',lwd=2.5)
    lines(XX,PI2[mm]*dnorm(XX,MUA(mm,K0,T0,T1,LINF,M01,M02),sqrt(VARA(mm,VV,K0,T0,T1,LINF,M01,M02))),col='green',lwd=2.5)
    lines(XX,(1-(PI1[mm]+PI2[mm]))*dnorm(XX,MUJ_AUX[mm],sqrt(VARJ(mm,VV,K0,T0,T1,LINF,M01,M02))),col='blue',lwd=2.5)
  }
  XX2 <- seq(0,1,length.out=100)
  YY <- K0 + T0*cos(2*pi*XX2) + T1*sin(2*pi*XX2)
  plot(XX2,((YY+abs(YY))/2),type='l',xlab="Fraction of a year since Jan 1",ylab="k (per year)",lwd=2.5)
  BLAH <- 50:160
  plot(BLAH,exp(VV[1])+VV[2]*BLAH+VV[3]*BLAH^2,type='l',xlab="Mean length (mm)", ylab="Variance",lwd=2.5,ylim=c(0,300))
  plot(density(rnorm(1000000,PARA[4],sqrt(PARA2[2]*PARA[4]+ PARA2[3]*PARA[4]^2))),xlab="Length (mm)",main="",lwd=2.5)
}






