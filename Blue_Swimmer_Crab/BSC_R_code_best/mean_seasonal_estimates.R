#######GETTING THE ESTIMATES FOR WAYNE#######GETTING THE ESTIMATES FOR WAYNE
#######GETTING THE ESTIMATES FOR WAYNE#######GETTING THE ESTIMATES FOR WAYNE
#######GETTING THE ESTIMATES FOR WAYNE#######GETTING THE ESTIMATES FOR WAYNE

MM<-seq(0,1,1/12)
K0=0.99101351; T0=1.17902223;   T1=0.68237758;  #MALES
K0=1.1056427; T0=-0.1631402 ;   T1=0.5722303;	#FEMALES
t<-seq(0,1,0.01)
G<-K0 +T1*sin(2*pi*t)+ T0*cos(2*pi*t)
plot(t,((G+abs(G))/2),type='l',xlab="Fraction of a year since Jan 1",ylab="k (1/year)",lwd=2.5)
abline(v=c(seq(0,1,1/12)))


MM<-seq(0,1,1/12)
MM<-seq(0,1,1/12)
mns<-array(0,12)
for (i in(seq(1:12)))
{
M1<-K0*(MM[i+1]-MM[i]) +   (T0/(2*pi))*(sin(2*pi*MM[i+1])-sin(2*pi*MM[i]))-(T1/(2*pi))*(cos(2*pi*MM[i+1])-cos(2*pi*MM[i]))
Kavg<-M1/(MM[i+1]-MM[i])
mns[i]<-max(c(M1/(MM[i+1]-MM[i])),0)
lineseg<-seq(MM[i],MM[i+1],0.005)
length(lineseg)
lines(seq(MM[i],MM[i+1],0.005),array(Kavg,length(lineseg)))
}

round(mns,3)