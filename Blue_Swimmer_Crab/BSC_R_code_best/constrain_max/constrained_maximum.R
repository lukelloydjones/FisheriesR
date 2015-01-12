t=seq(0,1,1/100)
g=0.841+0.616*cos(2*pi*t)+0.455*sin(2*pi*t)

plot(t,g)

# Calculate the maximum

theta1=0.616
theta2=0.455
M=(1/(2*pi))*acos((theta1/theta2)/(sqrt(1+(theta1/theta2)^2)))

#theta1=theta2*cos(2*pi*M)/(sqrt(1-cos(2*pi*M)^2))



val1=(1/(2*pi))*acos((theta1/(theta2))/(sqrt(1+((theta1/(theta2)))^2)))
val2=1-(1/(2*pi))*acos(-(theta1/(theta2))/(sqrt(1+((theta1/(theta2)))^2)))
abline(v=val1)
abline(v=val2)