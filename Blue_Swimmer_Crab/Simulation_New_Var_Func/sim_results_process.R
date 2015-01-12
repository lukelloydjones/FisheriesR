######################################################################################
######################################################################################
######################################################################################
###																				   ###
###			     CODE TO PROCESS OUTPUT OF CRAB SIMULATION 						   ###
###																				   ###	
######################################################################################
######################################################################################
######################################################################################


# Read in the file containing all the names of the simulation results

file_names<-read.table("~/Dropbox/AAUni/APhD/Blueswimmer/Simulation_New_Var_Func/crab_sim_results/file_names_crab_sim.txt",header=F)


# For each of the file read them in and process into distinct matrices for the parameters of interest

para_mat<-matrix(0,nrow=dim(file_names)[1],ncol=12)
pi1_mat<-matrix(0,nrow=dim(file_names)[1],ncol=12)
pi2_mat<-matrix(0,nrow=dim(file_names)[1],ncol=12)
pi3_mat<-matrix(0,nrow=dim(file_names)[1],ncol=12)
mu1_mat<-matrix(0,nrow=dim(file_names)[1],ncol=12)
mu2_mat<-matrix(0,nrow=dim(file_names)[1],ncol=12)
mu3_mat<-matrix(0,nrow=dim(file_names)[1],ncol=12)


# Run through all the file names and put each of the lines in the respective matrix above

# Set the working directory

setwd("~/Dropbox/AAUni/APhD/Blueswimmer/Simulation_New_Var_Func/crab_sim_results/")

for (i in seq(1,dim(file_names)[1]))
{
	#i=1
	fil_name=as.character(file_names[i,])			# Determine the file name
	results<-read.table(fil_name,header=F,fill=T)	# Read in the results for this file
	para_mat[i,]<-as.numeric(results[1,])
	pi1_mat[i,]<-as.numeric(results[2,])
	pi2_mat[i,]<-as.numeric(results[3,])
	pi3_mat[i,]<-as.numeric(results[4,])
	mu1_mat[i,]<-as.numeric(results[5,])
	mu2_mat[i,]<-as.numeric(results[6,])
	mu3_mat[i,]<-as.numeric(results[7,])
}


# Calculte the means and standard deviations of each of the matrices

colMeans(para_mat)
apply(para_mat,2,sd)
summary(para_mat[,1])













