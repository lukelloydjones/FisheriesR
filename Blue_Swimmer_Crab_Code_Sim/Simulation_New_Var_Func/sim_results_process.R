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

file_names <-  read.table("~/Dropbox/Git_Repos/Fisheries_R_Scripts/Blue_Swimmer_Crab_Code_Sim/Simulation_New_Var_Func/Sim_Move_Res_Old_2/Sim_Names.txt",header=F)


# For each of the files read them in and process into distinct matrices for the parameters of interest

para_mat<- matrix(0, nrow = dim(file_names)[1], ncol = 12)
pi1_mat <- matrix(0, nrow = dim(file_names)[1], ncol = 12)
pi2_mat <- matrix(0, nrow = dim(file_names)[1], ncol = 12)
pi3_mat <- matrix(0, nrow = dim(file_names)[1], ncol = 12)
mu1_mat <- matrix(0, nrow = dim(file_names)[1], ncol = 12)
mu2_mat <- matrix(0, nrow = dim(file_names)[1], ncol = 12)
mu3_mat <- matrix(0, nrow = dim(file_names)[1], ncol = 12)


# Set the working directory for the simulation results

setwd("~/Dropbox/Git_Repos/Fisheries_R_Scripts/Blue_Swimmer_Crab_Code_Sim/Simulation_New_Var_Func/Sim_Move_Res_Old_2")


# Run through all the file names for each simulation and put each of the lines in the respective matrix above

for (i in seq(1,dim(file_names)[1]))
	{
	
	fil_name      =  as.character(file_names[i, ])			    # Determine the file name
	results       <- read.table(fil_name, header = F, fill = T)	# Read in the results for this file
	para_mat[i, ] <- as.numeric(results[1, ])
	pi1_mat[i, ]  <- as.numeric(results[2, ])
	pi2_mat[i, ]  <- as.numeric(results[3, ])
	pi3_mat[i, ]  <- as.numeric(results[4, ])
	mu1_mat[i, ]  <- as.numeric(results[5, ])
	mu2_mat[i, ]  <- as.numeric(results[6, ])
	mu3_mat[i, ]  <- as.numeric(results[7, ])
	}


# Calculte the means and standard deviations of each of the matrices

# Model estimates

mean_ests<-colMeans(para_mat[-1,-c(10,11,12)])
std_errs<-apply(para_mat[-1,],2,sd)


# Mixing proportion estimates means and standard errors

pi1_mat_means <- colMeans(pi1_mat[-1,])
pi2_mat_means <- colMeans(pi2_mat[-1,])
pi3_mat_means <- colMeans(pi3_mat[-1,])

pi1_mat_means_se <- apply(pi1_mat[-1, ], 2, sd)
pi2_mat_means_se <- apply(pi2_mat[-1, ], 2, sd)
pi3_mat_means_se <- apply(pi3_mat[-1, ], 2, sd)


# Mean length estimates means and standard errors

mu1_mat_means <- colMeans(mu1_mat[-1, ])
mu2_mat_means <- colMeans(mu2_mat[-1, ])
mu3_mat_means <- colMeans(mu3_mat[-1, ])

mu1_mat_means_se <- apply(mu1_mat[-1, ], 2, sd)
mu2_mat_means_se <- apply(mu2_mat[-1, ], 2, sd)
mu3_mat_means_se <- apply(mu3_mat[-1, ], 2, sd)


# Print them all out with & in between to make it easier to transfer into latex


# Model estimates

paste(as.character(signif(mean_ests,4)),collapse=" & ")
paste(as.character(signif(std_errs,4)), collapse=" & ")


# Mixing proportion estimates

paste( as.character(signif(pi1_mat_means,3)), collapse = " & ")
paste( as.character(signif(pi2_mat_means,3)), collapse = " & ")
paste( as.character(signif(pi3_mat_means,3)), collapse = " & ")


# Mixing proportion SE estimates

paste( as.character( signif(pi1_mat_means_se, 2)), collapse = " & ")
paste( as.character( signif(pi2_mat_means_se, 2)), collapse = " & ")
paste( as.character( signif(pi3_mat_means_se, 2)), collapse = " & ")

 
# Mean length estimates

paste(as.character(signif(mu1_mat_means, 4)), collapse = " & ")
paste(as.character(signif(mu2_mat_means, 4)), collapse = " & ")
paste(as.character(signif(mu3_mat_means, 4)), collapse = " & ")


# Mixing proportion SE estimates

paste(as.character(signif(mu1_mat_means_se, 3)), collapse = " & ")
paste(as.character(signif(mu2_mat_means_se, 3)), collapse = " & ")
paste(as.character(signif(mu3_mat_means_se, 3)), collapse = " & ")


