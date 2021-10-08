## Attempt to infer the parameters of the model that includes health care workers
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("simulation_functions.r")
Rcpp::sourceCpp("mcmc.cpp")


t_a <- c(0,0,0,0,0,0)
t_c <- c(2,0,4,4,-1,2)
t_d <- c(8,8,8,8,8,8)
source <- c(2,-1,2,6,-2)
screening_matrix <- matrix(-1,nrow=6,ncol=9)
screening_matrix[1,3] <- 0
screening_matrix[1,7] <- 1
screening_matrix[2,c(1,5)] <- 1
screening_matrix[3,1] <- 0
screening_matrix[3,7] <- 1
screening_matrix[4,1] <- 0
screening_matrix[4,5] <- 1
screening_matrix[5,c(1,5)] <- 0
screening_matrix[6,1] <- 0
screening_matrix[6,4] <- 1


beta <- 0.01
alpha <- 0.005

-35*beta-8*alpha+2*log(1-exp(-(2*beta+alpha)))+log(1-exp(-(beta)))-log(2)

parameters <- c(0.9,0.1,0.01,0.005)
hcw_ind <- c(0,0,0,0,0,1)
CalculateLogLikelihood_R(parameters, t_a, t_c, t_d, source, hcw_ind, screening_matrix)


debug_flags <- c(0,0,0,0,0)
MCMC_options <- list("iterations" = 2,
                     "prior_parameters" = c(1.0,1.0,1.0,1.0,1e-5,1e-5),
                     "initial_chain_state" = c(0.9,0.1,0.01,0.005),
                     "proposal_variance" = c(0.01, 0.0005),
                     "output_file" = "output.dat",
                     "debug_flags" = debug_flags,
                     "source_file" = "source.dat",
                     "loglik_file" = "loglik.dat",
                     "coltime_file" = "coltimes.dat",
                     "num_updates" = 10,
                     "matrix_file" = "matrix.dat")

#MCMC_Epi_HCW(MCMC_options, data$t_a, data$t_c, data$t_d, data$source, data$hcw_ind, data$screening_matrix)
MCMC_EPI_SOURCE(MCMC_options, t_a, t_c, t_d, source, hcw_ind, screening_matrix)




