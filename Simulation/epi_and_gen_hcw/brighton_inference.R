## Attempt to infer the parameters of the model that includes health care workers
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("simulation_functions.r")
Rcpp::sourceCpp("mcmc.cpp")


data <- readRDS("hcw_mcmc_data.rds")


debug_flags <- c(0,0,0,0,1)
MCMC_options <- list("iterations" = 300,
                     "prior_parameters" = c(1.0,1.0,1.0,1.0,1e-5,1e-5),
                     "initial_chain_state" = c(0.7,0.2,0.01,0.001),
                     "proposal_variance" = c(0.001, 0.0001),
                     "output_file" = "output.dat",
                     "debug_flags" = debug_flags,
                     "source_file" = "source.dat",
                     "loglik_file" = "loglik.dat",
                     "coltime_file" = "coltimes.dat",
                     "num_updates" = 10,
                     "matrix_file" = "matrix.dat")

#MCMC_Epi_HCW(MCMC_options, data$t_a, data$t_c, data$t_d, data$source, data$hcw_ind, data$screening_matrix)
hcw_ind <- rep(0,length(data$t_a)); hcw_ind[1934:2068] <- 1
MCMC_Epi_HCW(MCMC_options, data$t_a, data$t_c, data$t_d, data$source_vector, hcw_ind, data$screening_matrix)



pos_swabs <- which(data$screening_matrix==1,arr.ind=T)
unique(pos_swabs[,1])
