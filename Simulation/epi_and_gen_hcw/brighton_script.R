## Attempt to infer the parameters of the model that includes health care workers
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("simulation_functions.R")
Rcpp::sourceCpp("mcmc.cpp")



data <- readRDS("mcmc_data.rds")

debug_flags <- c(0,0,0,0,0)
MCMC_options <- list("iterations" = 200,
                     "prior_parameters" = c(1.0,1.0,1.0,1.0,1e-5,1e-5),
                     "initial_chain_state" = c(0.7,0.2,0.01,0.001),
                     "proposal_variance" = c(0.0008, 0.0001),
                     "output_file" = "output.dat",
                     "debug_flags" = debug_flags,
                     "source_file" = "source.dat",
                     "loglik_file" = "loglik.dat",
                     "coltime_file" = "coltimes.dat",
                     "num_updates" = 10,
                     "matrix_file" = "matrix.dat")

hcw_ind <- rep(0,length(data$data$t_a))
#MCMC_Epi_HCW(MCMC_options, data$t_a, data$t_c, data$t_d, data$source, data$hcw_ind, data$screening_matrix)
MCMC_EPI_SOURCE(MCMC_options, data$data$t_a, data$data$t_c, data$data$t_d, data$data$source, hcw_ind, data$screening_matrix)



MCMC_options$proposal_variance[1] <- 0.001
MCMC_EPI_NS(MCMC_options, data$data$t_a, data$data$t_c, data$data$t_d, data$data$source, hcw_ind, data$screening_matrix)



hcw_data <- readRDS("mcmc_hcw_data.rds")
MCMC_EPI_NS(MCMC_options, hcw_data$t_a, hcw_data$t_c, hcw_data$t_d, hcw_data$source, hcw_data$hcw_ind, hcw_data$screening_matrix)

MCMC_options$proposal_variance[2] <- 0.00003
MCMC_EPI_SOURCE(MCMC_options, hcw_data$t_a, hcw_data$t_c, hcw_data$t_d, hcw_data$source, hcw_data$hcw_ind, hcw_data$screening_matrix)





