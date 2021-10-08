## Attempt to infer the parameters of the model that includes health care workers
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("simulation_functions.R")
Rcpp::sourceCpp("mcmc.cpp")

population_size <- 100
L <- 100
LOS <- 7
p <- 0.1
z <- 0.70
beta_p <- 0.01
beta_h <- 0.0005
num_col_hcw <- 5

patient_screening_interval <- 3
hcw_screening_interval <- 7


epi_data <- SimulateEpiDataMCMC(population_size,L,LOS,p,z,beta_p,beta_h,num_col_hcw,patient_screening_interval,hcw_screening_interval)

epi_data$source <- c(epi_data$source, rep(-1,sum(epi_data$hcw_ind)))
epi_data$true_source <- c(epi_data$true_source, rep(-1,sum(epi_data$hcw_ind)))

variant_numbers <- rep(1,sum(epi_data$screening_matrix==1))
average_import_distance <- 10
average_variant_distance <- 20
N <- 2800000
mutation_rate <- 3.3*10^(-6)/365
gen_data <- SimulateGeneticData_WHD(epi_data, variant_numbers, average_import_distance, average_variant_distance, N, mutation_rate) 



debug_flags <- c(0,0,0,0,0,0,0)
MCMC_options <- list("iterations" = 5000,
                     "prior_parameters" = c(1.0,1.0,1.0,1.0,1e-5,1e-5,1e-5,1.0,1.0),
                     "initial_chain_state" = c(0.7,0.2,0.01,0.001,3.3*10^(-6)/365,10),
                     "proposal_variance" = c(0.005, 0.001,0.000000005),
                     "output_file" = "output.dat",
                     "debug_flags" = debug_flags,
                     "source_file" = "source.dat",
                     "loglik_file" = "loglik.dat",
                     "coltime_file" = "coltimes.dat",
                     "num_updates" = 10,
                     "matrix_file" = "matrix.dat")

#MCMC_EPI_GEN_HCW(MCMC_options, N, epi_data$t_a, epi_data$true_coltimes, epi_data$t_d, epi_data$true_source, epi_data$hcw_ind, epi_data$screening_matrix, 
#                 gen_data$genetic_ids, gen_data$sample_times, variant_numbers, gen_data$genetic_matrix)

MCMC_EPI_GEN_HCW(MCMC_options, N, epi_data$t_a, epi_data$t_c, epi_data$t_d, epi_data$source, epi_data$hcw_ind, epi_data$screening_matrix, 
                 gen_data$genetic_ids, gen_data$sample_times, variant_numbers, gen_data$genetic_matrix)


MCMC_EPI_SOURCE(MCMC_options, epi_data$t_a, epi_data$t_c, epi_data$t_d, epi_data$source, epi_data$hcw_ind, epi_data$screening_matrix)

#void MCMC_EPI_SOURCE(List MCMC_options, IntegerVector t_a, IntegerVector t_c, IntegerVector t_d, IntegerVector source, 
#                     IntegerVector hcw_ind, IntegerMatrix screening_matrix)

res <- read.table("output.dat")

dev.new()
par(mfrow=c(4,2))
plot(res$V1,type="l", main="z")
abline(h=z,col=2)
plot(res$V2,type="l", main="p")
abline(h=p,col=2)
plot(res$V3,type="l", main="beta_p")
abline(h=beta_p,col=2)
plot(res$V4,type="l", main="beta_h")
abline(h=beta_h,col=2)
plot(res$V5,type="l", main="lambda")
abline(h=mutation_rate,col=2)
plot(res$V6,type="l", main="mu")
abline(h=average_import_distance,col=2)


coltimes <- read.table("coltimes.dat")
plot(coltimes$V1,type="l", main="colonisation time sum")
abline(h=sum(epi_data$true_coltimes[epi_data$true_coltimes != -1]),col=2)


patient_ever_infected <- which(epi_data$t_c != -1 & epi_data$hcw_ind == 0)
source <- read.table("source.dat")
source[source >= 0] <- source[source >= 0] + 1
source_proportion <- apply(source, 1, function(x) sum(x[patient_ever_infected]==epi_data$true_source[patient_ever_infected])/length(patient_ever_infected))
plot(source_proportion, type="l", main="Proportion of correctly identified sources")


