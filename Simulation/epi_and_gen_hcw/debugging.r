## This file is for debugging problematic data sets in the simulation study
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("simulation_functions.R")
Rcpp::sourceCpp("mcmc.cpp")



current_problem_data_set <- 487
current_res <- readRDS(paste("problem_data_sets/results_",current_problem_data_set,".rds",sep=""))
epi_data <- current_res$simulation_info$epi_data

debug_flags <- c(0,0,0,0,0,0,0)
MCMC_options <- list("iterations" = 1000,
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



MCMC_EPI_SOURCE(MCMC_options, epi_data$t_a, epi_data$t_c, epi_data$t_d, epi_data$source, epi_data$hcw_ind, epi_data$screening_matrix)

res_epi <- read.table("output.dat")
source_epi <- read.table("source.dat")

plot(res_epi$V1,type="l")

#variant_numbers <- ReturnVariantNumbersCutoff(current_res$simulation_info$gen_data$genetic_matrix, 5)

MCMC_EPI_GEN_HCW(MCMC_options, 2800000, current_res$simulation_info$epi_data$t_a, current_res$simulation_info$epi_data$t_c, current_res$simulation_info$epi_data$t_d, 
                 current_res$simulation_info$epi_data$source, current_res$simulation_info$epi_data$hcw_ind, current_res$simulation_info$epi_data$screening_matrix, 
                 current_res$simulation_info$gen_data$genetic_ids, current_res$simulation_info$gen_data$sample_times, current_res$simulation_info$variant_numbers, current_res$simulation_info$gen_data$genetic_matrix)


res <- read.table("output.dat")
source <- read.table("source.dat")

gen_source_data <- CalculateSourceClusterProportion(epi_data, source)
epi_source_data <- CalculateSourceClusterProportion(epi_data, source_epi)
#res <- read.table("output.dat")


res <- current_res
ever_infected <- which(res$simulation_info$epi_data$t_c != -1 & res$simulation_info$epi_data$hcw_ind == 0)
source_table <- data.frame("id" = numeric(length(ever_infected)),
                          "inferred_source" = numeric(length(ever_infected)),
                          "true_source" = numeric(length(ever_infected)),
                          "posterior_probability" = numeric(length(ever_infected)))


source <- read.table("source.dat")

for(i in 1:length(ever_infected)) {
  current_target <- ever_infected[i]
  posterior_source_distribution <- table(source[,current_target])
  most_likely_source <- as.numeric(names(posterior_source_distribution)[which.max(posterior_source_distribution)])
  if(most_likely_source >=0 ) {
    most_likely_source <- most_likely_source + 1
  }
  source_table[i,1] <- current_target
  source_table[i,2] <- most_likely_source
  source_table[i,3] <- current_res$simulation_info$epi_data$true_source[current_target]
  source_table[i,4] <- as.numeric(max(table(source[,current_target]))/sum(table(source[,current_target])))
}
sum(source_table[,2]==source_table[,3])/nrow(source_table)





t_a <- epi_data$t_a
t_c <- epi_data$true_coltimes
t_d <- epi_data$t_d
source <- epi_data$true_source


ever_infected <- which(t_c != -1)

plot(0,type="n", ylim=c(0,length(ever_infected)), xlim=c(0,max(t_d)))
for(i in 1:length(ever_infected)) {
  current_person <- ever_infected[i]
  current_source <- source[current_person]
  segments(t_a[current_person],i,t_c[current_person],i,lwd=2)
  segments(t_c[current_person],i,t_d[current_person],i,col=2, lwd=2)
  if(current_source != -1) {
    source_loc <- which(current_source == ever_infected)
    arrows(t_c[current_person], source_loc, t_c[current_person], i)
  }
  text(t_a[current_person]-3,i,current_person)
}

genetic_ids <- gen_data$genetic_ids
sample_times <- gen_data$sample_times
for(i in 1:length(genetic_ids)) {
  current_gen_id <- genetic_ids[i]
  current_sample_time <- sample_times[i]
  text(current_sample_time, which(current_gen_id==ever_infected)+0.4,i,cex=0.7)
}





string <- "1 1 2 2 4 4 6 6 6 7 7 8 10 14 15 16 18 20 20 21 23 24 26 26 26 28 28 29 30 30 30 32 32 32 33 34 34 38 38 40 42 44 46 46 47 48 48 48 50 51 51 51 51 51 54 54 57 59 60 62 63 63 64 64 65 68 68 68 68 68 69 69 70 71 72 73 75 75 76 76 77 78 79 81 84 84 84 87 87 89 90 90 90 94 95 95 96 98 98 100 0 0 0 0
-1 -1 -1 -1 -1 10 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 21 -1 -1 -1 -1 -1 -1 -1 -1 -1 33 -1 -1 -1 -1 36 34 -1 38 38 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 54 54 -1 -1 -1 -1 -1 -1 -1 -1 -1 68 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 76 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 70 14 21 84
6 6 9 10 12 11 12 11 10 14 17 11 16 17 19 25 21 24 29 34 28 33 33 34 36 36 36 34 31 37 37 37 37 44 42 36 39 46 47 44 51 51 54 53 51 57 58 53 60 62 57 58 58 59 59 58 62 68 71 69 71 66 71 67 76 75 75 71 73 75 78 76 82 81 76 77 78 85 83 82 82 84 84 86 93 91 96 92 91 96 95 93 102 101 102 108 104 105 106 107 108 108 108 108
-2 -2 -2 -2 -2 -1 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -1 -2 -2 -2 -2 -2 -2 -2 -2 -2 -1 -2 -2 -2 -2 -1 -1 -2 -1 -1 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -1 -1 -2 -2 -2 -2 -2 -2 -2 -2 -2 -1 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -1 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -1 -1 -1 -1
5 101 19 101 102 19 101 102 29 35 101 102 34 37 38 34 37 38 101 102 37 38 38 101 102 54 55 101 54 55 101 102 65 100 101 102 65 65 78 100 101 102 78 78 100 101 103 101 100 102 101 102 103
10 14 21 21 21 24 28 28 33 34 35 35 36 38 38 39 41 41 42 42 44 44 47 49 49 54 54 56 57 57 63 63 68 70 70 70 71 74 76 77 77 77 79 82 84 84 84 91 98 98 105 105 105"


PlotNetworkFromString(string)



epi_source <- CalculateSourceClusterProportion(epi_data, current_res$epi_res$source)
gen_source <- CalculateSourceClusterProportion(epi_data, current_res$gen_res$source)




## Construct ROC curve



ReturnTPR_FPR <- function(epi_data, source, x) {
  from_edge <- c()
  to_edge <- c()
  for(i in 1:ncol(source)) {
    posterior_source_distribution <- table(source[,i]) / sum(table(source[,i]))
    links <- which(posterior_source_distribution > x & names(posterior_source_distribution) != -2)
    links <- as.numeric(names(links))
    links[links >= 0] <- links[links >= 0] + 1
    #if(length(links)>0) browser()
    for(j in links) {
      if(j == -1) {
        from_edge <- c(from_edge, i)
        to_edge <- c(to_edge, i)
      } else {
        from_edge <- c(from_edge, j)
        to_edge <- c(to_edge, i)
      }
    }
  }
  
  unique_ids <- unique(c(from_edge,to_edge))
  TP <- 0
  FP <- 0
  FN <- 0
  TN <- 0
  #browser()
  for(i in 1:length(unique_ids)) {
    for(j in 1:length(unique_ids)) {
      person_i <- unique_ids[i]
      person_j <- unique_ids[j]
      source_j <- epi_data$true_source[person_j]
      if(source_j == -1) {
        # The current person j is an importation
        if(person_i == person_j) {
          true_edge <- TRUE
        } else {
          true_edge <- FALSE
        }
        
        # Check if the person is inferred to be an importation
        if(sum(from_edge==person_i & to_edge==person_j) > 0) {
          inferred_edge <- TRUE
        } else {
          inferred_edge <- FALSE
        }
        
      } else {
        # The current person is an acquisition
        if(person_i == source_j) {
          true_edge <- TRUE
        } else {
          true_edge <- FALSE
        }
        
        # Check if the person is inferred to be an acquisiton from the correct source
        if(sum(from_edge==person_i & to_edge==person_j) > 0) {
          inferred_edge <- TRUE
        } else {
          inferred_edge <- FALSE
        }
      }
      
      if(true_edge && inferred_edge) {
        TP <- TP + 1
      }
      
      if(true_edge && !inferred_edge) {
        FN <- FN + 1
      }
      
      if(!true_edge && inferred_edge) {
        #browser()
        FP <- FP + 1
      }
      
      if(!true_edge && !inferred_edge) {
        TN <- TN + 1
      }
      
    }
  }
  
  #browser()
  TPR <- TP/(TP+FN)
  FPR <- FP/(FP+TN)
  
  out <- list("TPR" = TPR,
              "FPR" = FPR)
  return(out)
}

x <- seq(from=0.01,to=0.99,len=101)
ReturnTPR_FPR(epi_data, source, 0.5)
tpr_fpr <- sapply(x, ReturnTPR_FPR, epi_data=epi_data, source=source)


p_L_values <- seq(from=0,to=1,len=101)
gen_resolution <- sapply(p_L_values, CalculateTrueTransmissionProportion, source=source, epi_data=epi_data)
epi_resolution <- sapply(p_L_values, CalculateTrueTransmissionProportion, source=source_epi, epi_data=epi_data)
plot(p_L_values,gen_resolution,type="s")
lines(p_L_values,epi_resolution, type="s", col="blue")





