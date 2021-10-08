rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("simulation_functions.r")
#Rcpp::sourceCpp("mcmc2.cpp")

N <- 300
D <- 300
length_of_stay <- 7
p <- 0.05
z <- 0.95
a0 <- 1e-10
a1 <- 0.01
a2 <- 0.0005
num_col_hcw <- 10


outbreak <- simulateWard(N,D,length_of_stay,p, z, a0, a1, a2, num_col_hcw)


importations <- sum(outbreak$source == -1,na.rm=T)
patient_acq <- sum(outbreak$source>=0 & outbreak$source <= 300, na.rm=T)
hcw_acq <- sum(outbreak$source > 300, na.rm=T)


max_iter <- 2000
importations_vec <- numeric(max_iter)
patient_acq_vec <- numeric(max_iter)
hcw_acq_vec <- numeric(max_iter)
for(i in 1:max_iter) {
  outbreak <- simulateWard(N,D,length_of_stay,p, z, a0, a1, a2, num_col_hcw)
  importations_vec[i] <- sum(outbreak$source == -1,na.rm=T)
  patient_acq_vec[i] <- sum(outbreak$source>=0 & outbreak$source <= 300, na.rm=T)
  hcw_acq_vec[i] <- sum(outbreak$source > 300, na.rm=T)
}


hist(importations_vec)
hist(patient_acq_vec)
hist(hcw_acq_vec)
hist(importations_vec+patient_acq_vec+hcw_acq_vec)



## Simulate outbreak and screen the data


patient_screening_interval <- 3
hcw_screening_interval <- 7
outbreak <- simulateWard(N,D,length_of_stay,p, z, a0, a1, a2, num_hcw, col_hcw_prop)
screening_matrix <- screenData(outbreak, patient_screening_interval, hcw_screening_interval, z)






