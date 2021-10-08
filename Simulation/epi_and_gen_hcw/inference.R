## Attempt to infer the parameters of the model that includes health care workers
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("simulation_functions.r")
Rcpp::sourceCpp("mcmc.cpp")



N <- 300
L <- 300
LOS <- 7
p <- 0.3
z <- 0.70
beta_p <- 0.01
beta_h <- 0.0005
num_col_hcw <- 0

patient_screening_interval <- 3
hcw_screening_interval <- 7


data <- SimulateEpiDataMCMC(N,L,LOS,p,z,beta_p,beta_h,num_col_hcw,patient_screening_interval,hcw_screening_interval) 





debug_flags <- c(0,0,0,0,1)
MCMC_options <- list("iterations" = 5000,
                     "prior_parameters" = c(1.0,1.0,1.0,1.0,1e-5,1e-5),
                     "initial_chain_state" = c(0.7,0.2,0.01,0.001),
                     "proposal_variance" = c(0.01, 0.0005),
                     "output_file" = "output.dat",
                     "debug_flags" = debug_flags,
                     "source_file" = "source.dat",
                     "loglik_file" = "loglik.dat",
                     "coltime_file" = "coltimes.dat",
                     "num_updates" = 10,
                     "matrix_file" = "matrix.dat")

#MCMC_Epi_HCW(MCMC_options, data$t_a, data$t_c, data$t_d, data$source, data$hcw_ind, data$screening_matrix)
MCMC_EPI_SOURCE(MCMC_options, data$t_a, data$true_coltimes, data$t_d, data$true_source, data$hcw_ind, data$screening_matrix)


res <- read.table("output.dat")
coltimes <- read.table("coltimes.dat")
source <- read.table("source.dat")

MCMC_EPI_NS(MCMC_options, data$t_a, data$true_coltimes, data$t_d, data$true_source, data$hcw_ind, data$screening_matrix)

res_ns <- read.table("output.dat")
coltimes_ns <- read.table("coltimes.dat")



### Source plots
plot(coltimes$V1,type="l")
#coltime_sums <- apply(coltimes,1,function(x) sum(x[x!=-1]))
#plot(coltime_sums,type="l")
abline(h=sum(data$true_coltimes[data$true_coltimes != -1]), col=2)
plot(res$V1,type="l")
abline(h=z,col=2)
plot(res$V2,type="l")
abline(h=p,col=2)
plot(res$V3,type="l")
abline(h=beta_p,col=2)
plot(res$V4,type="l")
abline(h=beta_h,col=2)

## No source plots
#plot(coltimes_ns$V1,type="l")
coltime_sums_ns <- apply(coltimes_ns,1,function(x) sum(x[x!=-1]))
plot(coltime_sums_ns,type="l")
abline(h=sum(data$true_coltimes[data$true_coltimes != -1]), col=2)
plot(res_ns$V1,type="l")
abline(h=z,col=2)
plot(res_ns$V2,type="l")
abline(h=p,col=2)
plot(res_ns$V3,type="l")
abline(h=beta_p,col=2)
plot(res_ns$V4,type="l")
abline(h=beta_h,col=2)



## source analysis
ever_infected <- which(data$true_coltimes != -1 & data$hcw_ind == 0)
posterior_source_matrix_summary <- matrix(nrow=length(ever_infected), ncol=4)
for(i in 1:length(ever_infected)) {
  current_target <- ever_infected[i]
  posterior_source_distribution <- table(source[,current_target])
  most_likely_source <- as.numeric(names(posterior_source_distribution)[which.max(posterior_source_distribution)])
  if(most_likely_source >=0 ) {
    most_likely_source <- most_likely_source + 1
  }
  true_source <- data$true_source[current_target]
  posterior_probability <- posterior_source_distribution[which.max(posterior_source_distribution)]/sum(posterior_source_distribution)
  posterior_source_matrix_summary[i,] <- c(current_target, true_source,most_likely_source, posterior_probability)
}
posterior_source_matrix_summary

sum(posterior_source_matrix_summary[,2]==posterior_source_matrix_summary[,3])/nrow(posterior_source_matrix_summary)


table(source[,ever_infected[7]])



z_source <- data.frame('z'=res$V1)
z_nosource <- data.frame('z'=res_ns$V1)
z_source$type <- 'source'
z_nosource$type <- 'no_source'

p_source <- data.frame('p'=res$V2)
p_nosource <- data.frame('p'=res_ns$V2)
p_source$type <- 'source'
p_nosource$type <- 'no_source'

beta_p_source <- data.frame('beta'=res$V3)
beta_p_nosource <- data.frame('beta'=res_ns$V3)
beta_p_source$type <- 'source'
beta_p_nosource$type <- 'no_source'

beta_h_source <- data.frame('beta'=res$V4)
beta_h_nosource <- data.frame('beta'=res_ns$V4)
beta_h_source$type <- 'source'
beta_h_nosource$type <- 'no_source'

coltime_source <- data.frame('coltime'=coltimes$V1)
coltime_nosource <- data.frame('coltime'=coltimes_ns$V1)
coltime_source$type <- 'source'
coltime_nosource$type <- 'no_source'

complete_z <- rbind(z_source,z_nosource)
complete_p <- rbind(p_source,p_nosource)
complete_beta_p <- rbind(beta_p_source,beta_p_nosource)
complete_beta_h <- rbind(beta_h_source,beta_h_nosource)
complete_coltime <- rbind(coltime_source,coltime_nosource)
true_coltimesum <- sum(data$true_coltimes[data$true_coltimes != -1])
library(ggplot2)


ggplot(complete_coltime, aes(coltime, fill=type)) + geom_density(alpha=0.2) + geom_vline(xintercept = true_coltimesum, color='red')
ggplot(complete_z, aes(z, fill=type)) + geom_density(alpha=0.2) + geom_vline(xintercept = 0.7, color='red')
ggplot(complete_p, aes(p, fill=type)) + geom_density(alpha=0.2) + geom_vline(xintercept = 0.05, color='red')
ggplot(complete_beta_p, aes(beta, fill=type)) + geom_density(alpha=0.2) + geom_vline(xintercept = 0.005, color='red')
ggplot(complete_beta_h, aes(beta, fill=type)) + geom_density(alpha=0.2) + geom_vline(xintercept = 0.0005, color='red')





ever_infected <- which(data$true_coltimes != -1)
par(mfrow=c(2,1))
i <- 18
plot(coltimes[,ever_infected[i]], type="l")
plot(coltimes_ns[,ever_infected[i]], type="l")
par(mfrow=c(1,1))


plot(apply(coltimes_ns, 1, function(x) sum(x!=-1)),type="l")


#ever_inf <- which(data$true_coltimes != -1)
#coltime_sum <- apply(coltimes[ever_inf],1,sum)

plot(coltimes$V1,type="l")
abline(h=sum(data$true_coltimes[data$true_coltimes != -1]), col=2)
plot(res$V1,type="l")
abline(h=z,col=2)
plot(res$V2,type="l")
abline(h=p,col=2)
plot(res$V3,type="l")
abline(h=beta_p,col=2)
plot(res$V4,type="l")
abline(h=beta_h,col=2)



