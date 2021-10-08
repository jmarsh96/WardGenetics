# simulates a ward of N individuals over D days with healthcare workers included
# length_of_stay is the average time in the ward
# p is the probability of being colonised on admission
# z is the test sensitivity 
# a0 is the background rate of infection
# a1 is the contact rate with other infectives on the ward
# a2 is the contact rate with health care workers
# num_col_hcw is the healthcare workers ever to be colonised on the ward


simulateWard <- function(N,D,length_of_stay=7,p=0.05, z=0.95, a0 = 0.000000001, a1 = 0.005, a2=0.001, num_col_hcw) {
  
  t_a <- numeric(N+num_col_hcw)
  t_c <- rep(NA,N+num_col_hcw)
  t_d <- numeric(N+num_col_hcw)
  source <- rep(NA,N)
  #browser()
  hcw_indicator <- c(rep(0,N),rep(1,num_col_hcw))
  
  
  # generate admission times
  t_a[1:N] <- sort(sample(1:D,N,replace=T))
  t_a[1:N] <- t_a[1:N] - min(t_a)
  
  if(num_col_hcw > 0) {
    t_a[(N+1):(N+num_col_hcw)] <- 0
  }
  
  
  # sample importations
  positive_on_admission <- sample(c(TRUE,FALSE),N,replace=T,prob=c(p,1-p))
  t_c[positive_on_admission] <- t_a[positive_on_admission]
  source[positive_on_admission] <- -1
  
  # generate discharge times
  t_d[1:N] <- t_a[1:N] + rpois(N,length_of_stay)
  maxD <- max(t_d)
  
  if(num_col_hcw > 0) {
    t_d[(N+1):(N+num_col_hcw)] <- maxD
    t_c[(N+1):(N+num_col_hcw)] <- sample(0:maxD,size=num_col_hcw,replace=T)
  }

  
  
  #browser()
  #go through each day
  for(i in 1:maxD) {
    #if(i == 28) browser()
    susceptibles <- countSusceptible(t_a[1:N],t_c[1:N],t_d[1:N],i)
    colonised <- countColonisedR(t_a[1:N],t_c[1:N],t_d[1:N],i)
    if(num_col_hcw > 0) {
      colonised_hcw <- which(t_c[(N+1):(N+num_col_hcw)] < i) + N
    } else {
      colonised_hcw <- c()
    }
    
    current_rate <- a0 + length(colonised)*a1 + length(colonised_hcw)*a2
    for(j in susceptibles) {
      if(runif(1)<(1-exp(-1*current_rate))) {
        # colonisation has occured
        t_c[j] <- i
        
        #choose a source
        U <- runif(1)
        if(U <(a0/current_rate)) {
          #background infection
          source[j] <- j
        } else if(U<((a0+length(colonised)*a1)/current_rate)) {
          # colonisation from a patient
          source[j] <- sampleWithoutSurprises(colonised)
        } else {
          # colonisation from a HCW
          source[j] <- sampleWithoutSurprises(colonised_hcw)
        }
      }
    }
    
  }
  importations <- rep(0,N)
  importations[positive_on_admission] <- 1
  
  
  
  out <- list("t_a" = t_a, "t_c" = t_c, "t_d" = t_d, "importations" = importations, "source" = source, "hcw_indicator" = hcw_indicator)
  return(out)
}



## Screens the outbreak simulated from the function 

screenData <- function(data, patient_screening_interval, hcw_screening_interval, z) {
  hcw_indicator <- data$hcw_indicator
  N <- length(data$t_a)
  maxD <- max(data$t_d)
  results_matrix <- matrix(NA, nrow=N, ncol=maxD+1)
  #browser()
  for(i in 1:N) {
    if(hcw_indicator[i]==0) {
      test_days <- seq(from=data$t_a[i], to=data$t_d[i], by=patient_screening_interval)
    } else {
      test_days <- seq(from=data$t_a[i], to=data$t_d[i], by=hcw_screening_interval)
    }
    
    if(is.na(data$t_c[i])){
      results_matrix[i,test_days+1] <- 0
    } else {
      for(t in test_days) {
        if(t >= data$t_c[i]) {
          results_matrix[i,t+1] <- rbinom(1,1,z)
        } else {
          results_matrix[i,t+1] <- 0
        }
      }
    }
  }
  return(results_matrix)
}


countSusceptible <- function(t_a,t_c,t_d,t) {
  currently_inward <- intersect(which(t_a <= t),which(t_d >= t))
  susceptible <- intersect(currently_inward, which(is.na(t_c)))
  return(susceptible)
}

countColonisedR <- function(t_a,t_c,t_d,t) {
  currently_inward <- intersect(which(t_a <= t),which(t_d >= t))
  colonised <- intersect(currently_inward, which(!is.na(t_c)))
  return(colonised)
}

sampleWithoutSurprises <- function(x) {
  if (length(x) <= 1) {
    return(x)
  } else {
    return(sample(x, 1))
  }
}

## Simulates epi data and returns in a format for MCMC
SimulateEpiDataMCMC <- function(N,L,LOS,p,z,beta_p,beta_h,num_col_hcw,patient_screening_interval,hcw_screening_interval) {
  #browser()
  outbreak <- simulateWard(N,L,LOS,p, z, 1e-15, beta_p, beta_h, num_col_hcw)
  screening_matrix <- screenData(outbreak, patient_screening_interval, hcw_screening_interval, z)
  
  outbreak$t_c[is.na(outbreak$t_c)] <- -1
  outbreak$source[is.na(outbreak$source)] <- -2
  screening_matrix[is.na(screening_matrix)] <- -1
  
  naive_coltimes <- rep(-1,length(outbreak$t_c))
  for(i in 1:length(outbreak$t_c)) {
    swab_results <- screening_matrix[i,]
    positive_days <- which(swab_results == 1) - 1
    if(length(positive_days)>0) {
      naive_coltimes[i] <- min(positive_days)
    }
  }
  #browser()
  naive_source <- rep(-2,length(outbreak$source))
  for(i in 1:length(outbreak$source)) {
    col_time <- naive_coltimes[i]
    if(col_time >= 0) {
      ## finite colonisation time, infer a source
      possible_infectors <- which(naive_coltimes < col_time & outbreak$t_d >= col_time & naive_coltimes != -1)
      if(length(possible_infectors) == 0) {
        naive_source[i] <- -1
      } else {
        naive_source[i] <- sampleWithoutSurprises(possible_infectors)
      }
    }
  }
  
  #browser()
  if(sum(1:length(outbreak$source) == naive_source) > 0) browser()
  
  out <- list("t_a" = outbreak$t_a,
              "t_c" = naive_coltimes,
              "t_d" = outbreak$t_d,
              "source" = naive_source,
              "screening_matrix" = screening_matrix,
              "true_coltimes" = outbreak$t_c,
              "true_source" = outbreak$source,
              "hcw_ind" = outbreak$hcw_ind)
  return(out)
}

## jc prob mutation
JC_prob_mutation <- function(rate,t) {
  prob <- 0.75*(1-exp(-4*rate*t))
  return(prob)
}



SimulateGeneticData_WHD <- function(epi_data, variant_numbers, average_import_distance, average_variant_distance, N, mutation_rate) {
  ## simulate genetic data
  positive_swabs <- which(epi_data$screening_matrix == 1, arr.ind = T)
  genetic_ids <- as.numeric(positive_swabs[,1])
  #browser()
  sample_times <- as.numeric(positive_swabs[,2])-1
  
  num_variants <- max(variant_numbers)
  
  
  ## now add everyone at the time of colonisation with each variant
  observed_length <- length(genetic_ids)
  ever_infected <- which(epi_data$true_coltimes != -1)
  for(i in 1:num_variants) {
    for(person in ever_infected) {
      gen_loc <- which(person == genetic_ids & i == variant_numbers & sample_times[i] == epi_data$true_coltimes[person])
      if(length(gen_loc)==0) {
        genetic_ids <- c(genetic_ids, person)
        sample_times <- c(sample_times, epi_data$true_coltimes[person])
        variant_numbers <- c(variant_numbers, i)
      }
    }
  }
  
  #browser()
  epi_data$true_source[epi_data$true_source > 0] <- epi_data$true_source[epi_data$true_source > 0] - 1
  temp_data_list <- list("genetic_ids" = genetic_ids-1,
                         "source" = epi_data$true_source,
                         "sample_times" = sample_times,
                         "variant_numbers" = variant_numbers,
                         "t_c" = epi_data$true_coltimes)
  gen_source <- ReturnGenSourceVector(temp_data_list)
  #browser()
  #length(which((gen_source == -1)))
  
  
  
  ## now simulate 
  
  full_distance_matrix <- matrix(NA, nrow=length(genetic_ids),ncol=length(genetic_ids))
  
  #browser()
  
  for(i in 1:length(gen_source)) {
    current_source <- gen_source[i] + 1
    full_distance_matrix[i,i] <- 0
    if(current_source == 0) {
      ## first introduction of the pathogen in this chain, compare to master distance
      
    } else {
      time_diff <- sample_times[i] - sample_times[current_source]
      if(time_diff != 0) {
        prob <- JC_prob_mutation(mutation_rate, time_diff)
        draw <- rbinom(1, N, prob)
      } else {
        draw <- 0
      }
      full_distance_matrix[current_source,i] <- draw
      full_distance_matrix[i,current_source] <- draw
    }
  }
  #browser()
  imports <- which(gen_source == -1)
  if(length(imports)>1) {
    for(i in 1:(length(imports)-1)) {
      for(j in (i+1):length(imports)) {
        variant_i <- variant_numbers[imports[i]]
        variant_j <- variant_numbers[imports[j]]
        if(variant_i == variant_j) {
          draw <- rpois(1,average_import_distance)
        } else {
          draw <- rpois(1,average_variant_distance)
        }
        
        full_distance_matrix[imports[j],imports[i]] <- draw
        full_distance_matrix[imports[i],imports[j]] <- draw
      }
    }
  }
  
  
  
  for(i in 1:nrow(full_distance_matrix)) 
  {
    for(j in 1:ncol(full_distance_matrix)) 
    {
      if(is.na(full_distance_matrix[i,j])) 
      {
        distance <- CalculateDistanceBetweenNodes(i,j,full_distance_matrix, gen_source)
        #print(c(i,j))
        if(is.null(distance)) return(NULL)
        full_distance_matrix[i,j] <- distance
      }
    }
  }
  
  out <- list("genetic_ids" = genetic_ids[1:observed_length],
              "sample_times" = sample_times[1:observed_length],
              "genetic_matrix" = full_distance_matrix[1:observed_length, 1:observed_length])
  
  
  return(out)
}
  
  
CalculateDistanceBetweenNodes <- function(i,j,distance_matrix,gen_source) {
  if(i==j) return(0)
  path1 <- ReturnPathToRoot(i,gen_source)
  path2 <- ReturnPathToRoot(j,gen_source)
  
  ## check if they are in the same pathway
  same_pathway <- F
  if(length(path1) > length(path2)) {
    ## look for j in path1
    j_loc <- which(path1==j)
    if(length(j_loc)>0) {
      same_pathway <- T
    }
  } else {
    ## look for i in path 2
    i_loc = which(path2==i)
    if(length(i_loc)>0) {
      same_pathway <- T
    }
  }
  #browser()
  if(same_pathway) {
    ## they are in the same pathway so just add nodes
    if(length(path1)>length(path2)) {
      j_loc = which(path1 == j)
      nodes_to_sum <- path1[j_loc:length(path1)]
    } else {
      i_loc = which(path2 == i)
      nodes_to_sum <- path2[i_loc:length(path2)]
    }
    
  } else {
    ## try and find a common node
    common_node <- ReturnCommonNode(path1,path2)
    if(!is.null(common_node)) {
      ## common nodes found
      common_node_path1_loc <- which(common_node == path1)
      common_node_path2_loc <- which(common_node == path2)
      nodes_to_sum <- c(rev(path1[(common_node_path1_loc+1):length(path1)]),path2[common_node_path2_loc:length(path2)])
    } else {
      ## completely separate pathway, sum all the way to the roots and then look at the distance between roots
      nodes_to_sum <- c(rev(path1),path2)
    }
  }
  dist <- 0
  for(i in 1:(length(nodes_to_sum)-1)) {
    node_i <- nodes_to_sum[i]
    node_j <- nodes_to_sum[i+1]
    dist <- dist + distance_matrix[node_i, node_j]
  }
  return(dist)
}
  
  
  
ReturnPathToRoot <- function(node, gen_source) {
  #browser()
  path <- c(node)
  source <- gen_source[node] + 1
  while(source != 0) {
    path <- c(source, path)
    source <- gen_source[source] + 1
  }
  return(path)
}

ReturnCommonNode <- function(path1, path2) {
  ## here the paths are backwards
  if(length(path2) > length(path1)) {
    temp <- path1
    path1 <- path2
    path2 <- temp
  }
  common_node <- NULL
  ## path 1 should be longer
  for(i in length(path2):1) {
    ## check all the nodes in the smaller path and determine if it is in the bigger path
    current_node = path2[i]
    if(current_node %in% path1) {
      common_node <- current_node
      break
    }
  }
  return(common_node)
}


ReturnVariantNumbersCutoff <- function(distance_matrix, cutoff) {
  adjacency_matrix <- matrix(NA, nrow=nrow(distance_matrix),ncol=ncol(distance_matrix))
  for(i in 1:nrow(distance_matrix)) {
    for(j in 1:nrow(distance_matrix)) {
      if(distance_matrix[i,j] < cutoff) {
        adjacency_matrix[i,j] <- 1
        adjacency_matrix[j,i] <- 1
      } else {
        adjacency_matrix[i,j] <- 0
        adjacency_matrix[j,i] <- 0
      }
    }
  }
  
  variant_list <- vector('list')
  k <- 1
  
  
  completed_nodes <- c()
  for(j in 1:nrow(distance_matrix)) {
    if(!(j %in% completed_nodes)) {
      node <- j
      all_nodes_found <- FALSE
      nodes_to_search <- node
      current_variants <- node
      while(!all_nodes_found) {
        new_nodes_to_search <- c()
        for(i in nodes_to_search) {
          node <- i
          children <- which(adjacency_matrix[node,]==1)
          children <- children[which(!(children %in% current_variants))]
          new_nodes_to_search <- c(new_nodes_to_search, children)
          current_variants <- c(current_variants, new_nodes_to_search)
        }
        nodes_to_search <- new_nodes_to_search
        if(length(nodes_to_search)==0) {
          variant_list[[k]] <- unique(current_variants)
          all_nodes_found <- TRUE
          completed_nodes <- c(completed_nodes, current_variants)
          k <- k + 1
        }
      }
    }
  }
  
  
  
  
  
  variant_numbers <- rep(NA,nrow(distance_matrix))
  for(i in 1:nrow(distance_matrix)) {
    variant_numbers[i] <- which(lapply(variant_list, `%in%`, x=i) == T)
  }
  return(variant_numbers)
  
}

## if epi_and_gen = 1, then epi and gen, otherwise just epi
PlotTraceEpiGenHCW <- function(results, burn_in, thin_factor, epi_and_gen) {
  if(epi_and_gen == 1) {
    res_full <- results$gen_res$res
    coltimes_full <- results$gen_res$coltimes
    source_full <- results$gen_res$source
  } else {
    res_full <- results$epi_res$res
    coltimes_full <- results$epi_res$coltimes
    source_full <- results$epi_res$source
  }

  true_z <- results$simulation_info$true_z
  true_p <- results$simulation_info$true_p
  true_beta_p <- results$simulation_info$true_beta_p
  true_beta_h <- results$simulation_info$true_beta_h
  true_lambda <- results$simulation_info$true_lambda
  true_mu <- results$simulation_info$true_mu
  
  
  if(burn_in >0 ) {
    res_after_burn <- res_full[-c(1:burn_in),]
    coltimes_after_burn <- coltimes_full$V1[-c(1:burn_in)]
    source_after_burn <- source_full[-c(1:burn_in),]
  } else {
    res_after_burn <- res_full
    coltimes_after_burn <- coltimes_full
    source_after_burn <- source_full
  }
  
  
  if(thin_factor > 0) {
    res <- res_after_burn[seq(from=1,to=nrow(res_after_burn),by=thin_factor),]
    coltimes <- coltimes_after_burn[seq(from=1,to=length(coltimes_after_burn),by=thin_factor)]
    source <- source_after_burn[seq(from=1,to=nrow(source_after_burn),by=thin_factor),]
  } else {
    res <- res_after_burn
    coltimes <- coltimes_after_burn
    source <- source_after_burn
  }
  
  true_coltime_sum <- sum(current_res$simulation_info$epi_data$true_coltimes[current_res$simulation_info$epi_data$true_coltimes != -1])
  
  dev.new()
  par(mfrow=c(4,2))
  
  plot(res$V1,type="l",main="Trace plot for z")
  abline(h=true_z, col=2)
  abline(h=quantile(res$V1,c(0.025,0.975)),col="blue")
  
  plot(res$V2,type="l", main="Trace plot for p")
  abline(h=true_p, col=2)
  abline(h=quantile(res$V2,c(0.025,0.975)),col="blue")
  
  plot(res$V3,type="l", main="Trace plot for beta_p")
  abline(h=true_beta_p, col=2)
  abline(h=quantile(res$V3,c(0.025,0.975)),col="blue")
  
  plot(res$V4, type="l", main="Trace plot for beta_h")
  abline(h=true_beta_h, col=2)
  abline(h=quantile(res$V4,c(0.025,0.975)),col="blue")
  
  plot(res$V5, type="l", main="Trace plot for lambda")
  abline(h=true_lambda, col=2)
  abline(h=quantile(res$V5,c(0.025,0.975)),col="blue")
  
  plot(res$V6, type="l", main="Trace plot for mu")
  abline(h=true_mu, col=2)
  abline(h=quantile(res$V6,c(0.025,0.975)),col="blue")
  
  plot(coltimes, type="l", main="Trace plot for sum of colonisation times")
  abline(h=true_coltime_sum,col=2)
  abline(h=quantile(coltimes,c(0.025,0.975)),col="blue")
  
  par(mfrow=c(1,1))
}







# Here set inferred = 1 if we return the inferred chains, set to 0 to return the true chains
ReturnInferredTransmissionChains <- function(epi_data, source_table) {
  # Identify transmission clusters
  importations <- source_table$id[which(source_table$inferred_source==-1)]
  hcw_id <- which(epi_data$hcw_ind==1)
  importations <- c(importations, hcw_id[sapply(which(epi_data$hcw_ind==1), function(x) x %in% source_table$inferred_source)])
  transmission_chains <- vector('list',length(importations))
  #browser()
  for(i in 1:length(importations)) {
    current_importation <- importations[i]
    
    chain_found <- FALSE
    primary_targets <- current_importation
    chain <- c()
    while(!chain_found) {
      secondary_targets <- c()
      for(j in primary_targets) {
        chain <- c(chain, j)
        current_infections <- which(source_table$id==j)
        if(length(current_infections)>0) {
          secondary_targets <- c(secondary_targets, source_table$id[which(source_table$inferred_source==j)])
        }
        
      }
      primary_targets <- secondary_targets
      if(length(secondary_targets) == 0) {
        chain_found <- TRUE
      }
    }
    transmission_chains[[i]] <- chain
    
  }
  return(transmission_chains)
}

ReturnTrueTransmissionChains <- function(epi_data) {
  #browser()
  importations <- which(epi_data$true_source == -1 & epi_data$hcw_ind == 0)
  observed_infected <- which(epi_data$t_c != -1 & epi_data$hcw_ind == 0)
  
  transmission_chains <- vector('list',length(importations))
  for(i in 1:length(importations)) {
    current_importation <- importations[i]
    
    chain_found <- FALSE
    primary_targets <- current_importation
    chain <- c()
    while(!chain_found) {
      secondary_targets <- c()
      for(j in primary_targets) {
        chain <- c(chain, j)
        current_infections <- which(epi_data$true_source==j)
        if(length(current_infections)>0) {
          secondary_targets <- c(secondary_targets, current_infections)
        }
        
      }
      primary_targets <- secondary_targets
      if(length(secondary_targets) == 0) {
        chain_found <- TRUE
      }
    }
    transmission_chains[[i]] <- chain
    
  }
  return(transmission_chains)
}



# Calculates the proportion of clusters correctly identified
CalculateSourceClusterProportion <- function(epi_data, source) {
  
  ever_infected <- which(epi_data$t_c != -1 & epi_data$hcw_ind == 0)
  source_table <- data.frame("id" = numeric(length(ever_infected)),
                             "inferred_source" = numeric(length(ever_infected)),
                             "true_source" = numeric(length(ever_infected)),
                             "posterior_probability" = numeric(length(ever_infected)))
  
  
  
  for(i in 1:length(ever_infected)) {
    current_target <- ever_infected[i]
    posterior_source_distribution <- table(source[,current_target])
    most_likely_source <- as.numeric(names(posterior_source_distribution)[which.max(posterior_source_distribution)])
    if(most_likely_source >=0 ) {
      most_likely_source <- most_likely_source + 1
    }
    source_table[i,1] <- current_target
    source_table[i,2] <- most_likely_source
    source_table[i,3] <- epi_data$true_source[current_target]
    source_table[i,4] <- as.numeric(max(table(source[,current_target]))/sum(table(source[,current_target])))
  }
  source_proportion <- sum(source_table[,2]==source_table[,3])/nrow(source_table)
  
  inferred_chains <- ReturnInferredTransmissionChains(epi_data, source_table)
  inferred_chains <- lapply(inferred_chains, sort)
  #browser()
  true_chains <- ReturnTrueTransmissionChains(epi_data)
  true_observed_chains <- sapply(true_chains, function(x) x[x %in% ever_infected])
  true_observed_chains <- sapply(true_observed_chains, sort)
  
  #browser()
  # correct_chains <- 0
  # for(i in 1:length(inferred_chains)) {
  #   browser()
  #   current_inferred_chain <- inferred_chains[[i]]
  #   for(j in 1:length(true_observed_chains)) {
  #     current_true_observed_chain <- true_observed_chains[[j]]
  #     counter <- 0
  #     for(ii in current_inferred_chain) {
  #       if(ii %in% current_true_observed_chain) {
  #         counter <- counter + 1
  #       }
  #     }
  #     if(counter == length(current_inferred_chain)) {
  #       correct_chains <- correct_chains + 1
  #     }
  #   }
  # }
  correct_chains <- 0
  for(i in 1:length(inferred_chains)) {
    current_inferred_chain <- inferred_chains[[i]]
    if(sum(true_observed_chains %in% list(current_inferred_chain))>0) {
      correct_chains <- correct_chains + 1
    }
  }
  
  
  
  
  
  cluster_proportion <- correct_chains / length(true_observed_chains)
  
  out <- list("source_table" = source_table,
              "source_propotion" = source_proportion,
              "true_chains" = true_observed_chains,
              "inferred_chains" = inferred_chains,
              "chain_proportion" = cluster_proportion)
  return(out)
}

PlotTransmissionNetwork <- function(t_a, t_c, t_d, source, genetic_ids, sample_times) {
  ever_infected <- which(t_c != -1)
  browser()
  source[source>=0] <- source[source>=0] + 1
  genetic_ids <- genetic_ids+1

  
  plot(0,type="n", ylim=c(0,length(ever_infected)), xlim=c(0,max(t_d)))
  for(i in 1:length(ever_infected)) {
    current_person <- ever_infected[i]
    current_source <- source[current_person]
    segments(t_a[current_person],i,t_c[current_person],i,lwd=2)
    segments(t_c[current_person],i,t_d[current_person],i,col=2, lwd=2)
    if(current_source != -1) {
      source_loc <- which(current_source == ever_infected)
      if(length(source_loc)==1) {
        arrows(t_c[current_person], source_loc, t_c[current_person], i)
      } else {
        browser()
      }
      
    }
    text(t_a[current_person]-3,i,current_person)
  }
  
  for(i in 1:length(genetic_ids)) {
    current_gen_id <- genetic_ids[i]
    current_sample_time <- sample_times[i]
    text(current_sample_time, which(current_gen_id==ever_infected)+0.4,i,cex=0.7)
  }
}



## tree plotting

PlotNetworkFromString <- function(string) {
  split <- strsplit(string,"\n")[[1]]
  
  
  data.list <- vector('list',6)
  for(i in 1:length(split)) {
    data.list[[i]] <- as.numeric(strsplit(split[i], " ")[[1]])
  }
  
  #browser()
  PlotTransmissionNetwork(data.list[[1]],data.list[[2]],data.list[[3]],data.list[[4]],data.list[[5]],data.list[[6]])
}

CalculateTrueTransmissionProportion <- function(epi_data, source, p_L) {
  ever_infected <- which(epi_data$t_c != -1)
  correct_source <- 0
  for(current_person in ever_infected) {
    true_source <- epi_data$true_source[current_person]
    posterior_source_distribution <- table(source[,current_person])/sum(table(source[,current_person]))
    true_source_loc <- which(names(posterior_source_distribution)==true_source)
    if(length(true_source_loc)>0) {
      posterior_probability <- posterior_source_distribution[true_source_loc]
      if(posterior_probability > p_L) {
        correct_source <- correct_source + 1
      }
    }
  }
  return(correct_source/length(ever_infected))
}




  