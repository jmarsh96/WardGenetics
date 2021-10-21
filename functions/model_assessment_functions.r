# simulates a ward according to the observed admission and discharge times provided


simulateWard_assessment <- function(epi_data,p=0.05, z=0.95, a0 = 0.000000001, a1 = 0.005, a2=0.001) {
  
  N <- length(epi_data$t_a)
  t_a <- epi_data$t_a
  t_c <- rep(NA,N)
  t_d <- epi_data$t_d
  source <- rep(NA,N)
  #browser()
  hcw_indicator <- epi_data$hcw_ind
  
  

  # sample importations
  positive_on_admission <- sample(c(TRUE,FALSE),N,replace=T,prob=c(p,1-p))
  t_c[positive_on_admission] <- t_a[positive_on_admission]
  source[positive_on_admission] <- -1
  
  maxD <- max(t_d)
  
  t_c[hcw_indicator] <- epi_data$t_c[hcw_indicator]
  
  
  
  #browser()
  #go through each day
  for(i in 1:maxD) {
    #if(i == 28) browser()
    susceptibles <- countSusceptible(t_a[1:N],t_c[1:N],t_d[1:N],i)
    colonised <- countColonisedR(t_a[1:N],t_c[1:N],t_d[1:N],i)
    if(sum(hcw_indicator) > 0) {
      colonised_hcw <- which(t_c[hcw_indicator] < i)
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


screen_data_assessment <- function(epi_data, screening_matrix, hcw_screening_interval, z) {
  hcw_indicator <- epi_data$hcw_indicator
  N <- length(epi_data$t_a)
  maxD <- max(epi_data$t_d)
  results_matrix <- matrix(NA, nrow=N, ncol=maxD+1)
  #browser()
  for(i in 1:N) {
    if(hcw_indicator[i]==0) {
      test_days <- test_days <- which(screening_matrix[i,] != -1) - 1
    } else {
      test_days <- seq(from=epi_data$t_a[i], to=epi_data$t_d[i], by=hcw_screening_interval)
    }
    
    
    if(is.na(epi_data$t_c[i])){
      results_matrix[i,test_days+1] <- 0
    } else {
      for(t in test_days) {
        if(t >= epi_data$t_c[i]) {
          results_matrix[i,t+1] <- rbinom(1,1,z)
        } else {
          results_matrix[i,t+1] <- 0
        }
      }
    }
  }
  return(results_matrix)
}


## Simulate a genetic matrix given the genetic source, sample times and mutation rate
genetic_matrix_assessment <- function(gen_source, sample_times, mutation_rate) {
  
  
  full_distance_matrix <- matrix(NA,nrow=length(gen_source),ncol=length(gen_source))
  
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
  
  return(full_distance_matrix)
}


ReturnGeneticTree <- function(gen_source) {
  
  require(igraph)
  #browser()
  import_loc <- which(gen_source==-1)
  gen_source <- gen_source + 1
  gen_source[import_loc] <- import_loc
  #browser()
  g <- make_graph(c(rbind(gen_source,(1:(length(gen_source))))), n=length(gen_source))
  
  #plot(simplify(g), vertex.size = 8, vertex.label = genetic_ids, directed=T, edge.size=0.5, edge.arrow.size=0.3, edge.arrow.width=0.5)
  plot(simplify(g), vertex.size = 8, directed=T, edge.size=0.5, edge.arrow.size=0.3, edge.arrow.width=0.5, vertex.label = 0:(length(gen_source)-1))
  
}







