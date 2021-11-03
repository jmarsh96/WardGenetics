# simulates a ward according to the observed admission and discharge times provided


simulateWard_assessment <- function(epi_data,p=0.05, z=0.95, a0 = 0.000000001, a1 = 0.005, a2=0.001) {
  
  hcw_indicator <- epi_data$hcw_ind
  N <- sum(hcw_indicator==0)
  t_a <- epi_data$t_a
  t_c <- rep(NA,length(t_a))
  t_d <- epi_data$t_d
  source <- rep(NA,length(t_a))
  #browser()
  
  source[hcw_indicator==1] <- -1
  
  

  # sample importations
  positive_on_admission <- sample(c(TRUE,FALSE),N,replace=T,prob=c(p,1-p))
  t_c[positive_on_admission] <- t_a[positive_on_admission]
  source[positive_on_admission] <- -1
  
  maxD <- max(t_d)
  
  #browser()
  t_c[hcw_indicator==1] <- epi_data$t_c[hcw_indicator==1]
  
  
  
  #browser()
  #go through each day
  for(i in 1:maxD) {
    #if(i == 28) browser()
    susceptibles <- countSusceptible(t_a[hcw_indicator==0],t_c[hcw_indicator==0],t_d[hcw_indicator==0],i)
    colonised <- countColonisedR(t_a[hcw_indicator==0],t_c[hcw_indicator==0],t_d[hcw_indicator==0],i)
    if(sum(hcw_indicator) > 0) {
      colonised_hcw <- which(t_c[hcw_indicator==1] < i)
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


screen_data_assessment <- function(new_epi_data, epi_data, screening_matrix, hcw_screening_interval, z) {
  #browser()
  hcw_indicator <- epi_data$hcw_ind
  N <- length(epi_data$t_a)
  maxD <- max(epi_data$t_d)
  results_matrix <- matrix(NA, nrow=N, ncol=maxD+1)
  #browser()
  for(i in 1:N) {
    if(hcw_indicator[i]==0) {
      test_days <- which(screening_matrix[i,] != -1) - 1
    } else {
      test_days <- seq(from=epi_data$t_a[i], to=epi_data$t_d[i], by=hcw_screening_interval)
    }
    
    
    if(is.na(new_epi_data$t_c[i])){
      results_matrix[i,test_days+1] <- 0
    } else {
      for(t in test_days) {
        if(t >= new_epi_data$t_c[i]) {
          results_matrix[i,t+1] <- rbinom(1,1,z)
        } else {
          results_matrix[i,t+1] <- 0
        }
      }
    }
  }
  return(results_matrix)
}

SimulateEpiData_ModelAssessment <- function(epi_data, z, p, beta_p, beta_h, hcw_screening_interval) {
  epi <- simulateWard_assessment(epi_data, p, z, 1e-9, beta_p, beta_h)
  #browser()
  screening_matrix <- screen_data_assessment(epi, epi_data, epi_data$screening_matrix, hcw_screening_interval, z)
  
  #browser()
  out <- list("t_a" = epi$t_a,
              "t_c" = epi$t_c,
              "t_d" = epi$t_d,
              "importations" = epi$importations,
              "source" = epi$source,
              "hcw_ind" = epi$hcw_ind,
              "screening_matrix" = screening_matrix)
  return(out)
}


CalculateCt_fromScreenMatrix <- function(epi_data) {
  num_patients <- sum(epi_data$hcw_ind==0)
  inferred_col_times <- rep(NA,num_patients)
  for(i in 1:num_patients) {
    positive_days <- which(epi_data$screening_matrix[i,] == 1) - 1
    if(length(positive_days)>0) {
      inferred_col_times[i] <- positive_days[1]
    }
  }
  
  maxD <- max(epi_data$t_d)
  C <- sapply(0:maxD, function(x) CalculateNumColonisedOnDay(inferred_col_times, epi_data$t_d, epi_data$hcw_ind, x))
  return(C)
}

CalculateNumberPatientsToEverHavePositiveSwab <- function(epi_data) {
  return(sum(apply(epi_data$screening_matrix, 1, function(x) ifelse(sum(x == 1, na.rm=T) > 0,T,F))))
}

CalculateNumColonisedOnDay <- function(t_c, t_d, hcw_ind, t) {
  #browser()
  t_d <- t_d[!hcw_ind]
  num_colonised <- sum(t_c < t & t_d > t & !is.na(t_c) & t_c != -1)
  return(num_colonised)
}


## Simulate a genetic matrix given the genetic source, sample times and mutation rate
genetic_matrix_assessment <- function(N, gen_source, sample_times, variant_numbers, mutation_rate, average_import_distance, num_observed_sequences) {
  
  
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
  

  
  import_idx <- which(gen_source == -1)
  for(i in 1:length(import_idx)) {
    for(j in 1:length(import_idx)) {
      cur_i <- import_idx[i]
      cur_j <- import_idx[j]
      if(cur_i != cur_j) {
        if(variant_numbers[cur_i] == variant_numbers[cur_j]) {
          full_distance_matrix[cur_i,cur_j] <- rpois(1, average_import_distance)
        }
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
  
  return(full_distance_matrix[1:num_observed_sequences,1:num_observed_sequences])
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


CalculateTotalEuclidianDistance <- function(x_points, y_points) {
  total_distance <- 0
  k <- length(x_points)
  for(j in 2:k) {
    for(i in 1:j) {
      #browser()
      total_distance <- total_distance + sqrt((x_points[i]-y_points[j])^2)
      if(is.na(total_distance)) browser()
    }
  }
  return(total_distance)
}

genetic_matrix_assessment_marginal_calculation <- function(true_matrix, simulated_matrices, variant_numbers) {
  #browser()
  k <- nrow(true_matrix)
  p_value_matrix <- matrix(NA,nrow=k, ncol=k)
  for(j in 2:k) {
    for(i in 1:j) {
      simulated_values <- sapply(simulated_matrices, function(x) x[i,j])
      true_value <- true_matrix[i,j]
      p_value_matrix[i,j] <- sum(simulated_values > true_value)
    }
  }
  #browser()
  diag(p_value_matrix) <- rep(NA,length(diag(p_value_matrix)))
  ordered_variants <- order(variant_numbers)
  p_value_matrix <- p_value_matrix[ordered_variants,ordered_variants]
  return(p_value_matrix/length(simulated_matrices))
}



GeneticModelAssessmentSummaries <- function(true_matrix, distance_matrices) {
  par(mfrow=c(2,2))
  upper_tri_sums <- sapply(distance_matrices, function(x) sum(x[upper.tri(x)]))
  hist(upper_tri_sums)
  abline(v=sum(true_matrix[upper.tri(true_matrix)]),col=2)
  
  matrix_means <- sapply(distance_matrices, function(x) mean(x[upper.tri(x)]))
  hist(matrix_means)
  abline(v=mean(true_matrix[upper.tri(true_matrix)]),col=2)
  
  matrix_medians <- sapply(distance_matrices, function(x) median(x[upper.tri(x)]))
  hist(matrix_medians)
  abline(v=median(true_matrix[upper.tri(true_matrix)]),col=2)
  
  matrix_one_norms <- sapply(distance_matrices, function(x) norm(x))
  hist(matrix_one_norms)
  abline(v=norm(true_matrix),col=2)
  
  matrix_frobenieus_norms <- sapply(distance_matrices, function(x) norm(x,type="f"))
  hist(matrix_frobenieus_norms)
  abline(v=norm(true_matrix,type="f"),col=2)
  
  matrix_m_norms <- sapply(distance_matrices, function(x) norm(x,type="m"))
  hist(matrix_m_norms)
  abline(v=norm(true_matrix,type="m"),col=2)
  
  matrix_2_norms <- sapply(distance_matrices, function(x) norm(x,type="2"))
  hist(matrix_2_norms)
  abline(v=norm(true_matrix,type="2"),col=2)
  
  GeneticMatrix_MDS(true_matrix, distance_matrices)
  par(mfrow=c(1,1))
  
}

GeneticMatrix_MDS <- function(true_matrix, distance_matrices) {
  
  d <- as.dist(true_matrix) # euclidean distances between the rows
  fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
  fit # view results
  
  # plot solution
  x <- fit$points[,1]
  y <- fit$points[,2]
  #plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
  #     main="Metric MDS", type="n")
  #text(x, y, cex=.7) 
  
  max_iter <- length(distance_matrices)
  ## perform MDS for the simulated distance matrices
  x_points <- vector('list', max_iter)
  y_points <- vector('list', max_iter)
  
  for(i in 1:max_iter) {
    cur_dist_matrix <- distance_matrices[[i]]
    d <- as.dist(cur_dist_matrix) # euclidean distances between the rows
    fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
    fit # view results
    
    # plot solution
    x_points[[i]] <- fit$points[,1]
    y_points[[i]] <- fit$points[,2]
  }
  
  
  total_distances <- rep(NA,max_iter)
  for(i in 1:max_iter) {
    total_distances[i] <- CalculateTotalEuclidianDistance(x_points[[i]],y_points[[i]])
  }
  hist(total_distances)
  abline(v=CalculateTotalEuclidianDistance(x,y),col=2)
}


