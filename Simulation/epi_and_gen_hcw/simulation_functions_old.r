# simulates a ward of N individuals over D days
# length_of_stay is the average time in the ward
# p is the probability of being colonised on admission
# z is the test sensitivity 
# a0 is the background rate of infection
# a1 is the contact rate with other infectives on the ward
# lambda is the JC mutation rate
# screening_interval is the number of days a screen is performed

simulateWard <- function(N,D,length_of_stay=7,p=0.05, z=0.95, a0 = 0.05, a1 = 0.000001, mutation_rate = 3.3*10^(-6)/365) {
  t_a <- numeric(N)
  t_c <- rep(NA,N)
  t_d <- numeric(N)
  source <- rep(NA,N)
  
  # generate admission times
  t_a <- sort(sample(1:D,N,replace=T))
  t_a <- t_a - min(t_a)
  
  # sample importations
  positive_on_admission <- sample(c(TRUE,FALSE),N,replace=T,prob=c(p,1-p))
  t_c[positive_on_admission] <- t_a[positive_on_admission]
  source[positive_on_admission] <- -1
  
  # generate discharge times
  t_d <- t_a + rpois(N,length_of_stay)
  
  maxD <- max(t_d)
  #browser()
  #go through each day
  for(i in 1:maxD) {
    susceptibles <- countSusceptible(t_a,t_c,t_d,i)
    colonised <- countColonisedR(t_a,t_c,t_d,i)
    current_rate <- a0 + length(colonised)*a1
    for(j in susceptibles) {
      if(runif(1)<(1-exp(-1*current_rate))) {
        # colonisation has occured
        t_c[j] <- i
        
        #choose a source
        if(runif(1)<(a0/current_rate)) {
          #background infection
          source[j] <- j
        } else {
          source[j] <- sampleWithoutSurprises(colonised)
        }
      }
    }
    
  }
  importations <- rep(0,N)
  importations[positive_on_admission] <- 1
  

  
  out <- list("t_a" = t_a, "t_c" = t_c, "t_d" = t_d, "importations" = importations, "source" = source)
  return(out)
}






# simulates a ward of N individuals over D days with healthcare workers included
# length_of_stay is the average time in the ward
# p is the probability of being colonised on admission
# z is the test sensitivity 
# a0 is the background rate of infection
# a1 is the contact rate with other infectives on the ward
# a2 is the contact rate with health care workers
# num_hcw is the number of healthcare workers on the ward
# col_hcw_prop is the proportion of healthcare workers that are colonised

simulateWard_hcw <- function(N,D,length_of_stay=7,p=0.05, z=0.95, a0 = 0.005, a1 = 0.000001, a2=0.001, num_hcw, col_hcw_prop=0.3) {
  t_a <- numeric(N+num_hcw)
  t_c <- rep(NA,N+num_hcw)
  t_d <- numeric(N+num_hcw)
  source <- rep(NA,N)
  
  HCW_indicator <- rep(0,N+num_hcw)
  HCW_indicator[(N+1):(N+num_hcw)] <- 1
  
  
  # generate admission times
  t_a[1:N] <- sort(sample(1:D,N,replace=T))
  t_a[1:N] <- t_a[1:N] - min(t_a)
  
  t_a[(N+1):(N+num_hcw)] <- 0
  
  # sample importations
  positive_on_admission <- sample(c(TRUE,FALSE),N,replace=T,prob=c(p,1-p))
  t_c[positive_on_admission] <- t_a[positive_on_admission]
  source[positive_on_admission] <- -1
  
  # generate discharge times
  t_d[1:N] <- t_a[1:N] + rpois(N,length_of_stay)
  maxD <- max(t_d)
  
  
  colonised_hcw <- sample(num_hcw, floor(num_hcw*col_hcw_prop))
  t_c[N+colonised_hcw] <- sample(0:maxD,size=length(colonised_hcw))
  

  #browser()
  #go through each day
  for(i in 1:maxD) {
    susceptibles <- countSusceptible(t_a[1:N],t_c[1:N],t_d[1:N],i)
    colonised <- countColonisedR(t_a[1:N],t_c[1:N],t_d[1:N],i)
    colonised_hcw <- which(t_c[(N+1):(N+num_hcw)] < i)
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
          source[j] <- sampleWithoutSurprises(N+colonised_hcw)
        }
      }
    }
    
  }
  importations <- rep(0,N)
  importations[positive_on_admission] <- 1
  
  
  
  out <- list("t_a" = t_a, "t_c" = t_c, "t_d" = t_d, "importations" = importations, "source" = source)
  return(out)
}



# simulates a ward of N individuals over D days with healthcare workers included
# length_of_stay is the average time in the ward
# p is the probability of being colonised on admission
# z is the test sensitivity 
# a0 is the background rate of infection
# a1 is the contact rate with other infectives on the ward
# a2 is the contact rate with health care workers
# num_hcw is the number of healthcare workers on the ward
# col_hcw_prop is the proportion of healthcare workers that are colonised

## this function has separate vectors for HCW and patients
simulateWard_hcw2<- function(N,D,length_of_stay=7,p=0.05, z=0.95, a0 = 0.005, a1 = 0.000001, a2=0.001, num_hcw, col_hcw_prop=0.3) {
  t_a <- numeric(N)
  t_c <- rep(NA,N)
  t_d <- numeric(N)
  source <- rep(NA,N)
  
  
  # generate admission times
  t_a <- sort(sample(1:D,N,replace=T))
  t_a <- t_a - min(t_a)
  
  # sample importations
  positive_on_admission <- sample(c(TRUE,FALSE),N,replace=T,prob=c(p,1-p))
  t_c[positive_on_admission] <- t_a[positive_on_admission]
  source[positive_on_admission] <- -1
  
  # generate discharge times
  t_d <- t_a + rpois(N,length_of_stay)
  maxD <- max(t_d)
  
  ## generate hcw data
  t_a_hcw <- numeric(num_hcw)
  t_d_hcw <- rep(maxD,num_hcw)
  t_c_hcw <- rep(NA,num_hcw)
  
  colonised_hcw <- sample(num_hcw, floor(num_hcw*col_hcw_prop))
  t_c_hcw[colonised_hcw] <- sample(0:maxD,size=length(colonised_hcw))
  
  
  #browser()
  #go through each day
  for(i in 1:maxD) {
    susceptibles <- countSusceptible(t_a,t_c,t_d,i)
    colonised <- countColonisedR(t_a,t_c,t_d,i)
    colonised_hcw <- which(t_c_hcw < i)
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
          source[j] <- paste("H",sampleWithoutSurprises(colonised_hcw), sep="")
        }
      }
    }
    
  }
  importations <- rep(0,N)
  importations[positive_on_admission] <- 1
  
  
  
  out <- list("t_a" = t_a, "t_c" = t_c, "t_d" = t_d, "importations" = importations, "source" = source)
  return(out)
}


simulateWard_old <- function(N,D,length_of_stay=7,p=0.05, z=0.95, a0 = 0.05, a1 = 0.01, mutation_rate = 3.3*10^(-6)/365) {
  t_a <- numeric(N)
  t_c <- rep(NA,N)
  t_d <- numeric(N)
  source <- rep(NA,N)
  pID <- 1:N
  
  # generate admission times
  t_a <- sort(sample(1:D,N,replace=T))
  t_a <- t_a - min(t_a)
  
  # sample importations
  positive_on_admission <- sample(c(TRUE,FALSE),N,replace=T,prob=c(p,1-p))
  t_c[positive_on_admission] <- t_a[positive_on_admission]
  source[positive_on_admission] <- -1
  
  # generate discharge times
  t_d <- t_a + rpois(N,length_of_stay)
  
  maxD <- max(t_d)
  #go through each day
  for(i in 1:maxD) {
    susceptibles <- countSusceptible(t_a,t_c,t_d,i)
    colonised <- countColonisedR(t_a,t_c,t_d,i)
    current_rate <- a0 + length(colonised)*a1
    for(j in susceptibles) {
      if(runif(1)<(1-exp(-1*current_rate))) {
        # colonisation has occured
        t_c[j] <- i
        
        #choose a source
        if(runif(1)<(a0/current_rate)) {
          #background infection
          source[j] <- j
        } else {
          source[j] <- sampleWithoutSurprises(colonised)
        }
      }
    }
    
  }
  importations <- rep(0,N)
  importations[positive_on_admission] <- 1
  
  
  
  out <- list("t_a" = t_a, "t_c" = t_c, "t_d" = t_d, "importations" = importations, "source" = source)
  return(out)
}

#### simulate ward with only 1 importation randomly selected
simulateWard2 <- function(N,D,length_of_stay=7, z=0.95, a0 = 0.05, a1 = 0.01, mutation_rate = 3.3*10^(-6)/365) {
  t_a <- numeric(N)
  t_c <- rep(NA,N)
  t_d <- numeric(N)
  source <- rep(NA,N)
  pID <- 1:N
  
  # generate admission times
  t_a <- sort(sample(1:D,N,replace=T))
  t_a <- t_a - min(t_a)
  
  # sample importations
  positive_on_admission <- sample(N,1)
  t_c[positive_on_admission] <- t_a[positive_on_admission]
  source[positive_on_admission] <- -1
  
  # generate discharge times
  t_d <- t_a + rpois(N,length_of_stay)
  
  maxD <- max(D)
  #go through each day
  for(i in 1:maxD) {
    susceptibles <- countSusceptible(t_a,t_c,t_d,i)
    colonised <- countColonisedR(t_a,t_c,t_d,i)
    current_rate <- a0 + length(colonised)*a1
    for(j in susceptibles) {
      if(runif(1)<(1-exp(-1*current_rate))) {
        # colonisation has occured
        t_c[j] <- i
        
        #choose a source
        if(runif(1)<(a0/current_rate)) {
          #background infection
          source[j] <- j
        } else {
          source[j] <- sampleWithoutSurprises(colonised)
        }
      }
    }
    
  }
  importations <- rep(0,N)
  importations[positive_on_admission] <- 1
  
  
  
  out <- list("t_a" = t_a, "t_c" = t_c, "t_d" = t_d, "importations" = importations, "source" = source)
  return(out)
}


## jc prob mutation
JC_prob_mutation <- function(rate,t) {
  prob <- 0.75*(1-exp(-4*rate*t))
  return(prob)
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

IsEveryoneAbleToBeInfected <- function(t_c, t_d, importations) {
  colonised <- which(t_c != -1)
  for(i in colonised) {
    
    able_to_infect <- which(t_c < t_c[i] & t_c != -1 & t_d >= t_c[i])
    if(length(able_to_infect)==0) {
      # check if an import can affect
      import_able_to_infect <- which(importations == 1 & t_c == t_c[i])
      if(length(import_able_to_infect)==0) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

returnTransmissionNetwork <- function(data) {
  if(!require(igraph)) install.packages('igraph')
  library(igraph)
  sources <- data$source[which(!is.na(data$t_c))]
  IDs <- which(!is.na(data$t_c))
  importations <- which(sources == -1)
  
  
  links <- data.frame("source" = sources[-importations], "target" = IDs[-importations])
  nodes <- data.frame("nodes" = IDs)
  net <- graph_from_data_frame(d=links,vertices=nodes)
  return(net)
}

## Return the true transmission network of simulated data
ReturnTransmissionNetwork <- function(epi_data)
{
  library(shape)
  use_truth <- T
  if(use_truth) {
    ever_infected <- which(epi_data$true_coltimes != -1)
  } else {
    ever_infected <- which(epi_data$t_c != -1)
  }
  source <- epi_data$true_source[ever_infected]
  N <- length(ever_infected)
  plot(0,type="n", xlim=c(-2,2),ylim=c(-2,2))
  for(i in 1:N) {
    if(source[i] == -1) {
      points(cos(2*pi*i/N),sin(2*pi*i/N),pch=16,col=2)
    } else {
      points(cos(2*pi*i/N),sin(2*pi*i/N),pch=16)
      source_loc <- which(source[i]==ever_infected)
      Arrows(cos(2*pi*source_loc/N),sin(2*pi*source_loc/N),cos(2*pi*i/N),sin(2*pi*i/N), 
             arr.type = "triangle",  arr.width = 0.1, arr.length = 0.1, arr.adj = 1, col="blue", lwd = 2)
    }
    text(cos(2*pi*i/N),sin(2*pi*i/N)+0.15, ever_infected[i])
  }
}



ReturnInferredTransmissionNetwork <- function(epi_data, source) {
  use_truth <- F
  if(use_truth) {
    ever_infected <- which(epi_data$true_coltimes != -1)
  } else {
    ever_infected <- which(epi_data$t_c != -1)
  }
  true_source <- epi_data$true_source[ever_infected]
  N <- length(ever_infected)
  plot(0,type="n", xlim=c(-2,2),ylim=c(-2,2))
  for(i in 1:N) {
    inferred_sources <- table(source[,ever_infected[i]])
    inferred_post_probability <- inferred_sources/sum(inferred_sources)
    import_loc <- which(names(inferred_post_probability) == "-1")
    if(inferred_post_probability[import_loc] > 0.5 && length(import_loc) != 0) {
      # they are inferred to be an importation
      points(cos(2*pi*i/N),sin(2*pi*i/N),pch=16,col=2)
    } else {
      # they are inferred to be an acquisition
      points(cos(2*pi*i/N),sin(2*pi*i/N),pch=16)
      for(j in 1:length(inferred_post_probability)) {
        current_probability <- inferred_post_probability[j]
        current_source <- as.numeric(names(inferred_post_probability)[j]) + 1
        if(current_source != 0) {
          source_loc <- which(current_source==ever_infected)
          if(current_probability > 0.75) {
            Arrows(cos(2*pi*source_loc/N),sin(2*pi*source_loc/N),cos(2*pi*i/N),sin(2*pi*i/N), 
                   arr.type = "triangle",  arr.width = 0.1, arr.length = 0.1, arr.adj = 1, lwd = 2,
                   col = rgb(0,0,1,1))
          } else if(current_probability > 0.5 && current_probability <= 0.75) {
            Arrows(cos(2*pi*source_loc/N),sin(2*pi*source_loc/N),cos(2*pi*i/N),sin(2*pi*i/N), 
                   arr.type = "triangle",  arr.width = 0.1, arr.length = 0.1, arr.adj = 1, lwd = 2,
                   col = rgb(0,0,1,0.75))
          } else if(current_probability > 0.25 && current_probability <= 0.50) {
            Arrows(cos(2*pi*source_loc/N),sin(2*pi*source_loc/N),cos(2*pi*i/N),sin(2*pi*i/N), 
                   arr.type = "triangle",  arr.width = 0.1, arr.length = 0.1, arr.adj = 1, lwd = 2,
                   col = rgb(0,0,1,0.50))
          }
        }
      }
    }
    text(cos(2*pi*i/N),sin(2*pi*i/N)+0.15, ever_infected[i])
  }
}



screenData <- function(data, screening_interval, z) {
  N <- length(data$t_a)
  maxD <- max(data$t_d)
  results_matrix <- matrix(NA, nrow=N, ncol=maxD+1)
  #browser()
  for(i in 1:N) {
    test_days <- seq(from=data$t_a[i], to=data$t_d[i], by=screening_interval)
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

screenData2 <- function(data, screening_interval, z) {
  resultsmatrix <- matrix(-1,n,maxD)
  
  for (t in testdays) {
    for (j in 1:n) {
      if (day_adm[j] <= t && day_dis[j] >= t) {
        if (col_t[j] <= t && col_t[j]!=0 && runif(1,0,1)<z) {
          resultsmatrix[j,t] <- 1
        } else {
          resultsmatrix[j,t] <- 0
        }
      }
    }
  }
}

screenData_old <- function(data, screening_interval, z) {
  N <- length(data$t_a)
  maxD <- max(data$t_d)
  results_matrix <- matrix(NA, nrow=N, ncol=maxD+1)
  for(i in 1:N) {
    test_days <- seq(from=data$t_a[i], to=data$t_d[i], by=screening_interval)
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


ReturnPossibleInfectorsR <- function(target, epi_data) {
  possible_infectors <- c()
  target_t_c <- epi_data$t_c[target]
  
  for(i in 1:length(epi_data$t_c)) {
    
  }
}


simulateGeneticData8 <- function(sequence_length, mutation_rate, data, screening_data, returnDistMat = TRUE) {
  positive_swabs <- which(screening_data==1,arr.ind=T)
  positive_swabs[,2] <- positive_swabs[,2] - 1 # adjust so that the column reflects time
  imported_case <- which(data$importations==1)
  ever_infected <- which(!is.na(data$t_c))
  sequences_on_infection <- matrix(NA,nrow=sequence_length,ncol=length(ever_infected))
  
  ## assign the imported case a random sequence
  sequences_on_infection[,which(ever_infected==imported_case)] <- sample(c("A","G","C","T"),sequence_length,replace=T)

  
  
  ## now assign the sequences observed
  
  sequence_matrix <- matrix(NA,nrow=sequence_length,ncol=nrow(positive_swabs))
  matrixFilled <- FALSE
  lastEventScreen <- FALSE
  primary_cases <- imported_case
  primary_times <- data$t_c[imported_case]
  while(!matrixFilled) {
    temp <- NULL
    for(i in primary_cases) {
      ## check what happens first, screening or infection (or both)
      screen_times <- positive_swabs[which(positive_swabs[,1]==i),2]
      secondary_cases <- which(data$source == i)
      if(length(secondary_cases)>0) temp <- c(temp,secondary_cases)
      secondary_times <- data$t_c[secondary_cases]
      combined_events <- rbind(c(rep(0,length(secondary_times)),rep(1,length(screen_times))),c(secondary_times,screen_times))
      if(ncol(combined_events)>1) {
        combined_events <- combined_events[,order(combined_events[2,])]
      }
      
      for(j in unique(combined_events[2,])) {
        screen <- FALSE
        event.locs <- which(combined_events[2,]==j)
        screen <- ifelse(sum(combined_events[1,event.locs]==1)==1,TRUE,FALSE)
        
        #check what the last event was
        if(min(event.locs)==1) {
          ### first event, copy sequence from infection
          time_diff <- j - data$t_c[i]
          current_sequence <- simulateNewSequence(sequences_on_infection[,which(ever_infected==i)], mutation_rate, time_diff)
        } else {
          ## not the first event, copy previous
          loc <- min(event.locs)-1
          time_diff <- j - combined_events[2,loc]
          if(combined_events[1,loc]==1) {
            #last event was a screen, copy from screen sequence matrix
            screen.time <- combined_events[2,loc]
            screen.loc <- which(positive_swabs[,2]==screen.time & positive_swabs[,1]==i)
            current_sequence <- simulateNewSequence(sequence_matrix[,screen.loc], mutation_rate, time_diff)
          } else {
            # last event was infection, copy from infection matrix
            inf_time <- combined_events[2,loc]
            person_infected <- min(which(data$source == i & data$t_c==inf_time))
            current_sequence <- simulateNewSequence(sequences_on_infection[,which(ever_infected==person_infected)], mutation_rate, time_diff)
          }
        }
        
        if(screen) {
          # if a screen has occured, add to sequence matrix
          #loc <- min(event.locs)
          #screen.time <- combined_events[2,loc]
          screen.loc <- which(positive_swabs[,2]==j & positive_swabs[,1]==i)
          sequence_matrix[,screen.loc] <- current_sequence
        }
        
        # check if anyone has been infected
        who_infected <- secondary_cases[secondary_times == j]
        for(k in who_infected) {
          current.loc <- which(ever_infected==k)
          sequences_on_infection[,current.loc] <- current_sequence
        }
        
      }
    }
    primary_cases <- temp
    if(length(primary_cases)==0) {
      matrixFilled <- TRUE
    }
  }
  
  
  
  
  
  if(returnDistMat) {
    distMat <- full_sequences_to_distance_matrix(sequence_matrix)
    colnames(distMat) <- positive_swabs[,1]
    rownames(distMat) <- positive_swabs[,1]
    return(distMat)
  } else {
    return(sequence_matrix)
  }
  
  

}


simulateGeneticDataImports <- function(sequence_length, mutation_rate, data, screening_data, returnDistMat = TRUE) {
  positive_swabs <- which(screening_data==1,arr.ind=T)
  positive_swabs[,2] <- positive_swabs[,2] - 1 # adjust so that the column reflects time
  imported_case <- which(data$importations==1)
  ever_infected <- which(!is.na(data$t_c))
  sequences_on_infection <- matrix(NA,nrow=sequence_length,ncol=length(ever_infected))
  
  ## assign the imported case a random sequence
  for(i in imported_case) {
    sequences_on_infection[,which(ever_infected==i)] <- sample(c("A","G","C","T"),sequence_length,replace=T)
  }
  
  
  
  ## now assign the sequences observed
  
  sequence_matrix <- matrix(NA,nrow=sequence_length,ncol=nrow(positive_swabs))
  matrixFilled <- FALSE
  lastEventScreen <- FALSE
  primary_cases <- imported_case
  primary_times <- data$t_c[imported_case]
  while(!matrixFilled) {
    temp <- NULL
    for(i in primary_cases) {
      ## check what happens first, screening or infection (or both)
      screen_times <- positive_swabs[which(positive_swabs[,1]==i),2]
      secondary_cases <- which(data$source == i)
      if(length(secondary_cases)>0) temp <- c(temp,secondary_cases)
      secondary_times <- data$t_c[secondary_cases]
      combined_events <- rbind(c(rep(0,length(secondary_times)),rep(1,length(screen_times))),c(secondary_times,screen_times))
      if(ncol(combined_events)>1) {
        combined_events <- combined_events[,order(combined_events[2,])]
      }
      
      for(j in unique(combined_events[2,])) {
        screen <- FALSE
        event.locs <- which(combined_events[2,]==j)
        screen <- ifelse(sum(combined_events[1,event.locs]==1)==1,TRUE,FALSE)
        
        #check what the last event was
        if(min(event.locs)==1) {
          ### first event, copy sequence from infection
          time_diff <- j - data$t_c[i]
          current_sequence <- simulateNewSequence(sequences_on_infection[,which(ever_infected==i)], mutation_rate, time_diff)
        } else {
          ## not the first event, copy previous
          loc <- min(event.locs)-1
          time_diff <- j - combined_events[2,loc]
          if(combined_events[1,loc]==1) {
            #last event was a screen, copy from screen sequence matrix
            screen.time <- combined_events[2,loc]
            screen.loc <- which(positive_swabs[,2]==screen.time & positive_swabs[,1]==i)
            current_sequence <- simulateNewSequence(sequence_matrix[,screen.loc], mutation_rate, time_diff)
          } else {
            # last event was infection, copy from infection matrix
            inf_time <- combined_events[2,loc]
            person_infected <- min(which(data$source == i & data$t_c==inf_time))
            current_sequence <- simulateNewSequence(sequences_on_infection[,which(ever_infected==person_infected)], mutation_rate, time_diff)
          }
        }
        
        if(screen) {
          # if a screen has occured, add to sequence matrix
          #loc <- min(event.locs)
          #screen.time <- combined_events[2,loc]
          screen.loc <- which(positive_swabs[,2]==j & positive_swabs[,1]==i)
          sequence_matrix[,screen.loc] <- current_sequence
        }
        
        # check if anyone has been infected
        who_infected <- secondary_cases[secondary_times == j]
        for(k in who_infected) {
          current.loc <- which(ever_infected==k)
          sequences_on_infection[,current.loc] <- current_sequence
        }
        
      }
    }
    primary_cases <- temp
    if(length(primary_cases)==0) {
      matrixFilled <- TRUE
    }
  }
  
  
  
  
  
  if(returnDistMat) {
    distMat <- full_sequences_to_distance_matrix(sequence_matrix)
    colnames(distMat) <- positive_swabs[,1]
    rownames(distMat) <- positive_swabs[,1]
    return(distMat)
  } else {
    return(sequence_matrix)
  }
  
  
  
}





## simulate a distance matrix rather than sequences
SimulateGeneticMatrix <- function(sequence_length, mutation_rate, data, results_matrix) {
  positive_swabs <- which(results_matrix==1,arr.ind=T)
  
}

## simulate genetic data of postive screen results and time of infection
simulateGeneticDataFull <- function(sequence_length, mutation_rate, data, screening_data, returnDistMat = TRUE) {
  positive_swabs <- which(screening_data==1,arr.ind=T)
  positive_swabs[,2] <- positive_swabs[,2] - 1 # adjust so that the column reflects time
  imported_case <- which(data$importations==1)
  ever_infected <- which(!is.na(data$t_c))
  sequences_on_infection <- matrix(NA,nrow=sequence_length,ncol=length(ever_infected))
  
  ## assign the imported case a random sequence
  sequences_on_infection[,which(ever_infected==imported_case)] <- sample(c("A","G","C","T"),sequence_length,replace=T)
  
  
  
  ## now assign the sequences observed
  
  sequence_matrix <- matrix(NA,nrow=sequence_length,ncol=nrow(positive_swabs))
  matrixFilled <- FALSE
  lastEventScreen <- FALSE
  primary_cases <- imported_case
  primary_times <- data$t_c[imported_case]
  while(!matrixFilled) {
    temp <- NULL
    for(i in primary_cases) {
      ## check what happens first, screening or infection (or both)
      screen_times <- positive_swabs[which(positive_swabs[,1]==i),2]
      secondary_cases <- which(data$source == i)
      if(length(secondary_cases)>0) temp <- c(temp,secondary_cases)
      secondary_times <- data$t_c[secondary_cases]
      combined_events <- rbind(c(rep(0,length(secondary_times)),rep(1,length(screen_times))),c(secondary_times,screen_times))
      if(ncol(combined_events)>1) {
        combined_events <- combined_events[,order(combined_events[2,])]
      }
      
      for(j in unique(combined_events[2,])) {
        screen <- FALSE
        event.locs <- which(combined_events[2,]==j)
        screen <- ifelse(sum(combined_events[1,event.locs]==1)==1,TRUE,FALSE)
        
        #check what the last event was
        if(min(event.locs)==1) {
          ### first event, copy sequence from infection
          time_diff <- j - data$t_c[i]
          current_sequence <- simulateNewSequence(sequences_on_infection[,which(ever_infected==i)], mutation_rate, time_diff)
        } else {
          ## not the first event, copy previous
          loc <- min(event.locs)-1
          time_diff <- j - combined_events[2,loc]
          if(combined_events[1,loc]==1) {
            #last event was a screen, copy from screen sequence matrix
            screen.time <- combined_events[2,loc]
            screen.loc <- which(positive_swabs[,2]==screen.time & positive_swabs[,1]==i)
            current_sequence <- simulateNewSequence(sequence_matrix[,screen.loc], mutation_rate, time_diff)
          } else {
            # last event was infection, copy from infection matrix
            inf_time <- combined_events[2,loc]
            person_infected <- min(which(data$source == i & data$t_c==inf_time))
            current_sequence <- simulateNewSequence(sequences_on_infection[,which(ever_infected==person_infected)], mutation_rate, time_diff)
          }
        }
        
        if(screen) {
          # if a screen has occured, add to sequence matrix
          #loc <- min(event.locs)
          #screen.time <- combined_events[2,loc]
          screen.loc <- which(positive_swabs[,2]==j & positive_swabs[,1]==i)
          sequence_matrix[,screen.loc] <- current_sequence
        }
        
        # check if anyone has been infected
        who_infected <- secondary_cases[secondary_times == j]
        for(k in who_infected) {
          current.loc <- which(ever_infected==k)
          sequences_on_infection[,current.loc] <- current_sequence
        }
        
      }
    }
    primary_cases <- temp
    if(length(primary_cases)==0) {
      matrixFilled <- TRUE
    }
  }
  
  noSequenceOnColonisation <- c()
  for(i in 1:length(ever_infected)) {
    #browser()
    #if(screening_data[ever_infected[i],data$t_c[ever_infected[i]]+1]!=1 || is.na(screening_data[ever_infected[i],data$t_c[ever_infected[i]]+1])) noSequenceOnColonisation <- c(noSequenceOnColonisation,i)
    if(is.na(screening_data[ever_infected[i],data$t_c[ever_infected[i]]+1])) noSequenceOnColonisation <- c(noSequenceOnColonisation,i)
    }
  #browser()
  sequence_matrix <- cbind(sequence_matrix, sequences_on_infection[,noSequenceOnColonisation])
  
  
  if(returnDistMat) {
    distMat <- full_sequences_to_distance_matrix(sequence_matrix)
    colnames(distMat) <- c(positive_swabs[,1],ever_infected[noSequenceOnColonisation])
    rownames(distMat) <- c(positive_swabs[,1],ever_infected[noSequenceOnColonisation])
    out <- list("distMat" = distMat, "sample_times" = c(positive_swabs[,2],data$t_c[ever_infected[noSequenceOnColonisation]]))
    return(out)
  } else {
    return(sequence_matrix)
  }
  
  
  
}


simulateGeneticDataFull2 <- function(sequence_length, mutation_rate, data, screening_data, returnDistMat = TRUE) {
  positive_swabs <- which(screening_data==1,arr.ind=T)
  positive_swabs[,2] <- positive_swabs[,2] - 1 # adjust so that the column reflects time
  imported_case <- which(data$importations==1)
  ever_infected <- which(!is.na(data$t_c))
  sequences_on_infection <- matrix(NA,nrow=sequence_length,ncol=length(ever_infected))
  
  ## assign the imported case a random sequence
  sequences_on_infection[,which(ever_infected==imported_case)] <- sample(c("A","G","C","T"),sequence_length,replace=T)
  
  
  
  ## now assign the sequences observed
  
  sequence_matrix <- matrix(NA,nrow=sequence_length,ncol=nrow(positive_swabs))
  matrixFilled <- FALSE
  lastEventScreen <- FALSE
  primary_cases <- imported_case
  primary_times <- data$t_c[imported_case]
  while(!matrixFilled) {
    temp <- NULL
    for(i in primary_cases) {
      ## check what happens first, screening or infection (or both)
      screen_times <- positive_swabs[which(positive_swabs[,1]==i),2]
      secondary_cases <- which(data$source == i)
      if(length(secondary_cases)>0) temp <- c(temp,secondary_cases)
      secondary_times <- data$t_c[secondary_cases]
      combined_events <- rbind(c(rep(0,length(secondary_times)),rep(1,length(screen_times))),c(secondary_times,screen_times))
      if(ncol(combined_events)>1) {
        combined_events <- combined_events[,order(combined_events[2,])]
      }
      
      for(j in unique(combined_events[2,])) {
        screen <- FALSE
        event.locs <- which(combined_events[2,]==j)
        screen <- ifelse(sum(combined_events[1,event.locs]==1)==1,TRUE,FALSE)
        
        #check what the last event was
        if(min(event.locs)==1) {
          ### first event, copy sequence from infection
          time_diff <- j - data$t_c[i]
          current_sequence <- simulateNewSequence(sequences_on_infection[,which(ever_infected==i)], mutation_rate, time_diff)
        } else {
          ## not the first event, copy previous
          loc <- min(event.locs)-1
          time_diff <- j - combined_events[2,loc]
          if(combined_events[1,loc]==1) {
            #last event was a screen, copy from screen sequence matrix
            screen.time <- combined_events[2,loc]
            screen.loc <- which(positive_swabs[,2]==screen.time & positive_swabs[,1]==i)
            current_sequence <- simulateNewSequence(sequence_matrix[,screen.loc], mutation_rate, time_diff)
          } else {
            # last event was infection, copy from infection matrix
            inf_time <- combined_events[2,loc]
            person_infected <- min(which(data$source == i & data$t_c==inf_time))
            current_sequence <- simulateNewSequence(sequences_on_infection[,which(ever_infected==person_infected)], mutation_rate, time_diff)
          }
        }
        
        if(screen) {
          # if a screen has occured, add to sequence matrix
          #loc <- min(event.locs)
          #screen.time <- combined_events[2,loc]
          screen.loc <- which(positive_swabs[,2]==j & positive_swabs[,1]==i)
          sequence_matrix[,screen.loc] <- current_sequence
        }
        
        # check if anyone has been infected
        who_infected <- secondary_cases[secondary_times == j]
        for(k in who_infected) {
          current.loc <- which(ever_infected==k)
          sequences_on_infection[,current.loc] <- current_sequence
        }
        
      }
    }
    primary_cases <- temp
    if(length(primary_cases)==0) {
      matrixFilled <- TRUE
    }
  }
  
  noSequenceOnColonisation <- c()
  for(i in 1:length(ever_infected)) {
    #browser()
    #if(screening_data[ever_infected[i],data$t_c[ever_infected[i]]+1]!=1 || is.na(screening_data[ever_infected[i],data$t_c[ever_infected[i]]+1])) noSequenceOnColonisation <- c(noSequenceOnColonisation,i)
    if(is.na(screening_data[ever_infected[i],data$t_c[ever_infected[i]]+1])) noSequenceOnColonisation <- c(noSequenceOnColonisation,i)
  }
  #browser()
  sequence_matrix <- cbind(sequence_matrix, sequences_on_infection[,noSequenceOnColonisation])
  
  
  if(returnDistMat) {
    distMat <- full_sequences_to_distance_matrix(sequence_matrix)
    colnames(distMat) <- c(positive_swabs[,1],ever_infected[noSequenceOnColonisation])
    rownames(distMat) <- c(positive_swabs[,1],ever_infected[noSequenceOnColonisation])
    out <- list("distMat" = distMat, "sample_times" = c(positive_swabs[,2],data$t_c[ever_infected[noSequenceOnColonisation]]))
    return(out)
  } else {
    return(sequence_matrix)
  }
  
  
  
}


## given an old sequence and time, generate a new sequence
simulateNewSequence <- function(sequence, mutation_rate, time) {
  mutation_probability <- JC_prob_mutation(rate = mutation_rate,t=time)
  sequence_length <- length(sequence)
  SNPs <- rbinom(1, sequence_length, mutation_probability)
  SNPs.loc <- sample(1:sequence_length, SNPs)
  for(j in SNPs.loc) {
    current_base <- sequence[j]
    sequence[j] <- sample(setdiff(c("A","C","G","T"),current_base),1)
  }
  return(sequence)
}

full_sequences_to_distance_matrix <- function(sequence_matrix) {
  number_sequences <- ncol(sequence_matrix)
  distance_matrix <- matrix(NA,ncol=number_sequences,nrow=number_sequences)
  for(i in 1:number_sequences) {
    for(j in 1:number_sequences) {
      distance_matrix[i,j] <- sum(sequence_matrix[,i]!=sequence_matrix[,j])
    }
  }
  return(distance_matrix)
}


marginal_genetic_prob <- function(N, x, rate, t) {
  prob <- (N-x)*log(1-JC_prob_mutation(rate, t)) + x*log(JC_prob_mutation(rate, t)/3)
  return(exp(prob))
}


returnNaiveSourceVector <- function(data) {
  source_vec <- data$true_source
  colonisedIDs <- which(data$t_c != -1)
  importID <- which(data$true_importations == 1)
  for(i in colonisedIDs) {
    #browser()
    if(!(i %in% importID)) {
      abletoinfect <- which(data$t_c < data$t_c[i] & data$t_d >= data$t_c[i] & data$t_c != -1)
      for(j in importID)
      {
        if(data$t_c[i]==data$t_c[j]) abletoinfect <- c(abletoinfect, j)
      }
      source_vec[i] <- sampleWithoutSurprises(abletoinfect)
    } else {
      source_vec[i] <- -1
    }

  }
  return(source_vec)
}


## Given perfect simulated data, return naive colonisation times, importation status and source
## Individuals have a source randomly selected at the colonisation time from when their first positive swab result is
## those without a source are considered an importation
ReturnNaiveSourceImportation <- function(epi_data, colonisation_times)
{
  colonisedIDs <- which(colonisation_times != -1)
  num_patients <- length(colonisation_times)
  source_vec <- rep(-2,num_patients)
  importations <- rep(0, num_patients)
  col_times <- rep(-1,num_patients)
  for(i in colonisedIDs)
  {
    abletoinfect <- which(colonisation_times < colonisation_times[i] & epi_data$t_d >= colonisation_times[i] & colonisation_times != -1)
    if(length(abletoinfect)==0)
    {
      source_vec[i] <- -1
      importations[i] <- 1
      col_times[i] <- epi_data$t_a[i]
    }
    else
    {
      source_vec[i] <- sampleWithoutSurprises(abletoinfect)
      col_times[i] <- colonisation_times[i]
    }
  }
  out <- list("source" = source_vec,
              "importations" = importations,
              "naive_coltimes" = col_times)
  return(out)
}

PlotParameter <- function(res, index, truth, cred_interval = FALSE, hist = FALSE) {
  if(hist) {
    hist(res[,index])
    abline(v=truth,col="red")
    if(cred_interval) abline(v=quantile(res[,index],c(0.05,0.95)),col="blue")
  } else {
    plot(res[,index],type="l")
    abline(h=truth,col="red")
    if(cred_interval) abline(h=quantile(res[,index],c(0.05,0.95)),col="blue")
  }
}

PlotTrace <- function(res,truth, cred_interval = FALSE) {
  dev.new()
  par(mfrow=c(4,1))
  PlotParameter(res,1,truth$true_z, cred_interval)
  PlotParameter(res,2,truth$true_p, cred_interval)
  PlotParameter(res,3,truth$true_beta, cred_interval)
  PlotParameter(res,4,truth$true_mutation_rate, cred_interval)
  par(mfrow=c(1,1))
}

PlotColoniationTimeSum <- function(colonisation_times, true_coltimes) {
  #colonisation_times <- read.table("coltimes.dat")
  col_times_sum <- apply(colonisation_times, 1, function(x) sum(x[x>0]))
  dev.new()
  plot(col_times_sum, type="l")
  abline(h=sum(true_coltimes[true_coltimes!=-1]),col=2)
}

SimulateEpiData <- function(true_z, true_p, true_beta, true_mutation_rate, N, num_patients, maxD, length_of_stay, screening_interval, constrain_final_size, assign_sequence_on_infection = FALSE, naive_coltimes = TRUE, naive_source = FALSE, impose_no_missing_infectives = TRUE) {
  
  dataGenerated <- FALSE
  while(!dataGenerated) {
    #data <- simulateWard(N=num_patients,D=maxD,p=true_p,length_of_stay=length_of_stay, z=0.5, a0 = 0, a1 = true_beta)
    data <- simulateWard2(N=num_patients,D=maxD,length_of_stay=length_of_stay, z=0.5, a0 = 0, a1 = true_beta)
    
    final_size <- sum(data$t_c>=0,na.rm=T)
    if(final_size > constrain_final_size[1] && final_size < constrain_final_size[2]) {
      #dataGenerated <- TRUE
      results_mat <- screenData(data,screening_interval,z=true_z)
      data$t_c[is.na(data$t_c)] <- -1
      if(naive_coltimes==TRUE) {
        col_times <- rep(-1,length(data$t_c))
        for(i in 1:length(data$t_c)) {
          positive_locs <- which(results_mat[i,]==1)
          if(length(positive_locs)>0)
          {
            col_times[i] <- min(positive_locs)-1
          }
        }
      } else {
        col_times <- data$t_c
      }
      if(impose_no_missing_infectives) {
        if(length(col_times[col_times>=0])==length(data$t_c[data$t_c>=0])) {
          dataGenerated <- TRUE
        }
      } else {
        dataGenerated <- TRUE
      }
      
      possible_epidemic <- IsEveryoneAbleToBeInfected(col_times, data$t_d, data$importations)
      if(possible_epidemic == FALSE) {
        dataGenerated <- FALSE
      }
    }
  }
  results_mat[is.na(results_mat)] <- -1
  #data$source[data$source != -1 & !is.na(data$source)] <- data$source[data$source != -1 & !is.na(data$source)] - 1
  data$source[is.na(data$source)] <- -2
  out <- list("t_a" = data$t_a,
              "t_c" = col_times,
              "t_d" = data$t_d,
              "results_matrix" = results_mat,
              "true_importations" = data$importations,
              "true_coltimes" = data$t_c,
              "true_source" = data$source)
  return(out)
}

SimulateDataMCMC <- function(true_z, true_p, true_beta, true_mutation_rate, N, num_patients, maxD, length_of_stay, screening_interval, constrain_final_size, assign_sequence_on_infection = FALSE, naive_coltimes = TRUE, naive_source = FALSE, impose_no_missing_infectives = TRUE) {
  
  dataGenerated <- FALSE
  while(!dataGenerated) {
    #data <- simulateWard(N=num_patients,D=maxD,p=true_p,length_of_stay=length_of_stay, z=0.5, a0 = 0, a1 = true_beta)
    data <- simulateWard2(N=num_patients,D=maxD,length_of_stay=length_of_stay, z=0.5, a0 = 0, a1 = true_beta)
    final_size <- sum(data$t_c>=0,na.rm=T)
    if(final_size > constrain_final_size[1] && final_size < constrain_final_size[2]) {
      #dataGenerated <- TRUE
      results_mat <- screenData(data,screening_interval,z=true_z)
      data$t_c[is.na(data$t_c)] <- -1
      if(naive_coltimes==TRUE) {
        col_times <- rep(-1,length(data$t_c))
        for(i in 1:length(data$t_c)) {
          positive_locs <- which(results_mat[i,]==1)
          if(length(positive_locs)>0)
          {
            col_times[i] <- min(positive_locs)-1
          }
        }
      } else {
        col_times <- data$t_c
      }
      if(impose_no_missing_infectives) {
        if(length(col_times[col_times>=0])==length(data$t_c[data$t_c>=0])) {
          dataGenerated <- TRUE
        }
      } else {
        dataGenerated <- TRUE
      }
      
      possible_epidemic <- IsEveryoneAbleToBeInfected(col_times, data$t_d, data$importations)
      if(possible_epidemic == FALSE) {
        dataGenerated <- FALSE
      }
    }
    
    
  }
  
  
  gen_mat <- simulateGeneticDataImports(N,true_mutation_rate,data,results_mat,T)
  sample_times <- as.numeric(which(results_mat==1,arr.ind=T)[,2]-1)
  genMatpID <- as.numeric(colnames(gen_mat)) - 1 # minus one for zero based indexing
  
  
  # adjust data for mcmc
  results_mat[is.na(results_mat)] <- -1
  data$source[data$source != -1 & !is.na(data$source)] <- data$source[data$source != -1 & !is.na(data$source)] - 1
  data$source[is.na(data$source)] <- -2
  
  if(naive_source) {
    source_vec <- returnNaiveSourceVector(data)
  } else {
    source_vec <- data$source
  }
  
  out <- list('t_a' = data$t_a,
              't_c' = col_times,
              't_d' = data$t_d,
              'importations' = data$importations,
              'source' = source_vec,
              'results_matrix' = results_mat,
              'genetic_id' = genMatpID,
              'sample_times' = sample_times,
              'genetic_matrix' = gen_mat,
              'true_coltimes' = data$t_c,
              'true_source' = data$source,
              'final_size' = final_size)
  return(out)
}



ReturnGeneticTree <- function(genetic_ids, source_vector, sample_times) {

  
  
  ## adjust for zero based indexing
  temp_source <- source_vector
  temp_source[temp_source >= 0] <- temp_source[temp_source >= 0] - 1
  gen_source <- ReturnGenSourceTree(genetic_ids - 1, temp_source, sample_times)
  
  ## put back to one for R
  gen_source[gen_source>=0] <- gen_source[gen_source>=0] + 1
  
  require(igraph)
  import_loc <- which(gen_source==-2)
  gen_source[import_loc] <- import_loc
  #browser()
  g <- make_graph(c(rbind(gen_source,(1:length(gen_source)))), n=length(genetic_ids))
  
  #plot(simplify(g), vertex.size = 8, vertex.label = genetic_ids, directed=T, edge.size=0.5, edge.arrow.size=0.3, edge.arrow.width=0.5)
  plot(simplify(g), vertex.size = 8, directed=T, edge.size=0.5, edge.arrow.size=0.3, edge.arrow.width=0.5)
  
}

ReturnGeneticTreeZeroBased <- function(genetic_ids, source_vector, sample_times) {
  
  
  

  gen_source <- ReturnGenSourceTree(genetic_ids, source_vector, sample_times) + 1


  require(igraph)
  import_loc <- which(gen_source==-1)
  gen_source[import_loc] <- import_loc
  #browser()
  g <- make_graph(c(rbind(gen_source,(1:(length(gen_source))))), n=length(genetic_ids))
  
  #plot(simplify(g), vertex.size = 8, vertex.label = genetic_ids, directed=T, edge.size=0.5, edge.arrow.size=0.3, edge.arrow.width=0.5)
  plot(simplify(g), vertex.size = 8, directed=T, edge.size=0.5, edge.arrow.size=0.3, edge.arrow.width=0.5, vertex.label = 0:(length(gen_source)-1))
  
}

CalculateDistanceBetweenNodes <- function(node1, node2, genetic_tree, genetic_matrix) {
  #if(node1==1 && node2 ==8) browser()
  dist <- 0
  if(node1 == node2) return(0)
  #browser()
  path1 <- ReturnPathToRoot(node1, genetic_tree)
  path2 <- ReturnPathToRoot(node2, genetic_tree)
  if(is.null(path1) || is.null(path2)) return(NULL)
  import_idx <- which(genetic_tree == -2)
  path1 <- path1[-length(path1)]
  path2 <- path2[-length(path2)]
  #path1[path1==-2] <- import_idx
  #path2[path2==-2] <- import_idx
  if(ReturnIsSamePathway(node1, node2, genetic_tree)) {
    # nodes are in the same pathway, just need to add between them
    if(length(path1) > length(path2)) {
      node2_loc <- which(path1 == node2)
      for(i in 1:(node2_loc-1)) {
        dist <- dist + genetic_matrix[path1[i],path1[i+1]]
      }
    } else {
      node1_loc <- which(path2 == node1)
      for(i in 1:(node1_loc-1)) {
        dist <- dist + genetic_matrix[path2[i],path2[i+1]]
      }
    }
  } else {
    ## nodes are in separate paths, find common root
    common_node <- ReturnCommonNode(path1,path2)
    if(is.null(common_node)) return(NULL)
    path1 <- path1[1:which(path1==common_node)]
    path2 <- path2[1:which(path2==common_node)]
    for(i in 1:(length(path1)-1)) {
      dist <- dist + genetic_matrix[path1[i],path1[i+1]]
    }
    for(j in 1:(length(path2)-1)) {
      dist <- dist + genetic_matrix[path2[j],path2[j+1]]
    }
  }
  return(dist)
}


CalculateDistanceBetweenNodes2 <- function(node1, node2, genetic_tree, genetic_matrix) {
  #if(node1==1 && node2 ==8) browser()
  dist <- 0
  if(node1 == node2) return(0)
  browser()
  path1 <- ReturnPathToRoot2(node1, genetic_tree)
  path2 <- ReturnPathToRoot2(node2, genetic_tree)
  if(is.null(path1) || is.null(path2)) return(NULL)
  import_idx <- which(genetic_tree == -2)
  path1 <- path1[-length(path1)]
  path2 <- path2[-length(path2)]
  #path1[path1==-2] <- import_idx
  #path2[path2==-2] <- import_idx
  if(ReturnIsSamePathway(node1, node2, genetic_tree)) {
    # nodes are in the same pathway, just need to add between them
    if(length(path1) > length(path2)) {
      node2_loc <- which(path1 == node2)
      for(i in 1:(node2_loc-1)) {
        dist <- dist + genetic_matrix[path1[i],path1[i+1]]
      }
    } else {
      node1_loc <- which(path2 == node1)
      for(i in 1:(node1_loc-1)) {
        dist <- dist + genetic_matrix[path2[i],path2[i+1]]
      }
    }
  } else {
    ## nodes are in separate paths, find common root
    common_node <- ReturnCommonNode(path1,path2)
    if(is.null(common_node)) return(NULL)
    path1 <- path1[1:which(path1==common_node)]
    path2 <- path2[1:which(path2==common_node)]
    for(i in 1:(length(path1)-1)) {
      dist <- dist + genetic_matrix[path1[i],path1[i+1]]
    }
    for(j in 1:(length(path2)-1)) {
      dist <- dist + genetic_matrix[path2[j],path2[j+1]]
    }
  }
  return(dist)
}



ReturnCommonNode <- function(path1,path2) {
  common_root_found <- FALSE
  i <- 0
  if(length(path1) < length(path2))
  {
    while(!common_root_found) {
      i <- i + 1
      if(i > length(path1))
      {
        print("inf loop")
        return(NULL)
      }
      if(path1[i] %in% path2) {
        common_root_found <- TRUE
        common_root <- path1[i]
      }
    }
  } else {
    while(!common_root_found) {
      i <- i + 1
      if(i > length(path2))
      {
        print("inf loop")
        return(NULL)
      }
      if(path2[i] %in% path1) {
        common_root_found <- TRUE
        common_root <- path2[i]
      }
    }
  }
  return(common_root)
}

ReturnIsSamePathway <- function(node1, node2, genetic_tree) {
  path1 <- ReturnPathToRoot(node1, genetic_tree)
  path2 <- ReturnPathToRoot(node2, genetic_tree)
  if(length(path1)>length(path2)) {
    if(node2 %in% path1) {
      return(TRUE)
    }
  } else {
    if(node1 %in% path2) return(TRUE)
  }
  return(FALSE)
}


ReturnPathToRoot <- function(node, genetic_tree) {
  source <- genetic_tree[node]
  path <- c(node, source)
  #browser()
  counter <- 1
  while(source != -2)
  {
    target <- source
    source <- genetic_tree[target]
    if(target == source)
    {
      print("infinite loop")
      return(NULL)
    }
    path <- c(path, source)
    counter <- counter + 1
    if(counter > length(genetic_tree)) return(NULL)
  }
  return(path)
}

ReturnPathToRoot2 <- function(node, genetic_tree) {
  source <- genetic_tree[node]
  path <- c(node, source)
  #browser()
  counter <- 1
  while(source != -1)
  {
    target <- source
    source <- genetic_tree[target]
    if(target == source)
    {
      print("infinite loop")
      return(NULL)
    }
    path <- c(path, source)
    counter <- counter + 1
    if(counter > length(genetic_tree)) return(NULL)
  }
  return(path)
}

SimulateGeneticData <- function(epi_data, sequence_length, true_mutation_rate) {
  positive_swabs <- which(epi_data$results_mat == 1, arr.ind = T)
  genetic_ids <- as.numeric(positive_swabs[,1])
  num_observed <- length(genetic_ids)
  sample_times <- as.numeric(positive_swabs[,2]) - 1
  colonised <- which(epi_data$t_c != -1)
  for(i in colonised) {
    min_swab_time <- min(sample_times[which(genetic_ids == i)])
    if(epi_data$true_coltimes[i] != min_swab_time) {
      genetic_ids <- c(genetic_ids, i)
      sample_times <- c(sample_times, epi_data$true_coltimes[i])
    }
  }
  
  ## adjust for zero based indexing
  #temp_source <- epi_data$true_source
  #temp_source[temp_source >= 0] <- temp_source[temp_source >= 0] - 1
  # gen_source <- ReturnGenSourceVector(genetic_ids - 1, temp_source, sample_times)
  #gen_test <- ReturnGenSourceTree(genetic_ids - 1, temp_source, sample_times)
  #gen_source
  #gen_test
  #ReturnGeneticTree(genetic_ids,epi_data$true_source, sample_times)
  
  #gen_source[gen_source>=0] <- gen_source[gen_source>=0] + 1 # adjust the zero based indexing for R
  
  ## modify the gen_source such that it describes the genet
  #gen_source
  
  temp_source <- epi_data$true_source
  temp_source[temp_source >= 0] <- temp_source[temp_source >= 0] - 1
  
  first_import <- min(which(epi_data$true_importations==1))-1
  
  genetic_tree <- ReturnGenSourceTree(genetic_ids - 1, temp_source, sample_times, first_import)
  genetic_matrix <- matrix(rep(NA,length(genetic_tree)^2),nrow=length(genetic_tree))
  genetic_tree[genetic_tree >= 0] <- genetic_tree[genetic_tree >= 0] + 1
  for(i in 1:length(genetic_tree)) {
    current_gen_source <- genetic_tree[i]
    if(current_gen_source != -2) {
      # imported strain
      time_diff <- abs(sample_times[current_gen_source]-sample_times[i])
      if(time_diff == 0) {
        genetic_matrix[i,current_gen_source] <- 0
        genetic_matrix[current_gen_source,i] <- 0
      } else {
        draw <- rbinom(1, sequence_length, JC_prob_mutation(true_mutation_rate, time_diff))
        genetic_matrix[i,current_gen_source] <- draw
        genetic_matrix[current_gen_source,i] <- draw
      }
    } else {
      genetic_matrix[i,i] <- 0
    }
  }
  
  for(i in 1:nrow(genetic_matrix)) 
  {
    for(j in 1:ncol(genetic_matrix)) 
    {
      if(is.na(genetic_matrix[i,j])) 
      {
        distance <- CalculateDistanceBetweenNodes(i,j,genetic_tree,genetic_matrix)
        if(is.null(distance)) return(NULL)
        genetic_matrix[i,j] <- distance
      }
    }
  }
  #browser()
  out <- list("genetic_ids" = genetic_ids[1:num_observed],
              "sample_times" = sample_times[1:num_observed],
              "genetic_matrix" = genetic_matrix[1:num_observed,1:num_observed])
  return(out)
}

SimulateDataMCMC2 <- function(true_z, true_p, true_beta, true_mutation_rate, N, num_patients, maxD, length_of_stay, screening_interval, 
                              constrain_final_size, assign_sequence_on_infection = FALSE, naive_coltimes = TRUE, naive_source = FALSE, 
                              impose_no_missing_infectives = TRUE) {
  
  epi_data <- SimulateEpiData(true_z, true_p, true_beta, true_mutation_rate, N, num_patients, maxD, 
                  length_of_stay, screening_interval, constrain_final_size, naive_coltimes = naive_coltimes, impose_no_missing_infectives = TRUE)
  gen_data <- SimulateGeneticData(epi_data, N, true_mutation_rate)
  
  if(naive_source) {
    source_vec <- returnNaiveSourceVector(epi_data)
  } else {
    source_vec <- epi_data$true_source
  }
  
  # adjust data for mcmc
  epi_data$true_source[epi_data$true_source>=0] <- epi_data$true_source[epi_data$true_source>=0] - 1
  source_vec[source_vec>=0] <- source_vec[source_vec>=0] - 1
  gen_data$genetic_ids <- gen_data$genetic_ids - 1
  
  

  
  final_size <- sum(epi_data$t_c>=0)
  
  out <- list('t_a' = epi_data$t_a,
              't_c' = epi_data$t_c,
              't_d' = epi_data$t_d,
              'importations' = epi_data$true_importations,
              'source' = source_vec,
              'results_matrix' = epi_data$results_matrix,
              'genetic_id' = gen_data$genetic_ids,
              'sample_times' = gen_data$sample_times,
              'genetic_matrix' = gen_data$genetic_matrix,
              'true_coltimes' = epi_data$true_coltimes,
              'true_source' = epi_data$true_source,
              'final_size' = final_size)
  return(out)
}
  
  
SimulateEpiData3 <- function(true_z, true_p, true_beta, true_mutation_rate, N, num_patients, maxD, length_of_stay, screening_interval, constrain_final_size, assign_sequence_on_infection = FALSE, naive_coltimes = TRUE, naive_source = FALSE, impose_no_missing_infectives = TRUE)
{
  dataGenerated <- FALSE
  while(!dataGenerated) {
    data <- simulateWard(N=num_patients,D=maxD,length_of_stay=length_of_stay, p=true_p, z=0.5, a0 = 0.0000001, a1 = true_beta)
    final_size <- sum(data$t_c>=0,na.rm=T)
    #browser()
    if(final_size > constrain_final_size[1] && final_size < constrain_final_size[2]) {
      #dataGenerated <- TRUE
      results_mat <- screenData(data,screening_interval,z=true_z)
      data$t_c[is.na(data$t_c)] <- -1
      if(naive_coltimes==TRUE) {
        col_times <- rep(-1,length(data$t_c))
        for(i in 1:length(data$t_c)) {
          positive_locs <- which(results_mat[i,]==1)
          if(length(positive_locs)>0)
          {
            col_times[i] <- min(positive_locs)-1
          }
        }
      } else {
        col_times <- data$t_c
      }
      if(impose_no_missing_infectives) {
        if(length(col_times[col_times>=0])==length(data$t_c[data$t_c>=0])) {
          dataGenerated <- TRUE
        }
      } else {
        dataGenerated <- TRUE
      }
      #browser()
      #possible_epidemic <- IsEveryoneAbleToBeInfected(col_times, data$t_d, data$importations)
      #if(possible_epidemic == FALSE) {
      #  dataGenerated <- FALSE
      #}
    }
  }
  results_mat[is.na(results_mat)] <- -1
  #data$source[data$source != -1 & !is.na(data$source)] <- data$source[data$source != -1 & !is.na(data$source)] - 1
  data$source[is.na(data$source)] <- -2
  out <- list("t_a" = data$t_a,
              "t_c" = col_times,
              "t_d" = data$t_d,
              "results_matrix" = results_mat,
              "true_importations" = data$importations,
              "true_coltimes" = data$t_c,
              "true_source" = data$source)
  return(out)
}

SimulateGeneticData3 <- function(epi_data, sequence_length, true_mutation_rate, p_seq) {
  positive_swabs <- which(epi_data$screening_matrix == 1, arr.ind = T)
  genetic_ids <- as.numeric(positive_swabs[,1])
  num_observed <- length(genetic_ids)
  sample_times <- as.numeric(positive_swabs[,2]) - 1
  colonised <- which(epi_data$t_c != -1)
  for(i in colonised) {
    min_swab_time <- min(sample_times[which(genetic_ids == i)])
    if(epi_data$true_coltimes[i] != min_swab_time) {
      genetic_ids <- c(genetic_ids, i)
      sample_times <- c(sample_times, epi_data$true_coltimes[i])
    }
  }

  unobserved_colonised <- setdiff(which(epi_data$true_coltimes!=-1),which(epi_data$t_c!=-1))
  for(i in unobserved_colonised)
  {
    genetic_ids <- c(genetic_ids, i)
    sample_times <-  c(sample_times, epi_data$true_coltimes[i])
  }
  
  ## adjust for zero based indexing
  #temp_source <- epi_data$true_source
  #temp_source[temp_source >= 0] <- temp_source[temp_source >= 0] - 1
  # gen_source <- ReturnGenSourceVector(genetic_ids - 1, temp_source, sample_times)
  #gen_test <- ReturnGenSourceTree(genetic_ids - 1, temp_source, sample_times)
  #gen_source
  #gen_test
  #ReturnGeneticTree(genetic_ids,epi_data$true_source, sample_times)
  
  #gen_source[gen_source>=0] <- gen_source[gen_source>=0] + 1 # adjust the zero based indexing for R
  
  ## modify the gen_source such that it describes the genet
  #gen_source
  #browser()
  genetic_ids <- c(genetic_ids, -2)
  sample_times <- c(sample_times, -1)
  first_import <- -3

  temp_source <- epi_data$true_source
  temp_source[temp_source >= 0] <- temp_source[temp_source >= 0] - 1
  genetic_tree <- ReturnGenSourceTree(genetic_ids - 1, temp_source, sample_times, first_import)
  genetic_matrix <- matrix(rep(NA,length(genetic_tree)^2),nrow=length(genetic_tree))
  genetic_tree[genetic_tree >= 0] <- genetic_tree[genetic_tree >= 0] + 1
  #epidemic_start_time <- min(epi_data$true_coltimes[which(epi_data$true_importations==1)])
  #browser()
  for(i in 1:length(genetic_tree)) {
    current_gen_source <- genetic_tree[i]
    person <- genetic_ids[i]
    if(current_gen_source == -2)
    {
      genetic_matrix[i,i] <- 0
    }
    else if(epi_data$true_importations[person]==1)
    {
        first_sequence <- FALSE
        person_sample_times <- sample_times[which(genetic_ids==person)]
        if(sample_times[i] == min(person_sample_times))
        {
          first_sequence <- TRUE
        }
        if(first_sequence)
        {
          draw <- rbinom(1, sequence_length, p_seq)
          genetic_matrix[i,current_gen_source] <- draw
          genetic_matrix[current_gen_source,i] <- draw
        }
        else
        {
          time_diff <- abs(sample_times[current_gen_source]-sample_times[i])
          if(time_diff == 0) {
            genetic_matrix[i,current_gen_source] <- 0
            genetic_matrix[current_gen_source,i] <- 0
          } else {
            draw <- rbinom(1, sequence_length, JC_prob_mutation(true_mutation_rate, time_diff))
            genetic_matrix[i,current_gen_source] <- draw
            genetic_matrix[current_gen_source,i] <- draw
          }
        }

    }
    else
    {
      time_diff <- abs(sample_times[current_gen_source]-sample_times[i])
      if(time_diff == 0) {
        genetic_matrix[i,current_gen_source] <- 0
        genetic_matrix[current_gen_source,i] <- 0
      } else {
        draw <- rbinom(1, sequence_length, JC_prob_mutation(true_mutation_rate, time_diff))
        genetic_matrix[i,current_gen_source] <- draw
        genetic_matrix[current_gen_source,i] <- draw
      }
    }
  }
  
  for(i in 1:nrow(genetic_matrix)) 
  {
    for(j in 1:ncol(genetic_matrix)) 
    {
      if(is.na(genetic_matrix[i,j])) 
      {
        distance <- CalculateDistanceBetweenNodes(i,j,genetic_tree,genetic_matrix)
        if(is.null(distance)) return(NULL)
        genetic_matrix[i,j] <- distance
      }
    }
  }
  #browser()
  out <- list("genetic_ids" = genetic_ids[1:num_observed],
              "sample_times" = sample_times[1:num_observed],
              "genetic_matrix" = genetic_matrix[1:num_observed,1:num_observed])
  return(out)
}

SimulateDataMCMC3 <- function(true_z, true_p, true_beta, true_mutation_rate, N, num_patients, maxD, length_of_stay, screening_interval, importation_p, 
                              constrain_final_size, assign_sequence_on_infection = FALSE, naive_coltimes = TRUE, naive_source = FALSE, 
                              impose_no_missing_infectives = TRUE) {
  
  dataGenerated <- FALSE
  while(!dataGenerated)
  {
    epi_data <- SimulateEpiData3(true_z, true_p, true_beta, true_mutation_rate, N, num_patients, maxD, 
                                 length_of_stay, screening_interval, constrain_final_size, naive_coltimes = naive_coltimes, impose_no_missing_infectives = impose_no_missing_infectives)
    gen_data <- SimulateGeneticData3(epi_data, N, true_mutation_rate, importation_p)
    if(!is.null(gen_data)) dataGenerated <- TRUE
  }

  
  if(naive_source) {
    #source_vec <- returnNaiveSourceVector(epi_data)
    naive_data <- ReturnNaiveSourceImportation(epi_data)
    source_vec <- naive_data$source
    importations <- naive_data$importations
    col_times <- naive_data$naive_coltimes
  } else {
    source_vec <- epi_data$true_source
    importations <- epi_data$true_importations
  }
  
  # adjust data for mcmc
  epi_data$true_source[epi_data$true_source>=0] <- epi_data$true_source[epi_data$true_source>=0] - 1
  source_vec[source_vec>=0] <- source_vec[source_vec>=0] - 1
  gen_data$genetic_ids <- gen_data$genetic_ids - 1
  
  
  
  
  final_size <- sum(epi_data$t_c>=0)
  
  out <- list('t_a' = epi_data$t_a,
              't_c' = col_times,
              't_d' = epi_data$t_d,
              'importations' = importations,
              'source' = source_vec,
              'results_matrix' = epi_data$results_matrix,
              'genetic_id' = gen_data$genetic_ids,
              'sample_times' = gen_data$sample_times,
              'genetic_matrix' = gen_data$genetic_matrix,
              'true_coltimes' = epi_data$true_coltimes,
              'true_source' = epi_data$true_source,
              'true_importations' = epi_data$true_importations,
              'final_size' = final_size)
  return(out)
}


SimulateGeneticData_WithinHostDiversity <- function(N, population_size, genetic_tree, sample_times, mutation_rate, 
                                                    stochastic_variant_p, dist_mat=TRUE, include_master_sequence=TRUE) {
  nucleotides <- c("A","G","C","T")
  master_sequence <- sample(nucleotides,size=N, replace=T)
  sequence_matrix <- matrix(NA,nrow=N,ncol=nrow(genetic_tree)*ncol(genetic_tree))
  #browser()
  ## sequences at time zero
  ## fill a matrix with the initial sequences, which are variances of the master sequence
  initial_sequence_matrix <- matrix(NA,nrow=N,ncol=population_size)
  initial_snps <- rbinom(population_size,N,stochastic_variant_p)
  SNP_locs <- sapply(initial_snps, function(x) sample(N,x))
  for(i in 1:length(SNP_locs)) {
    seq <- master_sequence
    seq[SNP_locs[[i]]] <- sapply(SNP_locs[[i]], function(x) sample(setdiff(nucleotides,seq[x]),1))
    initial_sequence_matrix[,i] <- seq
  }
  sequence_times <- c()
  #browser()
  
  for(i in 1:population_size) {
    
    initial_sequence <- initial_sequence_matrix[,i]
    
    
    for(j in 1:ncol(genetic_tree)) {
      sequence_source <- genetic_tree[i,j]
      if(sequence_source == -1) {
        ## the source sequence is the sequence at initialisation
        sequence_matrix[,(ncol(genetic_tree)*(i-1)+j)] <- initial_sequence
      } else {
        ## the source sequence is found in the genetic tree
        source_sequence <- genetic_tree[i,j]
        source_loc_in_matrix <- ((ncol(genetic_tree))*(i-1)+source_sequence)
        previous_seq <- sequence_matrix[,source_loc_in_matrix]
        new_seq <- simulateNewSequence(previous_seq, mutation_rate, sample_times[i,j]-sample_times[i,j-1])
        sequence_matrix[,(ncol(genetic_tree)*(i-1)+j)] <- new_seq
      }
      sequence_times <- c(sequence_times,sample_times[i,j])

      
    }
    
  }
  
  if(include_master_sequence) {
    sequence_matrix <- cbind(as.character(master_sequence),sequence_matrix)
    sequence_times <- c(-1,sequence_times)
  }
  sequence_matrix <- sequence_matrix[,order(sequence_times)]

  
  if(dist_mat) {
    dist_mat <- full_sequences_to_distance_matrix(sequence_matrix)
    out <- list("dist_mat" = dist_mat,
                "sample_times" = sequence_times[order(sequence_times)],
                "master_sequence" = master_sequence)
    return(out)
  } else {
    return(sequence_matrix)
  }
}


## simulates epi data (t_a, t_c, t_d, importations, source, swab matrix) and adjusts them suitable for MCMC
## including zero basedi ndexing for the source
SimulateEpiDataMCMC <- function(population_size, study_period, length_of_stay, screening_interval, transmission_rate, test_sensitivity, 
                    importation_probability) {
  if(importation_probability==0) {
    epi_data <- simulateWard2(N=population_size,D=study_period,length_of_stay=length_of_stay,z=test_sensitivity, a0 = 0.0000000000001, a1 = transmission_rate, mutation_rate = 3.3*10^(-6)/365)
  } else {
    epi_data <- simulateWard(N=population_size,D=study_period,length_of_stay=length_of_stay,p=importation_probability, z=test_sensitivity, a0 = 0.0000000000001, a1 = transmission_rate, mutation_rate = 3.3*10^(-6)/365)
    
  }
  transmission_network <- returnTransmissionNetwork(epi_data)
  #plot(transmission_network)
  screening_matrix <- screenData(epi_data, screening_interval, test_sensitivity)
  screening_matrix[is.na(screening_matrix)] <- -1
  
  
  ## Now generate naive estimates for coltimes based on sampling results
  ever_infected <- which(!is.na(epi_data$t_c))
  coltimes <- rep(NA, length(epi_data$t_c))
  for(i in 1:nrow(screening_matrix)) {
    positive_swabs <- which(screening_matrix[i,]==1)
    if(length(positive_swabs) > 0) {
      ## the person has a positive swab
      day_first_swab <- positive_swabs[1] - 1
      coltimes[i] <- day_first_swab
    } else {
      coltimes[i] <- -1
    }
  }
  
  ## now assign source and importation status
  source_vector <- rep(-2, length(epi_data$t_c))
  importations <- rep(0, length(epi_data$t_c))
  #browser()
  for(i in which(coltimes != -1)) {
    possible_infectors <- which(coltimes < coltimes[i] & epi_data$t_d >= coltimes[i] & coltimes != -1)
    if(length(possible_infectors)>0) {
      source_vector[i] <- sampleWithoutSurprises(possible_infectors) #- 1 ### -1 FOR ZERO BASED INDEXING
    } else {
      source_vector[i] <- -1
      importations[i] <- 1
    }
  }
  
  epi_data$source[is.na(epi_data$source)] <- -2
  #epi_data$source[epi_data$source >= 0] <- epi_data$source[epi_data$source >= 0] - 1

  epi_data$t_c[is.na(epi_data$t_c)] <- -1
  


  
  
  out <- list("t_a" = epi_data$t_a,
              "t_c" = coltimes,
              "t_d" = epi_data$t_d,
              "importations" = importations,
              "source" = source_vector,
              "screening_matrix" = screening_matrix,
              "transmission_network" = transmission_network,
              "true_coltimes" = epi_data$t_c,
              "true_importations" = epi_data$importations,
              "true_source" = epi_data$source)
  
  return(out)
}


SimulateGeneticDataWHD <- function(N, epi_data, genetic_ids, sample_times, number_variants, variant_numbers, source_vector, importation_p, mutation_rate) {
  #browser()

  ## generate the master sequence and then the master variants, and then append all sequences at the time of infection
  
  number_observed_sequences <- length(genetic_ids)
  
  
  ## add distances that are unobserved at the time of infection
  ever_infected <- which(epi_data$t_c != -1)
  for(i in ever_infected) {
    inf_time <- epi_data$t_c[i]
    for(j in 1:number_variants) {
      genetic_id_locs <- which(genetic_ids == i)
      inf_loc_in_gen_ids <- which(sample_times[genetic_id_locs] == inf_time & variant_numbers[genetic_id_locs] == j)
      if(length(inf_loc_in_gen_ids)==0) {
        genetic_ids <- c(genetic_ids, i)
        sample_times <- c(sample_times, inf_time)
        variant_numbers <- c(variant_numbers, j)
      }
    }
  }
  
  ## add an extra column for the master distance to be removed at the end
  full_distance_matrix <- matrix(NA, nrow=length(genetic_ids)+number_variants+1, ncol=length(genetic_ids)+number_variants+1)
  diag(full_distance_matrix) <- 0
  
  
  ## generate distances for each of the master variants
  for(i in 1:number_variants) {
    draw <- rbinom(1, N, importation_p)
    full_distance_matrix[1,i+1] <- draw
    full_distance_matrix[i+1,1] <- draw
  }
  
  gen_source <- sapply((0:(length(genetic_ids)-1)), function(x) ReturnLastSequenceID_variant(x, genetic_ids-1, source_vector, sample_times, variant_numbers))
  
  
  ## add the imported distances which are drawn to be different from the master sequence with probability importation_p
  imported_idx <- which(gen_source==-1)
  #browser()
  for(i in imported_idx) {
    time_diff <- sample_times[i]
    draw <- rbinom(1,N, JC_prob_mutation(mutation_rate, time_diff))
    current_variant <- variant_numbers[i]
    full_distance_matrix[1+current_variant,i+1+number_variants] <- draw
    full_distance_matrix[i+1+number_variants,1+current_variant] <- draw
  }
  
  browser()
  #new_distance_matrix <- matrix(NA, nrow=length(genetic_ids)+number_variants+1, ncol=length(genetic_ids)+number_variants+1)
  #new_distance_matrix[(1:(number_observed_sequences+number_variants+1)),(1:(number_observed_sequences+number_variants+1))] <- full_distance_matrix
  #full_distance_matrix <- new_distance_matrix
  
  ## now go through the directly links and simulate them from the JC model
  for(i in 1:length(gen_source)) {
    target <- i
    genetic_source <- gen_source[i]+1
    if(genetic_source > 0)
    {
      time_diff <- sample_times[target]-sample_times[genetic_source]
      draw <- rbinom(1,N,JC_prob_mutation(mutation_rate, time_diff))
      full_distance_matrix[target+1+number_variants,genetic_source+1+number_variants] <- draw
      full_distance_matrix[genetic_source+1+number_variants,target+1+number_variants] <- draw
    }
  }
  full_gen_tree <- gen_source
  for(i in 1:length(gen_source)) {
    target <- i
    genetic_source <- gen_source[i]+1
    if(genetic_source == -1)
    {
      ## find the source and copy
      patient_id <- genetic_ids[i] ## note one-based indexing
      patient_source <- source_vector[patient_id] + 1
      #browser()
      loc_to_copy <- which(genetic_ids == patient_source & sample_times == sample_times[i] & variant_numbers == variant_numbers[i])
      full_distance_matrix[i+1+number_variants,] <- full_distance_matrix[loc_to_copy+1+number_variants,]
      full_distance_matrix[,i+1+number_variants] <- full_distance_matrix[,loc_to_copy+1+number_variants]
      full_gen_tree[i] <- loc_to_copy+1+number_variants
    }
  }
  #browser()
  ## fill the rest of the matrix deterministically
  full_gen_tree <- c(-5,-4,-4,full_gen_tree)
  full_gen_tree[full_gen_tree >=0] <- full_gen_tree[full_gen_tree >=0] + 2 + number_variants
  full_gen_tree[which(full_gen_tree == -1)] <- variant_numbers[which(full_gen_tree == -1)-number_variants-1] + 1
  full_gen_tree[full_gen_tree==-4] <- 1
  

  full_gen_tree[1] <- -2
  #return(full_gen_tree)
  
  
  for(i in 1:length(full_gen_tree)) {
    for(j in 1:length(full_gen_tree)) {
      #browser()
      full_distance_matrix[i,j] <- CalculateDistanceBetweenNodes(i,j,full_gen_tree,full_distance_matrix)
    }
  }
  return(full_distance_matrix)
}

SimulateGeneticDataWHD2 <- function(epi_data, N, genetic_ids, sample_times, variant_numbers, importation_p, mutation_rate) {
  ## suppose there are n observed samples (number_observed_samples), k variants (num_variants)
  num_observed_sampled <- length(genetic_ids)
  num_variants <- length(unique(variant_numbers))
  
  ## for each person infected, add their variant and a genetic id at the time of infection, and note that genetic ids are one-based for the minute
  ever_infected <- which(epi_data$t_c != -1)
  for(i in ever_infected) {
    time_of_infection <- epi_data$t_c[i]
    for(j in 1:num_variants) {
      ## check if their variant is stored in the genetic ids and sample times
      gen_loc <- which(genetic_ids == i & sample_times == time_of_infection & variant_numbers == j)
      if(length(gen_loc)==0) {
        ## there is no id at the time of infection, therefore append each of the vectors and move on
        genetic_ids <- c(genetic_ids, i)
        sample_times <- c(sample_times, time_of_infection)
        variant_numbers <- c(variant_numbers, j)
      }
    }
  }
  

  ## create a matrix to store the one master sequence, k variants and the genetic ids
  distance_matrix <- matrix(NA, nrow=1+num_variants+length(genetic_ids), ncol=1+num_variants+length(genetic_ids))
  diag(distance_matrix) <- 0
  
  ## fill in the links from the master sequence to each of the master variants
  for(i in 1:num_variants) {
    draw <- rbinom(1, N, importation_p)
    distance_matrix[1,1+i] <- draw
    distance_matrix[1+i,1] <- draw
  }
  
  ## fill in the rest of the direct links from the JC model
  gen_source <- sapply((0:(length(genetic_ids)-1)), function(x) ReturnLastSequenceID_variant(x, genetic_ids-1, epi_data$source, sample_times, variant_numbers))
  #browser()
  for(i in 1:length(gen_source)) {
    current_gen_source <- gen_source[i] + 1 ## note that gen source is zero based, so we add one for R
    if(current_gen_source > 0) {
      ## the genetic source is from another id, rather than a master variant
      location_idx <- i + num_variants + 1
      source_idx <- current_gen_source + num_variants + 1
      time_diff <- abs(sample_times[i]-sample_times[current_gen_source])
      mutation_probability <- JC_prob_mutation(mutation_rate,time_diff)
      draw <- rbinom(1,N,mutation_probability)
      distance_matrix[source_idx,location_idx] <- draw
      distance_matrix[location_idx,source_idx] <- draw
    } else if(current_gen_source == 0) {
      ## look at the master variant
      current_variant <- variant_numbers[i]
      source_idx <- current_variant + 1
      location_idx <- i + num_variants + 1
      time_diff <- sample_times[i]
      mutation_probability <- JC_prob_mutation(mutation_rate,time_diff)
      draw <- rbinom(1,N,mutation_probability)
      distance_matrix[source_idx,location_idx] <- draw
      distance_matrix[location_idx,source_idx] <- draw
    }
  }
  #browser()
  
  ## now copy any duplicated sequences (sequences that were taken at the same time as infection)
  gen_tree <- gen_source
  for(i in which(gen_tree == -2)) {
    patient_id <- genetic_ids[i]
    source <- epi_data$source[patient_id] + 1 ## note source is zero based indexing, so add one
    loc_to_copy <- which(genetic_ids==source & sample_times == sample_times[i] & variant_numbers == variant_numbers[i])
    distance_matrix[i+num_variants+1,which(!is.na(distance_matrix[loc_to_copy+num_variants+1,]))] <- distance_matrix[loc_to_copy+num_variants+1,which(!is.na(distance_matrix[loc_to_copy+num_variants+1,]))]
    distance_matrix[which(!is.na(distance_matrix[loc_to_copy+num_variants+1,])),i+num_variants+1] <- distance_matrix[which(!is.na(distance_matrix[loc_to_copy+num_variants+1,])),loc_to_copy+num_variants+1]
    distance_matrix[i+num_variants+1,loc_to_copy+num_variants+1] <- 0
    distance_matrix[loc_to_copy+num_variants+1,i+num_variants+1] <- 0
    gen_tree[i] <- loc_to_copy-1
  }

  
  ## adjust gen tree so that is corresponds to an actual path, with -2 at the root node (master sequence)
  ## first shift everyone such that they point to the correct location
  ## then adjust the variants
  ## then append the master and the variants
  gen_tree[gen_tree>=0] <- gen_tree[gen_tree>=0] + 1+ 1 + num_variants
  gen_tree[gen_tree==-1] <- variant_numbers[which(gen_tree==-1)] + 1
  gen_tree <- c(-2,1,1,gen_tree)
  
  ## now fill in the rest of the direct links
  
  for(i in 1:length(gen_tree)) {
    for(j in 1:length(gen_tree)) {
      #if(i == 1 && j == 7) browser()
      distance_matrix[i,j] <- CalculateDistanceBetweenNodes(i,j,gen_tree,distance_matrix)
    }
  }
  return(distance_matrix[(1+num_variants+1):(num_observed_sampled+1+num_variants),(1+num_variants+1):(num_observed_sampled+1+num_variants)])
  
}




PlotGeneticTree <- function(gen_source) {
  library(igraph)
  full_gen_tree <- gen_source
  full_gen_tree[1] <- 1
  tree <- make_graph(c(rbind(full_gen_tree,(1:length(full_gen_tree)))), n=length(full_gen_tree));
  plot(tree, edge.arrow.size=0.2, vertex.size = 1, edge.arrow.width = 0.5)
}

ReturnGeneticTree <- function(gen_source) {
  library(igraph)
  full_gen_tree <- c()
  for(i in 1:length(gen_source)) {
    if(gen_source[i] != -1) {
      full_gen_tree <- c(full_gen_tree, gen_source[i]+1,i)
    }
  }
  g <- make_graph(full_gen_tree)
  return(g)
}

## finds the pathway between a common source
ReturnTransmissionPathwayToSource <- function(originator, target, source_vector) {
  path <- c(target)
  source <- target
  #browser()
  while(source != originator) {
    source <- source_vector[target]
    path <- c(source, path)
    target <- source
  }
  return(path)
}

ReturnNodeAndTimeOfDivergence <- function(originator, target1, target2, source_vector, t_c) {
  path1 <- ReturnTransmissionPathwayToSource(originator, target1, source_vector+1)
  path2 <- ReturnTransmissionPathwayToSource(originator, target2, source_vector+1)
  for(i in 2:min(c(length(path1),length(path2)))) {
    if(path1[i] != path2[i]) {
      node_div_idx <- i
    }
  }
  
  col_times <- t_c[c(path1[node_div_idx],path2[node_div_idx])]
  earlier_node <- which.min(col_times)
  if(earlier_node == 1) {
    node_div <- path1[node_div_idx]
  } else {
    node_div <- path2[node_div_idx]
  }
  time_div <- t_c[node_div]
  out <- list("node" = node_div,
              "time" = time_div)
}

#ReturnTransmissionPathwayToSource(1475, 1781, data$source_vector+1)
#ReturnTransmissionPathwayToSource(1475, 1451, data$source_vector+1)


CalculateGeneticLikelihood_WHD2 <- function(N, mutation_rate, genetic_ids, sample_times, t_c, variant_numbers, source_vector, genetic_matrix) {
  gen_source_vector <- ReturnGenSourceVector_WHD(genetic_ids, source_vector, sample_times, variant_numbers, t_c)
  loglik <- 0
  for(i in 1:length(gen_source_vector)) {
    genetic_source <- gen_source_vector[i] + 1
    if(genetic_source > 0) {
      distance <- genetic_matrix[i,genetic_source]
      time_diff <- sample_times[i] - sample_times[genetic_source]
      if(time_diff > 0) {
        probability <- JC_prob_mutation(mutation_rate, time_diff)
        loglik <- loglik + (N-distance)*log(1-probability) + distance*log(probability/3);
      }

    }
  }
  return(loglik)
}

CalculateGeneticLikelihood_gen_source_WHD2 <- function(N, mutation_rate, genetic_ids, sample_times, t_c, variant_numbers, source_vector, gen_source_vector, genetic_matrix) {
  loglik <- 0
  for(i in 1:length(gen_source_vector)) {
    genetic_source <- gen_source_vector[i] + 1
    if(genetic_source > 0) {
      distance <- genetic_matrix[i,genetic_source]
      time_diff <- sample_times[i] - sample_times[genetic_source]
      if(time_diff > 0) {
        probability <- JC_prob_mutation(mutation_rate, time_diff)
        loglik <- loglik + (N-distance)*log(1-probability) + distance*log(probability/3);
      }
      
    }
  }
  return(loglik)
}



ReturnKmeansClustering <- function(matrix, num_clusters) {
  library(magrittr)
  library(dplyr)
  

  # Cmpute MDS
  mds <- matrix %>%
    dist() %>%          
    cmdscale() %>%
    as_tibble()
  colnames(mds) <- c("Dim.1", "Dim.2")
  # Plot MDS
  
  # K-means clustering
  j <- 2
  clust <- kmeans(mds, num_clusters)$cluster %>%
    as.factor()
  mds2 <- mds %>% mutate(groups = clust)
  # Plot and color by groups
  #dev.new()

  
  
  variant_numbers <- mds2$groups
  return(variant_numbers)
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



ReturnGenSourceChain <- function(gen_source, original_length, target) {
  #browser()
  path <- c(target)
  target <- gen_source[target]
  path <- c(path, target)
  while(target > original_length) {
    target <- gen_source[target]
    path <- c(path, target)
  }
  return(path)
}


## Note genetic IDs should be one based indexing, and source vector is zero based indexing
## (I DONT KNOW WHY)
ReturnNodesToImpute <- function(t_c, source_vector, genetic_ids, sample_times, variant_numbers) {
  ever_infected <- which(data$t_c != -1)
  
  nodes_found <- FALSE
  while(!nodes_found) {
    #browser()
    additional_genetic_ids <- genetic_ids
    additional_sample_times <- sample_times
    additional_variant_numbers <- variant_numbers
    current_variant_to_add <- 1
    ## Append IDs and sample times at the time of infection for the test vector
    for(i in ever_infected) {
      gen_loc <- which(genetic_ids == i & sample_times == t_c[i] & variant_numbers == current_variant_to_add)
      if(length(gen_loc)==0) {
        ## There is no observed genetic sequence at the time of infection, therefore append them to the additional vector
        additional_genetic_ids <- c(additional_genetic_ids, i)
        additional_sample_times <- c(additional_sample_times, t_c[i])
        additional_variant_numbers <- c(additional_variant_numbers, current_variant_to_add)
      }
    }
    
    ## Return the gen source vector and add one for one based indexing
    additional_gen_source_vector <- ReturnGenSourceVector_WHD(additional_genetic_ids-1, source_vector, additional_sample_times, t_c, additional_variant_numbers) + 1
    
    ## Isolate which notes have sources that are nodes to later be imputed
    nodes_to_keep <- which(additional_gen_source_vector[1:length(genetic_ids)] > length(genetic_ids))
    #browser()
    if(length(nodes_to_keep)==0) {
      nodes_found <- TRUE
      break
    }
    
    ## Return the genetic source chain for each of these individuals, and keep the immediate previous link
    updated_genetic_ids <- genetic_ids
    updated_sample_times <- sample_times
    updated_variant_numbers <- variant_numbers
    for(i in 1:length(nodes_to_keep)) {
      current_chain <- ReturnGenSourceChain(additional_gen_source_vector, length(genetic_ids), nodes_to_keep[i])
      imp_node_to_keep <- current_chain[length(current_chain)-1]
      loc_in_updated_list <- which(updated_genetic_ids == additional_genetic_ids[imp_node_to_keep] & updated_sample_times == additional_sample_times[imp_node_to_keep] & updated_variant_numbers == additional_variant_numbers[imp_node_to_keep])
      if(length(loc_in_updated_list)==0) {
        updated_genetic_ids <- c(updated_genetic_ids, additional_genetic_ids[imp_node_to_keep])
        updated_sample_times <- c(updated_sample_times, additional_sample_times[imp_node_to_keep])
        updated_variant_numbers <- c(updated_variant_numbers, additional_variant_numbers[imp_node_to_keep])
      }
    }
    
    genetic_ids <- updated_genetic_ids
    sample_times <- updated_sample_times
    variant_numbers <- updated_variant_numbers
  }
  
  out <- list("genetic_ids" = genetic_ids,
              "sample_times" = sample_times,
              "variant_numbers" = variant_numbers)
  return(out)

}


SimulateGeneticData_WHD <- function(epi_data, variant_numbers, p_master, p_variant, N, mutation_rate) {
  
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
  gen_source <- ReturnGenSourceVector_WHD(temp_data_list)
  #length(which((gen_source == -1)))
  
  
  
  ## now simulate 
  
  full_distance_matrix <- matrix(NA, nrow=length(genetic_ids)+num_variants+1,ncol=length(genetic_ids)+num_variants+1)
  
  #browser()
  for(i in 1:num_variants) {
    draw <- rbinom(1,N,p_master)
    full_distance_matrix[1,i+1] <- draw
    full_distance_matrix[i+1,1] <- draw
  }
  
  for(i in 1:length(gen_source)) {
    current_source <- gen_source[i] + 1
    full_distance_matrix[i,i] <- 0
    if(current_source == 0) {
      ## first introduction of the pathogen in this chain, compare to master distance
      draw <- rbinom(1,N,p_variant)
      current_variant <- variant_numbers[i]
      full_distance_matrix[1+current_variant,i+1+num_variants] <- draw
      full_distance_matrix[i+1+num_variants,1+current_variant] <- draw
    } else {
      time_diff <- sample_times[i] - sample_times[current_source]
      if(time_diff != 0) {
        prob <- JC_prob_mutation(mutation_rate, time_diff)
        draw <- rbinom(1, N, prob)
      } else {
        draw <- 0
      }
      full_distance_matrix[1+num_variants+current_source,i+1+num_variants] <- draw
      full_distance_matrix[i+1+num_variants,1+num_variants+current_source] <- draw
    }
  }
  
  
  
  #browser()
  genetic_tree <- c(-2, rep(1,num_variants), gen_source+1+1+num_variants)
  
  
  import_variant_numbers <- variant_numbers[which(gen_source == -1)]
  genetic_tree[which(gen_source == -1)+1+num_variants] <- import_variant_numbers + 1
  
  
  for(i in 1:nrow(full_distance_matrix)) 
  {
    for(j in 1:ncol(full_distance_matrix)) 
    {
      if(is.na(full_distance_matrix[i,j])) 
      {
        distance <- CalculateDistanceBetweenNodes(i,j,genetic_tree,full_distance_matrix)
        if(is.null(distance)) return(NULL)
        full_distance_matrix[i,j] <- distance
      }
    }
  }
  #browser()
  
  output <- 1
  
  if(output == 1) {
    genetic_matrix <- full_distance_matrix[(2+num_variants):(1+num_variants+observed_length),(2+num_variants):(1+num_variants+observed_length)]
    out <- list("genetic_ids" = genetic_ids[1:observed_length],
                "sample_times" = sample_times[1:observed_length],
                "genetic_matrix" = genetic_matrix)
  } else if(output == 2) {
    genetic_matrix <- full_distance_matrix[(2+num_variants):nrow(full_distance_matrix),(2+num_variants):nrow(full_distance_matrix)]
    
    out <- list("genetic_ids" = genetic_ids,
                "sample_times" = sample_times,
                "genetic_matrix" = genetic_matrix)
  }
  
  
  
  return(out)
}

SimulateGeneticData2_WHD <- function(epi_data, variant_numbers, p_master, p_variant, N, mutation_rate) {
  
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
  gen_source <- ReturnGenSourceVector_WHD(temp_data_list)
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
  
  imports <- which(gen_source == -1)
  for(i in 1:(length(imports)-1)) {
    for(j in (i+1):length(imports)) {
      draw <- rpois(1,40)
      full_distance_matrix[imports[j],imports[i]] <- draw
      full_distance_matrix[imports[i],imports[j]] <- draw
    }
  }
  
  browser()
  #genetic_tree <- c(-2, gen_source+1)
  genetic_tree <- gen_source[gen_source>=0] + 1
  #gen_source[gen_source>=0]
  #import_variant_numbers <- variant_numbers[which(gen_source == -1)]
  #genetic_tree[which(gen_source == -1)+1+num_variants] <- import_variant_numbers + 1
  
  
  
  for(i in 1:nrow(full_distance_matrix)) 
  {
    for(j in 1:ncol(full_distance_matrix)) 
    {
      if(is.na(full_distance_matrix[i,j])) 
      {
        distance <- CalculateDistanceBetweenNodes2(i,j,genetic_tree,full_distance_matrix)
        if(is.null(distance)) return(NULL)
        full_distance_matrix[i,j] <- distance
      }
    }
  }
  #browser()
  
  output <- 1
  
  if(output == 1) {
    genetic_matrix <- full_distance_matrix[(2+num_variants):(1+num_variants+observed_length),(2+num_variants):(1+num_variants+observed_length)]
    out <- list("genetic_ids" = genetic_ids[1:observed_length],
                "sample_times" = sample_times[1:observed_length],
                "genetic_matrix" = genetic_matrix)
  } else if(output == 2) {
    genetic_matrix <- full_distance_matrix[(2+num_variants):nrow(full_distance_matrix),(2+num_variants):nrow(full_distance_matrix)]
    
    out <- list("genetic_ids" = genetic_ids,
                "sample_times" = sample_times,
                "genetic_matrix" = genetic_matrix)
  }
  
  
  
  return(out)
}

colpop_R <- function(epi_data, day) {
  #browser()
  acq <- which(epi_data$t_c < day & epi_data$t_c != -1 & epi_data$t_d >= day & epi_data$source != -1)
  imp <- which(epi_data$t_c <= day & epi_data$t_c != -1 & epi_data$t_d >= day & epi_data$source == -1)
  return(length(acq)+length(imp))
}

ReturnPossibleInfectors_R <- function(epi_data, target) {
  day <- epi_data$t_c[target]
  #browser()
  acq <- which(epi_data$t_c < day & epi_data$t_c != -1 & epi_data$t_d >= day & epi_data$source != -1)
  imp <- which(epi_data$t_c <= day & epi_data$t_c != -1 & epi_data$t_d >= day & epi_data$source == -1)
  out <- setdiff(c(acq,imp),target)
  if(length(c(acq,imp))>0) {
    return(out)
  } else {
    return(c())
  }
}

ReturnNadd0 <- function(epi_data, Va)
{
  nadd0 <- c()
  for(i in 1:length(Va)) {
    person <- Va[i]
    loc <- which((person-1)==epi_data$source)
    if(length(loc)==0) {
      nadd0 <- c(nadd0, person)
    }
  }
  return(nadd0)
}

CalculateFalseNegatives_R <- function(epi_data) {
  screening_matrix <- epi_data$screening_matrix
  t_d <- epi_data$t_d
  t_c <- epi_data$t_c
  ever_infected <- which(t_c != -1)
  false_negatives <- 0
  for(i in ever_infected) {
    person <- i
    screen_results <- screening_matrix[i,(t_c[person]:t_d[person])+1]
    false_negatives <- false_negatives + sum(screen_results == 0)
    
  }
  return(false_negatives)
}


MCMC_R <- function(MCMC_options, t_a, t_c, t_d, source, screening_matrix)
{
  epi_data <- list("t_a" = t_a,
                   "t_c" = t_c,
                   "t_d" = t_d,
                   "source" = source,
                   "screening_matrix" = screening_matrix)
  epi_data$source[epi_data$source > 0] <- epi_data$source[epi_data$source > 0] - 1
  
  Vq <- which(epi_data$t_c != -1)
  Vs <- which(epi_data$t_c == -1)
  Va <- c()
  Va_vec <- numeric()
  Va_vec[1] <- 0
  
  res <- matrix(NA,nrow=MCMC_options$iterations, ncol=3)
  res[1,] <- MCMC_options$initial_chain_state
  beta_cur <- res[1,3]
  TP <- sum(epi_data$screening_matrix == 1)
  nacc <- 0
  
  coltime_sum <- numeric()
  coltime_sum[1] <- sum(epi_data$t_c[epi_data$t_c != -1])
  loglik_cur <- CalculateLogLikelihood(epi_data, MCMC_options$initial_chain_state)
  counter <- 1
  
  
  for(i in 2:MCMC_options$iterations) {
    
    #if(i %% 5000 == 0) browser()
    
    
    z_cur <- rbeta(1,MCMC_options$prior_parameters[1] + TP, MCMC_options$prior_parameters[2] + CalculateFalseNegatives(epi_data))
    
    p_cur <- rbeta(1,MCMC_options$prior_parameters[3] + sum(epi_data$source == -1), MCMC_options$prior_parameters[4] + length(epi_data$source) - sum(epi_data$source == -1))
    
    loglik_cur <- CalculateLogLikelihood(epi_data, c(z_cur,p_cur,beta_cur))
    
    if(MCMC_options$debug_flags[3]==0) {
      beta_can <- rnorm(1, mean=beta_cur, sd=MCMC_options$proposal_variance[1])
      if(beta_can>0) {
        logpi_cur <- loglik_cur + dexp(beta_cur,rate=MCMC_options$prior_parameters[5])
        loglik_can <- CalculateLogLikelihood(epi_data, c(z_cur,p_cur,beta_can))
        logpi_can <- loglik_can + dexp(beta_can,rate=MCMC_options$prior_parameters[5])
        if(log(runif(1)) < logpi_can - logpi_cur) {
          loglik_cur <- loglik_can
          beta_cur <- beta_can
          nacc <- nacc + 1
        }
      }
    }

    
    
    ### augmented data updates
    num_updates <- 5
    w <- 0.3
    for(j in 1:num_updates) {
      #move <- sample(3,1)
      move <- 1
      epi_data_can <- epi_data
      if(move==1)
      {
        # move a colonisation time
        #browser()
        target <- sample(c(Vq,Va),1)
        last_day <- epi_data$t_d[target]
        pos_days <- which(epi_data$screening_matrix[target,] == 1) - 1
        if(length(pos_days)>0) {
          last_day <- min(last_day, pos_days)
        }
        offspring <- which(epi_data$source == (target-1))
        offspring_coltime <- epi_data$t_c[offspring]
        last_day <- min(last_day, offspring_coltime)
        
        U <- runif(1)
        #browser()
        if(U < w)
        {
          # propose an importation
          epi_data_can$t_c[target] <- epi_data$t_a[target]
          epi_data_can$source[target] <- -1
          if(epi_data$source[target] == -1) {
            # Imp -> Imp
            log_prop_ratio <- 0
          } else {
            log_prop_ratio <- log(1-w) - log(w) - log(last_day-epi_data$t_a[target]+1) - log(colpop_R(epi_data,epi_data$t_c[target]))
            if(last_day-t_a[target]+1 <= 0) browser()
          }
        } else {
          # propose an acquisition
          epi_data_can$t_c[target] <- sampleWithoutSurprises(epi_data$t_a[target]:last_day)
          possible_infectors <- ReturnPossibleInfectors_R(epi_data_can, target)
          if(length(possible_infectors)>0) {
            epi_data_can$source[target] <- sampleWithoutSurprises(possible_infectors) - 1
            if(epi_data$source[target] == -1) {
              # Imp -> Acq
              log_prop_ratio <- log(w) + log(last_day-t_a[target]+1) + log(colpop_R(epi_data_can,epi_data_can$t_c[target])) - log(1-w)
              if(last_day-t_a[target]+1 <= 0) browser()
            } else {
              log_prop_ratio <- log(colpop_R(epi_data_can,epi_data_can$t_c[target])) - log(colpop_R(epi_data,epi_data$t_c[target]))
            }
            
            loglik_can <- CalculateLogLikelihood(epi_data_can, c(z_cur,p_cur,beta_cur))
            #browser()
            if(is.nan(loglik_can - loglik_cur + log_prop_ratio)) browser()
            if(log(runif(1)) < loglik_can - loglik_cur + log_prop_ratio) {
              epi_data <- epi_data_can
              loglik_cur <- loglik_can
            }
            
          }
        }

      } else if(move==2) {
        ## add a colonisation
        #browser()
        susceptible_pop <- setdiff(Vs,Va)
        if(length(susceptible_pop)==0) break
        target <- sampleWithoutSurprises(susceptible_pop)
        nadd0 <- ReturnNadd0(epi_data, Va)
        if(runif(1)<w) {
          ## propose they are an importation
          epi_data_can$t_c[target] <- epi_data$t_a[target]
          epi_data_can$source[target] <- -1
          log_prop_ratio <- log(length(susceptible_pop)) - log(w) - log(1+length(nadd0))
        } else {
          ## propose they are an acquisition
          epi_data_can$t_c[target] <- sampleWithoutSurprises(epi_data$t_a[target]:epi_data$t_d[target])
          possible_infectors <- ReturnPossibleInfectors_R(epi_data_can, target)
          if(length(possible_infectors)==0) break 
          epi_data_can$source[target] <- sampleWithoutSurprises(possible_infectors) - 1
          log_prop_ratio <- log(colpop_R(epi_data_can, epi_data_can$t_c[target])) + log(length(susceptible_pop)) 
                            + log(epi_data$t_d[target] - epi_data$t_a[target] + 1) - log(1-w) - log(length(nadd0)+1)
          
        }
        loglik_can <- CalculateLogLikelihood(epi_data_can, c(z_cur,p_cur,beta_cur))
        #browser()
        if(is.nan(loglik_can - loglik_cur + log_prop_ratio)) browser()
        if(log(runif(1)) < loglik_can - loglik_cur + log_prop_ratio) {
          epi_data <- epi_data_can
          loglik_cur <- loglik_can
          Va <- c(Va,target)
        }

        
      } else if(move == 3) {
        ## remove a colonisation
        #browser()
        if(length(Va)==0) break
        #target <- sampleWithoutSurprises(Va)
        #if(sum((target-1)==epi_data$source)!=0) break
        nadd0 <- ReturnNadd0(epi_data, Va)
        if(length(nadd0)==0) break
        #browser()
        target <- sampleWithoutSurprises(nadd0)
        
        epi_data_can$t_c[target] <- -1
        epi_data_can$source[target] <- -2
        #browser()
        if(is.na(epi_data$source[target])) browser()
        if(epi_data$source[target] == -1) {
          # we are removing an importation
          log_prop_ratio <- log(length(nadd0)) + log(w) - log(length(Vs)-length(Va)+1)
        } else {
          log_prop_ratio <- log(length(nadd0)) + log(1-w) - log(epi_data$t_d[target]-epi_data$t_a[target]+1) - 
                            log(length(Vs)-length(Va)+1) - log(colpop_R(epi_data, epi_data$t_c[target]))
        }
        
        loglik_can <- CalculateLogLikelihood(epi_data_can, c(z_cur,p_cur,beta_cur))
        #browser()
        if(is.nan(loglik_can - loglik_cur + log_prop_ratio)) browser()
        if(log(runif(1)) < loglik_can - loglik_cur + log_prop_ratio) {
          epi_data <- epi_data_can
          loglik_cur <- loglik_can
          Va_loc <- which(Va == target)
          Va <- Va[-Va_loc]
        }
      }
      
      
      
      
      Va_vec[counter] <- length(Va)
      coltime_sum[counter] <- sum(epi_data$t_c[epi_data$t_c != -1])
      counter <- counter + 1
    }
    
    res[i,] <- c(z_cur,p_cur,beta_cur)
    
    out <- list('res' = res,
                'coltime_sum' = coltime_sum,
                'Va' = Va_vec)
    
    if(i%%10) {
      print(i)
    }
  }
  return(out)
}

MCMC_NS_R <- function(MCMC_options, t_a, t_c, t_d, source, screening_matrix)
{
  epi_data <- list("t_a" = t_a,
                   "t_c" = t_c,
                   "t_d" = t_d,
                   "source" = source,
                   "screening_matrix" = screening_matrix)
  epi_data$source[epi_data$source > 0] <- epi_data$source[epi_data$source > 0] - 1
  
  Vq <- which(epi_data$t_c != -1)
  Vs <- which(epi_data$t_c == -1)
  Va <- c()
  Va_vec <- numeric()
  Va_vec[1] <- 0
  
  res <- matrix(NA,nrow=MCMC_options$iterations, ncol=3)
  res[1,] <- MCMC_options$initial_chain_state
  beta_cur <- res[1,3]
  TP <- sum(epi_data$screening_matrix == 1)
  nacc <- 0
  
  coltime_sum <- numeric()
  coltime_sum[1] <- sum(epi_data$t_c[epi_data$t_c != -1])
  loglik_cur <- CalculateLogLikelihood_NS(epi_data, MCMC_options$initial_chain_state)
  counter <- 1
  
  
  for(i in 2:MCMC_options$iterations) {
    
    #if(i %% 5000 == 0) browser()
    
    
    #z_cur <- rbeta(1,MCMC_options$prior_parameters[1] + TP, MCMC_options$prior_parameters[2] + CalculateFalseNegatives(epi_data))
    false_negatives_R <- CalculateFalseNegatives_R(epi_data)
    false_negatives_C <- CalculateFalseNegatives(epi_data)
    if(false_negatives_R != false_negatives_C) browser()
    z_cur <- rbeta(1,MCMC_options$prior_parameters[1] + TP, MCMC_options$prior_parameters[2] + CalculateFalseNegatives_R(epi_data))
    
    
    p_cur <- rbeta(1,MCMC_options$prior_parameters[3] + sum(epi_data$source == -1), MCMC_options$prior_parameters[4] + length(epi_data$source) - sum(epi_data$source == -1))
    
    loglik_cur <- CalculateLogLikelihood_NS(epi_data, c(z_cur,p_cur,beta_cur))
    
    if(MCMC_options$debug_flags[3]==0) {
      beta_can <- rnorm(1, mean=beta_cur, sd=MCMC_options$proposal_variance[1])
      if(beta_can>0) {
        logpi_cur <- loglik_cur + dexp(beta_cur,rate=MCMC_options$prior_parameters[5])
        loglik_can <- CalculateLogLikelihood_NS(epi_data, c(z_cur,p_cur,beta_can))
        logpi_can <- loglik_can + dexp(beta_can,rate=MCMC_options$prior_parameters[5])
        if(log(runif(1)) < logpi_can - logpi_cur) {
          loglik_cur <- loglik_can
          beta_cur <- beta_can
          nacc <- nacc + 1
        }
      }
    }
    
    
    
    ### augmented data updates
    num_updates <- 5
    w <- 0.3
    for(j in 1:num_updates) {
      #move <- sample(3,1)
      move <- 1
      epi_data_can <- epi_data
      if(move==1)
      {
        # move a colonisation time
        #browser()
        target <- sample(c(Vq,Va),1)
        last_day <- epi_data$t_d[target]
        pos_days <- which(epi_data$screening_matrix[target,] == 1) - 1
        if(length(pos_days)>0) {
          last_day <- min(last_day, pos_days)
        }
        #offspring <- which(epi_data$source == (target-1))
        #offspring_coltime <- epi_data$t_c[offspring]
        #last_day <- min(last_day, offspring_coltime)
        
        U <- runif(1)
        #browser()
        if(U < w)
        {
          # propose an importation
          epi_data_can$t_c[target] <- epi_data$t_a[target]
          epi_data_can$source[target] <- -1
          if(epi_data$source[target] == -1) {
            # Imp -> Imp
            log_prop_ratio <- 0
          } else {
            log_prop_ratio <- log(1-w) - log(w) - log(last_day-epi_data$t_a[target]+1)
            if(last_day-t_a[target]+1 <= 0) browser()
          }
        } else {
          # propose an acquisition
          epi_data_can$t_c[target] <- sampleWithoutSurprises(epi_data$t_a[target]:last_day)
          possible_infectors <- ReturnPossibleInfectors_R(epi_data_can, target)
          if(length(possible_infectors)>0) {
            epi_data_can$source[target] <- sampleWithoutSurprises(possible_infectors) - 1
            if(epi_data$source[target] == -1) {
              # Imp -> Acq
              log_prop_ratio <- log(w) + log(last_day-t_a[target]+1) - log(1-w)
              if(last_day-t_a[target]+1 <= 0) browser()
            } else {
              #log_prop_ratio <- log(colpop_R(epi_data_can,epi_data_can$t_c[target])) - log(colpop_R(epi_data,epi_data$t_c[target]))
              log_prop_ratio <- 0
            }
            
            loglik_can <- CalculateLogLikelihood_NS(epi_data_can, c(z_cur,p_cur,beta_cur))
            #browser()
            if(is.nan(loglik_can - loglik_cur + log_prop_ratio)) browser()
            if(log(runif(1)) < loglik_can - loglik_cur + log_prop_ratio) {
              epi_data <- epi_data_can
              loglik_cur <- loglik_can
            }
            
          }
        }
        
      } else if(move==2) {
        ## add a colonisation
        #browser()
        susceptible_pop <- setdiff(Vs,Va)
        if(length(susceptible_pop)==0) break
        target <- sampleWithoutSurprises(susceptible_pop)
        nadd0 <- ReturnNadd0(epi_data, Va)
        if(runif(1)<w) {
          ## propose they are an importation
          epi_data_can$t_c[target] <- epi_data$t_a[target]
          epi_data_can$source[target] <- -1
          log_prop_ratio <- log(length(susceptible_pop)) - log(w) - log(1+length(nadd0))
        } else {
          ## propose they are an acquisition
          epi_data_can$t_c[target] <- sampleWithoutSurprises(epi_data$t_a[target]:epi_data$t_d[target])
          possible_infectors <- ReturnPossibleInfectors_R(epi_data_can, target)
          if(length(possible_infectors)==0) break 
          epi_data_can$source[target] <- sampleWithoutSurprises(possible_infectors) - 1
          log_prop_ratio <- log(colpop_R(epi_data_can, epi_data_can$t_c[target])) + log(length(susceptible_pop)) 
          + log(epi_data$t_d[target] - epi_data$t_a[target] + 1) - log(1-w) - log(length(nadd0)+1)
          
        }
        loglik_can <- CalculateLogLikelihood(epi_data_can, c(z_cur,p_cur,beta_cur))
        #browser()
        if(is.nan(loglik_can - loglik_cur + log_prop_ratio)) browser()
        if(log(runif(1)) < loglik_can - loglik_cur + log_prop_ratio) {
          epi_data <- epi_data_can
          loglik_cur <- loglik_can
          Va <- c(Va,target)
        }
        
        
      } else if(move == 3) {
        ## remove a colonisation
        #browser()
        if(length(Va)==0) break
        #target <- sampleWithoutSurprises(Va)
        #if(sum((target-1)==epi_data$source)!=0) break
        nadd0 <- ReturnNadd0(epi_data, Va)
        if(length(nadd0)==0) break
        #browser()
        target <- sampleWithoutSurprises(nadd0)
        
        epi_data_can$t_c[target] <- -1
        epi_data_can$source[target] <- -2
        #browser()
        if(is.na(epi_data$source[target])) browser()
        if(epi_data$source[target] == -1) {
          # we are removing an importation
          log_prop_ratio <- log(length(nadd0)) + log(w) - log(length(Vs)-length(Va)+1)
        } else {
          log_prop_ratio <- log(length(nadd0)) + log(1-w) - log(epi_data$t_d[target]-epi_data$t_a[target]+1) - 
            log(length(Vs)-length(Va)+1) - log(colpop_R(epi_data, epi_data$t_c[target]))
        }
        
        loglik_can <- CalculateLogLikelihood(epi_data_can, c(z_cur,p_cur,beta_cur))
        #browser()
        if(is.nan(loglik_can - loglik_cur + log_prop_ratio)) browser()
        if(log(runif(1)) < loglik_can - loglik_cur + log_prop_ratio) {
          epi_data <- epi_data_can
          loglik_cur <- loglik_can
          Va_loc <- which(Va == target)
          Va <- Va[-Va_loc]
        }
      }
      
      
      
      
      Va_vec[counter] <- length(Va)
      coltime_sum[counter] <- sum(epi_data$t_c[epi_data$t_c != -1])
      counter <- counter + 1
    }
    
    res[i,] <- c(z_cur,p_cur,beta_cur)
    
    out <- list('res' = res,
                'coltime_sum' = coltime_sum,
                'Va' = Va_vec)
    
    if(i%%10) {
      print(i)
    }
  }
  return(out)
}


ReturnCorrectEdgeProportion <- function(res, p_L) {
  correct_edges <- 0
  ever_infected <- which(res$simulation_info$epi_data$true_coltimes != -1)
  for(j in ever_infected) {
    posterior_probability_table <- table(res$source[,j])/length(res$source[,j])
    truth <- res$simulation_info$epi_data$true_source[j]
    if(truth > 0) truth <- truth - 1
    truth_loc <- which(names(posterior_probability_table) == truth)
    if(length(truth_loc) > 0) {
      if(posterior_probability_table[truth_loc] > p_L) {
        correct_edges <- correct_edges + 1
      }
    }
  }
  correct_edges
  correct_edges/length(ever_infected)
  #out <- list ("correct_edges" = correct_edges/length(ever_infected),
  #             "false_edges" = 1-correct_edges/length(ever_infected))
  #return(out)
}

SimulateGeneticData3_WHD <- function(epi_data, variant_numbers, average_import_distance, N, mutation_rate) {
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
  gen_source <- ReturnGenSourceVector_WHD(temp_data_list)
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
  
  imports <- which(gen_source == -1)
  if(length(imports)>1) {
    for(i in 1:(length(imports)-1)) {
      for(j in (i+1):length(imports)) {
        draw <- rpois(1,average_import_distance)
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
        distance <- CalculateDistanceBetweenNodes3(i,j,full_distance_matrix, gen_source)
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

SimulateGeneticData4_WHD <- function(epi_data, variant_numbers, average_import_distance, average_variant_distance, N, mutation_rate) {
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
  gen_source <- ReturnGenSourceVector_WHD(temp_data_list)
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
        distance <- CalculateDistanceBetweenNodes3(i,j,full_distance_matrix, gen_source)
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


ReturnPathToRoot3 <- function(node, gen_source) {
  #browser()
  path <- c(node)
  source <- gen_source[node] + 1
  while(source != 0) {
    path <- c(source, path)
    source <- gen_source[source] + 1
  }
  return(path)
}


CalculateDistanceBetweenNodes3 <- function(i,j,distance_matrix,gen_source) {
  if(i==j) return(0)
  path1 <- ReturnPathToRoot3(i,gen_source)
  path2 <- ReturnPathToRoot3(j,gen_source)
  
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
    common_node <- ReturnCommonNode3(path1,path2)
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

ReturnCommonNode3 <- function(path1, path2) {
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



