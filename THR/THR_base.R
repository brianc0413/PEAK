# comparisons with existing methods

# 4 arms - one above, one slightly above, one slightly below, one below
# HDoC(1), LUCB-G(2), APT-G(3)
#n_sims <- 1

epsilon <- 0.5
alpha <- 0.05
k <- 4
n_sims <- 100
data <- data.frame(method = NA,
                   stopping_time =NA,
                   arm_1_good = NA,
                   arm_2_good = NA,
                   arm_3_good = NA,
                   arm_4_good = NA)

# k = number of arms
# t = total time index
# S_list = data collected so far on the arm, indexed by arm
# alpha = confidence parameter
# epsilon = threshold
# Method: 1 = HDoC, 2 = LUCB-G, APT-G
# tells which arm to pull
sample_function <- function(k, t, S_list, alpha, epsilon, method, active_set){
  upper_bounds <- rep(NA, length(S_list))
  for (i in 1:length(S_list)){ # i.e. for each arm, find the upper bounds
    S <- S_list[[i]]
    mu_hat <- mean(S)
    upper_bounds[i] <- mu_hat + radius_pull(k, t, S, alpha, epsilon, method) 
  }
  return(which.max(upper_bounds[active_set]))
}

# S is a single vector of recorded observations
radius_pull <- function(k, t, S, alpha, epsilon, method){
  radius <- NA
  if (method == 1){
    radius <- sqrt(log(t)/2/length(S))
  }
  if (method == 2){
    radius <- sqrt( log(4*k*(length(S)^2)/alpha) /2/length(S))
  }
  if (method == 3){
    radius <- sqrt(length(S)) * abs(epsilon - mean(S))
  }
  return(radius)
}

# arm elimination function
arm_eliminator <- function(k, t, S_list, alpha, epsilon){
  good_or_bad <- rep(0, length(S_list))
  for (i in 1:length(S_list)){
    S <- S_list[[i]]
    radius <- sqrt( log(4*k*(length(S)^2)/alpha) /2/length(S))
    mu_hat <- mean(S)
    
    # check if a good arm
    if (mu_hat - radius > epsilon){
      good_or_bad[i] <- 1
    }
    
    # bad arms get a 0
    if (mu_hat + radius < epsilon){
      good_or_bad[i] <- -1
    }
    
  }
   return(good_or_bad)
}

for (sim in 1:n_sims){
  print(sim)
  set.seed(sim)
  
  # large variance setting
  #S_1 <- rbinom(1,1,0.29)
  #S_2 <- rbinom(1,1,0.43)
  #S_3 <- rbinom(1,1,0.57)
  #S_4 <- rbinom(1,1,0.71)
    
  # for small variance setting
  S_1 <- rbeta(1,1,0.71/.29)
  S_2 <- rbeta(1,1,0.57/.43)
  S_3 <- rbeta(1,1,0.43/.57)
  S_4 <- rbeta(1,1,0.29/.79)
    
  S_list <- list(S_1,S_2,S_3,S_4)
  t <- 4
    
  active_set <- 1:4
  arm_above_below <- rep(NA, length(S_list))
    
    
  while( length(active_set) != 0 ){
    # determine which arm to sample
    which_arm <- active_set[sample_function(k, t, S_list, alpha, epsilon, method=1, active_set)]
    #which_arm <- sample(1:4, size=1, prob = rep(1/4, 4))
    t <- t+1
    #S_list[[which_arm]] <- c(S_list[[which_arm]], rbinom(1,1, 0.29 + 0.14*(which_arm - 1)))
    S_list[[which_arm]] <- c(S_list[[which_arm]], 
                             rbeta(1,1, (0.71 - 0.14*(which_arm - 1))/(0.29 + 0.14*(which_arm-1))) )
      
    # determine if an arm can be considered above/below a threshold
    good_or_bad <- arm_eliminator(k, t, S_list, alpha, epsilon)
    good_or_bad_relevant <- good_or_bad[which_arm]
      
    if (good_or_bad_relevant == 1 && which_arm %in% active_set){
      arm_above_below[which_arm] <- 1
      active_set <- active_set[active_set != which_arm]
      data <- rbind(data, c(method, t, arm_above_below))
    }
    if (good_or_bad_relevant == -1 && which_arm %in% active_set){
      arm_above_below[which_arm] <- 0
      active_set <- active_set[active_set != which_arm]
      data <- rbind(data, c(1, t, arm_above_below))
    }
  }
}

# check for correctness
data[rowSums(is.na(data))==0, 3:6]

# get stopping times and standard deviations
taus <- rep(0,4)
std_devs <- rep(0,4)

for (i in 0:3){
  taus[4-i] <- mean(data$stopping_time[rowSums(is.na(data)) == i])
  std_devs[4-i] <- sd(data$stopping_time[rowSums(is.na(data)) == i])
}

print(taus)
print(std_devs)

