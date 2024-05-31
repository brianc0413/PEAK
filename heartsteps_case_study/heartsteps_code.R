dat <- read.csv("data_use.csv")
source("heartsteps_helper.R")


head(dat)


users <- unique(dat$user.index)


taus <- data.frame(Sub_Bernoulii = rep(NA, length(users)),
                   Hedged = rep(NA, length(users)),
                   PEAK = rep(NA, length(users)))
best_arms <- data.frame(Sub_Bernoulii = rep(NA, length(users)),
                   Hedged = rep(NA, length(users)),
                   PEAK = rep(NA, length(users)))

alpha <- 0.3


for (user in 6:37){
  
  rel_data <- dat[dat$user.index == user, ]
  A <- rel_data$send.sedentary
  X <-  rel_data$steps
  X <- rel_data$steps / max(X) # normalize to between 0 and 1
  stopped <- rep(0,3)
  partitions <- rep(0,2)
  
  # initialize pseudosamples
  S_list <- list(arm_1 = mean(X[A==1]),  
                 arm_0 = mean(X[A==0]))
  
  mus <- list(1/2, 1/2)
  var_hats <- list(1/4, 1/4)
  
  t <- 2
  
  while (sum(stopped) < 3 && t < length(X) + 1){
    
    h <- A[t]+1
    
    # update variances and means
    var_hats[[h]] <- c(var_hats[[h]], 
                       (1/4 + sum( (S_list[[h]]-mean(S_list[[h]]))^2 ) )/(length(S_list[[h]])+1) )
    mus[[h]] <- c(mus[[h]], (1/2+sum(S_list[[h]])) / (length(S_list[[h]]) +1) )
    
    
    S_list[[h]] <- c(S_list[[h]], X[t])
    t <- t+1
    
    
    # bernoulli confidence sequence
    
    if (stopped[1] == 0){
      terms <- term_condition(S_list, var_hats, 2, alpha)
      
      if (terms[1] == 1){
        stopped[1] <- 1
        taus[user, 1] <- t
        best_arms[user, 1] <- terms[2]
      }
    }
    
  
    # Hedged stopping criteria
    if (stopped[2] == 0){
      terms <- term_condition(S_list, var_hats, 3, alpha)
      if (terms[1] == 1){
        stopped[2] <- 1
        taus[user, 2] <- t
        best_arms[user, 2] <- terms[2]
      }
      
    }
    
    
    
    
    # PEAK stopping criteria
    
    if (stopped[3] == 0){
      min_values <- rep(NA, length(S_list))
      for (i in 1:length(S_list)){
        S <- S_list[[i]]
        mu_hat <- mus[[i]]
        f1 <- function(x){sum((2*x-mu_hat-S)/(0.26 + (S-x)*(mu_hat-x)))}
        val <- findZeros(f1(x)~x, xlim=c(0,1))
        min_values[i] <- val[1,1]
      }
      
      # determine which arm is currently best
      ranks <- rank(min_values, ties.method = "random") 
      best_arm <- which(ranks == 2)
      
      # minimum e-value globally corresponds to partition with best arm.
      test_values <- rep(0, length(S_list))
      test_values[best_arm] <- f(S_list, mus, min_values)
      
      # for the second best arm, we only need to project second best, best
      second_best <- which(ranks == 1)
      f2 <- function(x){
        coeffs(x, S_list, mus, best_arm) * sum((2*x-mus[[best_arm]]-S_list[[best_arm]])/(0.26 + (S_list[[best_arm]]-x)*(mus[[best_arm]]-x))) + 
          coeffs(x, S_list, mus, second_best) * sum((2*x-mus[[second_best]]-S_list[[second_best]])/(0.26 + (S_list[[second_best]]-x)*(mus[[second_best]]-x)))
      }
      val <- findZeros(f2(x) ~ x, xlim=c(0,1))
      temp <- min_values
      temp[best_arm] <- val[1,1]
      temp[second_best] <- val[1,1]
      test_values[second_best] <- f(S_list, mus, temp)
      
      reject <- (test_values >= 1/alpha) # 0 if eliminated 
      partitions <- pmax(reject, partitions) # if eliminated before, doesnt matter
      
      if (sum(reject) == 1){
        stopped[3] <- 1
        taus[user, 3] <- t
        best_arms[user, 3] <- which.min(test_values)
      }
      
      print(test_values)
    }
    print(taus[user,])
    print(best_arms[user, ])
    print(t)
  }
  
  
}

diffs <- c()
for (user in users){
  diffs <- c(diffs, (mean(dat$steps[dat$user.index==user & dat$send.sedentary==1])
             - mean(dat$steps[dat$user.index==user & dat$send.sedentary==0])) /
               max(dat$steps[dat$user.index == user]) )
}


# how many samples could we have saved?

total_samples <- 0
total_samples_vec <- c()
for (i in which(!is.na(taus$PEAK))){
  total_samples <- total_samples +sum(dat$user.index == i)
  total_samples_vec <- c(total_samples_vec, sum(dat$user.index == i))
}

# table
stop_times_peak <- taus[which(!is.na(taus$PEAK)), 3]
tab_dat <- data.frame(Unit_Index = which(!is.na(taus$PEAK)),
                      Stop_Time = stop_times_peak,
                      Horizon = total_samples_vec, 
                      Prop = stop_times_peak/total_samples_vec)
xtable(tab_dat)


# sample savings
print(total_samples - sum(taus$PEAK, na.rm = T))




