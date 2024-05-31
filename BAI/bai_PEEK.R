library(mosaic)

# LUCB algorithm

arm_selector <- function(S_list, alpha){
  n <- length(S_list)
  upp_bounds <- rep(1, n)
  mus <- rep(0.5, n)
  
  for (i in 1:n){
    mus[i] <- mean(S_list[[i]])
    t <- length(S_list[[i]])
    rad <- sqrt( log(405.5 * n * t^(1.1)/alpha * log(405.5 * n * t^(1.1) / alpha) ) / (2*t) )
    upp_bounds[i] <- mus[i] + rad
  }
  
  # highest mean
  arm_1 <- which.max(mus)
  upp_bounds[arm_1] <- -1E10
  arm_2 <- which.max(upp_bounds)
  
  return(c(arm_1, arm_2))
}

# find_opt
f <- function(S_list, mus, m){
  val <- 0
  for (i in 1:4){
    val <- val + 0.25 * prod(1+ (S_list[[i]]-m[i])*(mus[[i]]-m[i])/0.26 )
  }
  return(val)
}

coeffs <- function(x, S_list, mus, a){
  return(prod(1+(S_list[[a]] - x)*(mus[[a]] - x)/0.26))
}


#f(S_list, mus, c(0,0,0,0))



# simulation for best arm identification using our approach
best_arm_list <- c()
stop_times <- c()

mu <- c(0.29, 0.43, 0.57, 0.71) # for bernoulli setting
#mu <- (0.71 - 0.14*(1:4 - 1))/(0.29 + 0.14*(1:4-1))

for (j in 1:100){ # number of simulations
  # means
  partitions <- rep(0,4)
  set.seed(j)
  S_list <- list(rbinom(1,1,mu[1]), rbinom(1,1,mu[2]), rbinom(1,1,mu[3]), rbinom(1,1,mu[4])) # for bernoulli setting
  #S_list <- list(rbeta(1,1,mu[1]), rbeta(1,1,mu[2]), rbeta(1,1,mu[3]), rbeta(1,1,mu[4])) # for beta setting
  mus <- list(1/2, 1/2, 1/2, 1/2)
  t <- 4
  alpha <- 0.05
  
  
  while (sum(partitions) < 3){
  #for (i in 1:300){
    which_arm <- arm_selector(S_list, alpha)
    h <- which_arm[1]
    l <- which_arm[2]
    
    # update mus
    mus[[h]] <- c(mus[[h]], (1/2+sum(S_list[[h]])) / (length(S_list[[h]]) +1) )
    mus[[l]] <- c(mus[[l]], (1/2+sum(S_list[[l]])) / (length(S_list[[l]]) +1) )
    
    # update samples
    S_list[[h]] <- c(S_list[[h]], rbinom(1,1,mu[h]))
    S_list[[l]] <- c(S_list[[l]], rbinom(1,1,mu[l]))
    
    #S_list[[h]] <- c(S_list[[h]], rbeta(1,1,mu[h]))
    #S_list[[l]] <- c(S_list[[l]], rbeta(1,1,mu[l]))
  
    
    # find the global minima in our set
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
    best_arm <- which(ranks == 4)
    
    # minimum e-value globally corresponds to partition with best arm.
    test_values <- rep(0, length(S_list))
    test_values[best_arm] <- f(S_list, mus, min_values)
    
    # for the second best arm, we only need to project second best, best
    second_best <- which(ranks == 3)
    f2 <- function(x){
      coeffs(x, S_list, mus, best_arm) * sum((2*x-mus[[best_arm]]-S_list[[best_arm]])/(0.26 + (S_list[[best_arm]]-x)*(mus[[best_arm]]-x))) + 
        coeffs(x, S_list, mus, second_best) * sum((2*x-mus[[second_best]]-S_list[[second_best]])/(0.26 + (S_list[[second_best]]-x)*(mus[[second_best]]-x)))
    }
    val <- findZeros(f2(x) ~ x, xlim=c(0,1))
    temp <- min_values
    temp[best_arm] <- val[1,1]
    temp[second_best] <- val[1,1]
    test_values[second_best] <- f(S_list, mus, temp)
  
    # for the third best arm, we need to project third best, second best, best 
    third_best <- which(ranks == 2)
    f3 <- function(x){
      coeffs(x, S_list, mus, best_arm) * sum((2*x-mus[[best_arm]]-S_list[[best_arm]])/(0.26 + (S_list[[best_arm]]-x)*(mus[[best_arm]]-x))) + 
        coeffs(x, S_list, mus, second_best) * (x < min_values[second_best]) * sum((2*x-mus[[second_best]]-S_list[[second_best]])/(0.26 + (S_list[[second_best]]-x)*(mus[[second_best]]-x))) + 
        coeffs(x, S_list, mus, third_best) * sum((2*x-mus[[third_best]]-S_list[[third_best]])/(0.26 + (S_list[[third_best]]-x)*(mus[[third_best]]-x))) 
    }
    val <- findZeros(f3(x) ~ x, xlim=c(0,1))
    temp <- min_values
    temp[best_arm] <- val[1,1]
    temp[second_best] <- val[1,1]
    temp[third_best] <- val[1,1]
    test_values[third_best] <- f(S_list, mus, temp)
    
    # for the fourth best arm, we project all values
    fourth_best <- which(ranks==1)
    f4 <- function(x){
      coeffs(x, S_list, mus, best_arm) * sum((2*x-mus[[best_arm]]-S_list[[best_arm]])/(0.26 + (S_list[[best_arm]]-x)*(mus[[best_arm]]-x))) + 
        coeffs(x, S_list, mus, second_best) * (x < min_values[second_best]) * sum((2*x-mus[[second_best]]-S_list[[second_best]])/(0.26 + (S_list[[second_best]]-x)*(mus[[second_best]]-x))) + 
        coeffs(x, S_list, mus, third_best) * (x < min_values[third_best]) * sum((2*x-mus[[third_best]]-S_list[[third_best]])/(0.26 + (S_list[[third_best]]-x)*(mus[[third_best]]-x))) + 
        coeffs(x, S_list, mus, fourth_best) * sum((2*x-mus[[fourth_best]]-S_list[[fourth_best]])/(0.26 + (S_list[[fourth_best]]-x)*(mus[[fourth_best]]-x)))
    }
    val <- findZeros(f4(x) ~ x, xlim=c(0,1))
    temp <- rep(val[1,1], 4)
    test_values[fourth_best] <- f(S_list, mus, temp)
    
    #print(test_values)
    
    # check to see which arms can be eliminated
    reject <- (test_values >= 1/alpha) # 0 if eliminated 
    partitions <- pmax(reject, partitions) # if eliminated before, doesnt matter
    
    # update number of pulls
    t <- t+2
  
  }
  best_arm_list <- c(best_arm_list, which(partitions == 0))
  stop_times <- c(stop_times, t)
  
  print(stop_times)
  print(best_arm_list)
}



## THIS
print(mean(stop_times))
print(sd(stop_times))




  