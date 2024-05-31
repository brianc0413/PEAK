library(confseq)


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

# Confidence Intervals to test Against
# S = S_list[[i]], var_hat = vector/list of running variance estimates (offset by 1)
# 
conf_int_hedged <- function(S, var_hat, alpha, npts = 1000){
  m <- 0:(npts-1)/npts + 1/(2*npts)
  theta <- 1/2
  t <- length(S)
  lambda <- sqrt( 2*log(2/alpha) / (var_hat * 1:t * log(1+1:t)) ) 
  
  conf_int <- c()
  
  for (mean in m){
    lambda_m_plus <- pmin(lambda,0.5/mean)
    lambda_m_minus <- pmin(lambda, 0.5/(1-mean))
    k_plus <- theta * prod( 1 + lambda_m_plus * (S - mean))
    k_minus <- (1-theta) * prod( 1 - lambda_m_minus * (S - mean))
    
    
    if (max(k_plus, k_minus) < 1/alpha){
      conf_int <- c(conf_int, mean)
    }
  }
  
  return( c( min(conf_int),  max(conf_int) ) )
}




conf_int_bernoulli <- function(S, alpha){
  bounds <- bernoulli_confidence_interval(sum(S), length(S), alpha, t_opt = 2500, alpha_opt = alpha)
  return(c(bounds$lower, bounds$upper))
}



conf_int_base <- function(S, alpha){
  n <- 4
  
  mu <- mean(S)
  t <- length(S)
  rad <- sqrt( log(405.5 * n * t^(1.1)/alpha * log(405.5 * n * t^(1.1) / alpha) ) / (2*t) )
  
  return(c(mu - rad, mu + rad))
}



# termination condition
term_condition <- function(S_list, var_hats, method, alpha){
  
  term <- 0
  best <- NULL
  
  upp_bounds <- rep(1, 4)
  low_bounds <- rep(0, 4)
  
  for (i in 1:4){
    S <- S_list[[i]]
    var_hat <- var_hats[[i]]
    
    #print(c(length(S), length(var_hat)))
    
    if (method == 1){ # base
      bounds <- conf_int_base(S, alpha) # base already accounts for number of arms
      upp_bounds[i] <- bounds[2]
      low_bounds[i] <- bounds[1]
    }
    if (method == 2){ # bernoulli
      bounds <- conf_int_bernoulli(S, alpha/4) # divide by 4 to union bound
      upp_bounds[i] <- bounds[2]
      low_bounds[i] <- bounds[1]
    }
    if (method == 3){ # hedged
      bounds <- conf_int_hedged(S, var_hat, alpha/4, npts = 100) # divide by 4 to union bound
      upp_bounds[i] <- bounds[2]
      low_bounds[i] <- bounds[1]
      
    }
    
  }
  
  best_low_bounds <- which.max(low_bounds)
  
  if (low_bounds[best_low_bounds] >= max(upp_bounds[1:4 != best_low_bounds])){
    term <- 1
    best <- best_low_bounds
  }
  return(c(term, best))
}







data <- data.frame(method = NA,
                   stopping_time =NA,
                   best_arm = NA)

alpha <- 0.05

#mu <- c(0.29, 0.43, 0.57, 0.71)
mu <- (0.71 - 0.14*(1:4 - 1))/(0.29 + 0.14*(1:4-1))

nsims <- 100
for (method in 1:3){
  for (i in 1:nsims){
    
    set.seed(i) # set starting seed at beginning of simulation to make all sample paths the same
    
    #S_list <- list(rbinom(1,1,mu[1]), rbinom(1,1,mu[2]), rbinom(1,1,mu[3]), rbinom(1,1,mu[4]))
    #S_list <- list(rbeta(1,1,mu[1]), rbeta(1,1,mu[2]), rbeta(1,1,mu[3]), rbeta(1,1,mu[4]))
    S_list <- list(rbeta(1,1,mu[1]), 
                   rbeta(1,1,mu[2]), 
                   generate_sample(1, 0.43/0.52, 0.05), 
                   generate_sample_2(1, 0.24/0.71, 0.05)) # for beta setting
    t <- 4
    var_hats <- list(1/4, 1/4, 1/4, 1/4)
    
    term <- 0
    best <- NULL
    
    while (term == 0){
      which_arm <- arm_selector(S_list, alpha)
      h <- which_arm[1]
      l <- which_arm[2]
      
      var_hats[[h]] <- c(var_hats[[h]], (1/4 + sum( (S_list[[h]]-mean(S_list[[h]]))^2 ) )/(length(S_list[[h]])+1))
      var_hats[[l]] <- c(var_hats[[l]], (1/4 + sum( (S_list[[l]]-mean(S_list[[l]]))^2 ) )/(length(S_list[[l]])+1))
      
      
      #S_list[[h]] <- c(S_list[[h]], rbeta(1,1,mu[h]))
      #S_list[[l]] <- c(S_list[[l]], rbeta(1,1,mu[l]))
      
      #S_list[[h]] <- c(S_list[[h]], rbinom(1,1,mu[h]))
      #S_list[[l]] <- c(S_list[[l]], rbinom(1,1,mu[l]))
      
      
      # update samples
      if (h == 1 || h == 2){
        S_list[[h]] <- c(S_list[[h]], rbeta(1,1,mu[h]))
      }
      if (h == 3){
        S_list[[h]] <- c(S_list[[h]], generate_sample(1, 0.43/0.52, 0.05))
      }
      if (h == 4){
        S_list[[h]] <- c(S_list[[h]], generate_sample_2(1, 0.24/0.71, 0.05))
      }
      
      
      if (l == 1 || l == 2){
        S_list[[l]] <- c(S_list[[l]], rbeta(1,1,mu[l]))
      }
      if (l == 3){
        S_list[[l]] <- c(S_list[[l]], generate_sample(1, 0.43/0.52, 0.05))
      }
      if (l == 4){
        S_list[[l]] <- c(S_list[[l]], generate_sample_2(1, 0.24/0.71, 0.05))
      }
      
      #print(lapply(S_list, length))
      #print(lapply(S_list, length))
      
      ####### start timing
      terms <- term_condition(S_list, var_hats, method, alpha)
      #######  end timing
      
      term <- terms[1]
      best <- terms[2]
      
      t <- t+2
    }
    
    data <- rbind(data, c(method, t, best))
    print(data)
    
  }
  
}


mean(data$stopping_time[data$method == 1], na.rm = T)
mean(data$stopping_time[data$method == 2], na.rm = T)
mean(data$stopping_time[data$method == 3], na.rm = T)

sd(data$stopping_time[data$method == 1], na.rm = T)
sd(data$stopping_time[data$method == 2], na.rm = T)
sd(data$stopping_time[data$method == 3], na.rm = T)











