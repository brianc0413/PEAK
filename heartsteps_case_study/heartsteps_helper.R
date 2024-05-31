# helper for heartsteps dataset


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




term_condition <- function(S_list, var_hats, method, alpha){
  
  term <- 0
  best <- NULL
  
  upp_bounds <- rep(1, 2)
  low_bounds <- rep(0, 2)
  
  for (i in 1:2){
    S <- S_list[[i]]
    var_hat <- var_hats[[i]]
    
    #print(c(length(S), length(var_hat)))
    
    if (method == 1){ # base
      bounds <- conf_int_base(S, alpha) # base already accounts for number of arms
      upp_bounds[i] <- bounds[2]
      low_bounds[i] <- bounds[1]
    }
    if (method == 2){ # bernoulli
      bounds <- conf_int_bernoulli(S, alpha/2) # divide by 4 to union bound
      upp_bounds[i] <- bounds[2]
      low_bounds[i] <- bounds[1]
    }
    if (method == 3){ # hedged
      bounds <- conf_int_hedged(S, var_hat, alpha/2, npts = 100) # divide by 4 to union bound
      upp_bounds[i] <- bounds[2]
      low_bounds[i] <- bounds[1]
      
    }
    
  }
  
  best_low_bounds <- which.max(low_bounds)
  
  if (low_bounds[best_low_bounds] >= max(upp_bounds[1:2 != best_low_bounds])){
    term <- 1
    best <- best_low_bounds
  }
  return(c(term, best))
}


f <- function(S_list, mus, m){
  val <- 0
  for (i in 1:2){
    val <- val + 0.5 * prod(1+ (S_list[[i]]-m[i])*(mus[[i]]-m[i])/0.26 )
  }
  return(val)
}

coeffs <- function(x, S_list, mus, a){
  return(prod(1+(S_list[[a]] - x)*(mus[[a]] - x)/0.26))
}



