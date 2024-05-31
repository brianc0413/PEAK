library(mosaic)
library(confseq)

t <- 5000
alpha <- 0.05
upp_bounds <- rep(1, t)
low_bounds <- rep(0, t)

# our confidence interval
conf_int_peek <- function(S,mu_hat, alpha){
  
  if (length(S) < 10){
    return(c(0,1))
  }
  
  val <-findZeros( prod( 1+ (mu_hat - x) * (S - x)/0.26 ) - 1/alpha ~ x,
                  xlim = c(0,1))
  if (length(val) == 0){
    return(c(0,1))
  }
  if (length(val[,1]) == 1){
    return(c(0,1))
  }
  
  return( c( min(val[,1]), max(val[,1]) ) )
}

# bernoulli confidence interval (sub-bernoulli condition)
conf_int_bernoulli <- function(S, alpha){
  bounds <- bernoulli_confidence_interval(sum(S), length(S), alpha, t_opt = 2500, alpha_opt = alpha)
  return(c(bounds$lower, bounds$upper))
}

# empirical bernstein bounds
conf_int_bernstein <- function(S, mu_hat, var_hat, alpha){
  t <- length(S)
  lambda <- pmin(1/2, sqrt( 2*log(2/alpha) / (var_hat * 1:t * log(1+1:t)) ) )
  
  center <- sum(lambda * S) / sum(lambda)
  radius <- (log(2/alpha) + sum( (S - mu_hat)^2 * (-log(1-lambda) - lambda) ) ) / sum(lambda)
  
  return( c( max(0, center - radius) , min(1, center + radius)) )
}


# capital process bounds
conf_int_hedged <- function(S, mu_hat, var_hat, alpha){
  m <- 0:999/1000 + 0.0005
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


# mixture distirbution generator 1
generate_sample <- function(alpha, beta, epsilon){
  
  sample <- NA
  
  if (runif(1) < epsilon){
    sample <- 1
  }
  else{
    sample <- rbeta(1, alpha, beta)
  }
  
  return(sample)
}

# mixture distirbution generator 2
generate_sample_2 <- function(alpha, beta, epsilon){
  
  sample <- NA
  
  if (runif(1) < epsilon){
    sample <- 0
  }
  else{
    sample <- rbeta(1, alpha, beta)
  }
  
  return(sample)
}

generate_sample_weird <- function(n){
  random_thing <- runif(1)
  
  if (random_thing <= 1/3){
    return(runif(1))
  }
  
  if (random_thing > 1/3 && random_thing <= 2/3){
    return(rbeta(1,1,1))
  }
  
  if (random_thing > 2/3){
    return(rbinom(1, 1, 0.5))
  }
}



# naive bounds (use subgaussian factor = 1/4 here), predictable plug-in with correct factor
conf_int_prplh <- function(S, alpha){
  t <- length(S)
  lambda <- pmin(1, sqrt(8 * log(2/alpha) / (1:t * log(1:t+1)) ) )
  
  center <- sum(lambda * S)/sum(lambda)
  radius <- (log(2/alpha) + sum(lambda^2/8))/sum(lambda)
  
  return(c(center -radius, center + radius))
}

update_bound <- function(a, b){
  lower <- max(a[1], b[1])
  upper <- min(a[2], b[2])
  
  return(c(lower, upper))
}


sum_peek_conf_int <- data.frame(lower = rep(0, t+1), upper = rep(1, t+1))
sum_bernoulli_conf_int <- data.frame(lower = rep(0,t+1), upper = rep(1,t+1))
sum_bernstein_conf_int <- data.frame(lower = rep(0,t+1), upper = rep(1,t+1))
sum_hedged_conf_int <- data.frame(lower = rep(0,t+1), upper = rep(1,t+1))
sum_prplh_conf_int <- data.frame(lower = rep(0,t+1), upper = rep(1,t+1))


set.seed(2024)
for (j in 1:5){
  
  S <- c()
  var_hat <- c(1/4)
  mu_hat <- c(1/2)
  t <- 5000
  
  peek_conf_int <- data.frame(lower = rep(0, t+1), upper = rep(1, t+1))
  bernoulli_conf_int <- data.frame(lower = rep(0,t+1), upper = rep(1,t+1))
  bernstein_conf_int <- data.frame(lower = rep(0,t+1), upper = rep(1,t+1))
  hedged_conf_int <- data.frame(lower = rep(0,t+1), upper = rep(1,t+1))
  prplh_conf_int <- data.frame(lower = rep(0,t+1), upper = rep(1,t+1))
  
  
  for (i in 1:t){
    #S <- c(S, rbeta(1,1/4,1))
    #S <- c(S, generate_sample(1, 4, 0.05))
    S <- c(S, generate_sample_weird(1))
    print(length(S))
    print(length(mu_hat))
    print(length(var_hat))
    
    # peek conf int
    new_int <- update_bound(c(peek_conf_int[i,1], peek_conf_int[i,2]), 
                            conf_int_peek(S, mu_hat, alpha) )
    peek_conf_int[i+1,1] <- new_int[1]
    peek_conf_int[i+1, 2] <- new_int[2]
    
    # bernoulli conf int
    new_int <- update_bound(c(bernoulli_conf_int[i,1], bernoulli_conf_int[i,2]), 
                            conf_int_bernoulli(S, alpha) )
    bernoulli_conf_int[i+1,1] <- new_int[1]
    bernoulli_conf_int[i+1, 2] <- new_int[2]
    
    # bernstein
    new_int <- update_bound(c(bernstein_conf_int[i,1], bernstein_conf_int[i,2]), 
                            conf_int_bernstein(S, mu_hat, var_hat, alpha) )
    bernstein_conf_int[i+1,1] <- new_int[1]
    bernstein_conf_int[i+1, 2] <- new_int[2]
    
    # hedged
    new_int <- update_bound(c(hedged_conf_int[i,1], hedged_conf_int[i,2]), 
                            conf_int_hedged(S, mu_hat, var_hat, alpha) )
    hedged_conf_int[i+1,1] <- new_int[1]
    hedged_conf_int[i+1, 2] <- new_int[2]
    
    
    # prplh
    new_int <- update_bound(c(prplh_conf_int[i,1], prplh_conf_int[i,2]), 
                            conf_int_prplh(S, alpha) )
    prplh_conf_int[i+1,1] <- new_int[1]
    prplh_conf_int[i+1, 2] <- new_int[2]
    
    
    # at the very end
    mu_hat <- c(mu_hat, (sum(S) + 1/4)/(i+1) )
    var_hat <- c(var_hat,  (1/4 + sum( (S-mu_hat[i+1])^2 ) )/(i+1)  )
    
  }
  
  sum_peek_conf_int <- sum_peek_conf_int + peek_conf_int
  sum_bernoulli_conf_int <- sum_bernoulli_conf_int + bernoulli_conf_int
  sum_bernstein_conf_int <- sum_bernstein_conf_int + bernstein_conf_int
  sum_hedged_conf_int <- sum_hedged_conf_int + hedged_conf_int
  sum_prplh_conf_int <- sum_prplh_conf_int + prplh_conf_int
  
}


avg_peek <- sum_peek_conf_int/5
avg_bernoulli <- sum_bernoulli_conf_int/5
avg_bernstein<- sum_bernstein_conf_int/5
avg_hedged <- sum_hedged_conf_int/5
avg_prplh <- sum_prplh_conf_int/5

plot(x = 0:5000, y = avg_peek$lower, type = "l", lwd=2, col = "blue", ylim= c(0.1, 0.8))+
  lines(x=0:5000, y = avg_peek$upper, type = "l", lwd = 2, col = "blue")+
  lines(x=0:5000, y = avg_bernoulli$lower, type = "l", lwd = 2, col = "red")+
  lines(x=0:5000, y = avg_bernoulli$upper, type = "l", lwd = 2, col = "red")+
  lines(x=0:5000, y = avg_bernstein$lower, type = "l", lwd = 2, col = "orange")+
  lines(x=0:5000, y = avg_bernstein$upper, type = "l", lwd = 2, col = "orange")+
  lines(x=0:5000, y = avg_hedged$lower, type = "l", lwd = 2, col = "purple")+
  lines(x=0:5000, y = avg_hedged$upper, type = "l", lwd = 2, col = "purple")+
  lines(x=0:5000, y = avg_prplh$lower, type = "l", lwd = 2, col = "green")+
  lines(x=0:5000, y = avg_prplh$upper, type = "l", lwd = 2, col = "green") + 
  legend("topright", legend = c("PEEK", "Sub-Bernoulli", "Emp. Bern.", "Hedged", "PrPlH" ), col = c("blue", "red", "orange", "purple", "green"), cex=1, lwd=3)

plot(x = 0:5000, y = log(avg_peek$upper -avg_peek$lower), type = "l", lwd=2, col = "blue", ylim = c(-2, 0), xlab = "t", ylab = "Log Volume of Confidence Sequence", main = "Mixture Distribution Confidence Sequence Volume")+
  lines(x=0:5000, y = log(avg_bernoulli$upper -avg_bernoulli$lower), type = "l", lwd = 2, col = "red")+
  lines(x=0:5000, y = log(avg_bernstein$upper - avg_bernstein$lower), type = "l", lwd = 2, col = "orange")+
  lines(x=0:5000, y = log(avg_hedged$upper-avg_hedged$lower), type = "l", lwd = 2, col = "purple")+
  lines(x=0:5000, y = log(avg_prplh$upper-avg_prplh$lower), type = "l", lwd = 2, col = "green") +
  legend("topright", legend = c("PEEK", "Sub-Bernoulli", "Emp. Bern.", "Hedged", "PrPlH" ), col = c("blue", "red", "orange", "purple", "green"), cex=1, lwd=3)



