## things to try:

grid <- 501
pi <- seq(0, 1, length.out = grid)
better_arm <- rep(0, grid)
mu <- c(0.1, 0.5)
qs <- rep(0, grid)
#qs <- seq(mus[1], mus[2], length.out = 100)





c <- 0.26
# objective function

for (i in 1:grid){  
  #x_1 <- rbinom(1E5, 1, prob = mus[1])
  #x_2 <- rbinom(1E5, 1, prob = mus[2])
  
  # objective function for pi = pi[i]
  obj_new <- function(q){
    return( 10 * pi[i] * (mu[1] * ((2 * mu[1]-q - 1)/(c + (mu[1] - q) * (1-q))) + (1-mu[1]) * ((2 * mu[1]-q - 0)/(c + (mu[1] - q) * (-q)))) + 
            10  * (1-pi[i]) *  (mu[2] * (2 * mu[2]-q - 1)/(c + (mu[2] - q) * (1 - q)) + (1-mu[2]) * (2 * mu[2]-q - 0)/(c + (mu[2] - q)*(-q) ) ) )
  }
  
  # optimal q
  qs[i] <- findZeros(obj_new(x) ~ x, xlim = c(mu[1],mu[2]))[1,1]
  
  # optimal rate
  min_growth_rate[i] <- pi[i] * (mu[1]*(log(1+(mu[1]-qs[i])*(1 - qs[i])/c)) + (1-mu[1])*(log(1+(mu[1]-qs[i])*(-qs[i])/c)) ) + 
    (1-pi[i]) * (mu[2]* log(1 + (mu[2]-qs[i])*(1-qs[i])/c) + (1-mu[2]) * log(1 + (mu[2]-qs[i])*(-qs[i])/c ) ) 
  
  better_arm[i] <- ((mu[1]*(log(1+(mu[1]-qs[i])*(1 - qs[i])/c)) + (1-mu[1])*(log(1+(mu[1]-qs[i])*(-qs[i])/c)) ) > 
     (mu[2]* log(1 + (mu[2]-qs[i])*(1-qs[i])/c) + (1-mu[2]) * log(1 + (mu[2]-qs[i])*(-qs[i])/c ) )) 
  
}


pi[which.max(min_growth_rate[1:grid-1])]
mean(better_arm)
qs
min_growth_rate


# tells me to just pull arm 1 the entire time
plot(x= pi[1:grid], y= min_growth_rate)
plot(qs)


mu <- seq(0,1, length.out = 100)
m <- seq(0, 1, length.out = 100)

lr <- function(mu, m){
  
  return( mu*log(mu/m) + (1-mu)*log((1-mu)/(1-m)) )
  
}

lr_2 <- function(mu,m){
  return( mu*log(1 + (mu -m)*(1-m)/0.26) + (1-mu) * log(1 + (mu- m)*(0-m)/0.26)  )
}


dat <- expand.grid(mu, m)
colnames(dat) <- c("mu", "m")
dat$lr <- lr(dat$mu, dat$m)
dat$lr_2 <- lr_2(dat$mu, dat$m)
dat$ratio <- dat$lr_2 / dat$lr

ggp <- ggplot(dat, aes(mu, m)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = ratio)) + theme_bw()

ggp 

ggp2 <- ggplot(dat, aes(mu, m)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = lr_2))


ggp2


full_dat <- data.frame(mu, m)
lr_2(full_dat$mu, full_dat$m)

opt_policy_thr <- function(mu, m){
  lrs <- lr_2(mu, m)
  unnormalized <- 1/lrs
  return(unnormalized/(sum(unnormalized)))
}



opt_policy <- opt_policy_thr(mu = c(0.29, 0.43,.57, .71), m = c(0.2,0.3, 0.5, 0.6))
opt_policy


plot(lr_2(mu = 1:100/100, m = 1:100/100))

lr_2


