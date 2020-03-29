library(mvtnorm)
library(tmvtnorm)
library(MultiRNG)
library(doSNOW)
library(purrr)

##### input -----

sigma11 <- c(2, 4, 8)
sigma22 <- c(2, 4, 8)
sigma12 <- c(0, 2)
mu <- c(2, 4, 8)
mu_star <- c(2, 4, 8)
N <- c(20, 30, 50, 100, 500, 1000)
D <- c("n", "t", "u")

# simulazione ----
simulations <- expand.grid(sigma11, sigma22, sigma12, mu, mu, mu_star, mu_star, N, D)
names(simulations) <- c("s11", "s22", "s12", "m1", "m2", "m_star1", "m_star2", "n", "d")

# remove where sigma11 + sigma 22 < sigma12 * 2
simulations <- simulations[simulations$s11 + simulations$s22 > simulations$s12*2, ]


sim <- function(list){
  s11 <- list$s11
  s22 <- list$s22
  s12 <- list$s12
  m1 <- list$m1
  m2 <- list$m2
  m_star1 <- list$m_star1
  m_star2 <- list$m_star2
  n <- list$n
  d <- list$d
  
  sigma <- matrix(c(s11, s12, s12, s22), nrow=2)
  mu <- c(m1, m2)
  mu_star <- c(m_star1, m_star2)
  x <- rmvnorm(n, mean=mu, sigma=sigma)
  
  if(d == "n") x_star <- rmvnorm(n, mean=mu_star, sigma=sigma)
  if(d == "t") x_star <- rtmvnorm(n, mean=mu_star, sigma=sigma)
  if(d == "u") x_star <- draw.d.variate.uniform(n,2,sigma)
  
  T <- ((x_star[, 1] - colMeans(x)[1])/sqrt(var(x[, 1]))) - ((x_star[, 2] - colMeans(x)[2])/sqrt(var(x[, 2])))
  
  return(list("s11" = s11, 
              "s22" = s22, 
              "s12" = s12, 
              "m1" = m1, 
              "m2" = m2, 
              "m_star1" = m_star1,
              "m_star2" = m_star2, 
              "n" = n, 
              "d" = d, 
              "x" = x, 
              "x_star" = x_star,
              "T" = T))
  
}


run_sim <- function(simulations, fun = sim, multicore = TRUE){
  if(multicore){
    cores <- parallel::detectCores()
    cl <- makeSOCKcluster(cores)
    registerDoSNOW(cl)
    
    pb <- txtProgressBar(min=1, max=nrow(simulations), style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    
    result <- foreach(i=1:nrow(simulations), .packages=c("mvtnorm", "tmvtnorm", "MultiRNG"), .options.snow=opts)  %dopar% {
                fun(simulations[i, ])
              }
    close(pb)
    stopCluster(cl)
    
  }else{
    result <- lapply(c(1:nrow(simulations)), function(i) fun(simulations[i,])) 
  }
  
  return(result)
}


result <- run_sim(simulations, multicore = TRUE)

# risultati ----

tresult <- transpose(result)
results <- data.frame(s11 = unlist(tresult$s11), 
                      s22 = unlist(tresult$s22), 
                      s12 = unlist(tresult$s12), 
                      m1 = unlist(tresult$m1), 
                      m2 = unlist(tresult$m2), 
                      m_star1 = unlist(tresult$m_star1), 
                      m_star2 = unlist(tresult$m_star2), 
                      n = unlist(tresult$n), 
                      d = unlist(tresult$d),
                      mean_x1 = sapply(tresult$x, function(s) mean(s[,1])), 
                      mean_x2 = sapply(tresult$x, function(s) mean(s[,2])), 
                      mean_x_star1 = sapply(tresult$x_star, function(s) mean(s[,1])), 
                      mean_x_star2 = sapply(tresult$x_star, function(s) mean(s[,2])), 
                      mean_T = sapply(tresult$T, function(s) mean(s)), 
                      sd_T = sapply(tresult$T, function(s) sd(s))
                      )
results
