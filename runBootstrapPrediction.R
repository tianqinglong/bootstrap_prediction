source("simulate_censored_data.R")
source("parallel.R")
source("get_mle_pure_r_code.R")
# source("findMLEWrapper.R")

conditional_prob <- function(beta, eta, t_c, t_w){
  return( (pweibull(t_w, beta, eta) - pweibull(t_c, beta, eta))/
            (1-pweibull(t_c, beta, eta)) )
  }

compute_coverage_prob <- function(quan, num_at_risk, p_true){
  p4 <- pbinom(quan[4], num_at_risk, p_true)
  p3 <- pbinom(quan[3], num_at_risk, p_true)
  
  p2 <- dbinom(quan[2], num_at_risk, p_true)+pbinom(quan[2], num_at_risk, p_true, lower.tail = F)
  p1 <- dbinom(quan[1], num_at_risk, p_true)+pbinom(quan[1], num_at_risk, p_true, lower.tail = F)
  
  return(c(p1, p2, p3, p4))
}

compute_gpq_conditional_prob <- function(MLEs, BT_MLEs, t_c, t_w)
{
  MLE_mu <- log (MLEs[[2]])
  MLE_sigma <- 1/MLEs[[1]]
  
  BT_mu <- log (BT_MLEs[[2]])
  BT_sigma <- 1/BT_MLEs[[1]]
  
  GPQ_Beta <- 1 / ( MLE_sigma*MLE_sigma / BT_sigma )
  GPQ_Eta <- exp( MLE_mu + (MLE_mu - BT_mu)/BT_sigma*MLE_sigma )
  
  p_star_star <- conditional_prob(GPQ_Beta, GPQ_Eta, t_c, t_w)
  
  return (p_star_star)
}

N <- 5000
B <- 5000
beta <- 2
eta <- 1

pf1 <- 0.1
delta <- 0.2
r <- 5

n <- r/pf1
t_c <- qweibull(pf1, beta, eta)
t_w <- qweibull(pf1+delta, beta, eta)

p_true <- conditional_prob(beta, eta, t_c, t_w)

censored_data_list <- lapply(1:N, function(x) {return(simulate_data(r, pf1, beta, eta))})

t1 <- Sys.time()
mclapply(censored_data_list, function(x) {
  censored_data <- x
  
  mles <- get_Weibull_mle_R(censored_data)
  r_hat <- censored_data$Number_of_Failures
  
  p_hat <- conditional_prob(mles[[1]], mles[[2]], t_c, t_w)
  
  bootstrap_samples <- lapply(1:B , function(x) {simulate_data(r, pf1, mles[[1]], mles[[2]])})
  bootstrap_mles <- lapply(bootstrap_samples, function(x) {get_Weibull_mle_R(x)})
  
  value_vector <- NULL
  value_vector_2 <- NULL
  prob_vector <- NULL
  
  PB_predictive <- numeric( n-r_hat+1 )
  GPQ_predictive <- numeric( n-r_hat+1 )
  
  for (j in 1:B) {
    r_star <- bootstrap_samples[[j]][[1]]
    num_at_risk_star <- n - r_star
    
    p_star_star <- compute_gpq_conditional_prob(mles, bootstrap_mles[[j]], t_c, t_w)
    
    p_star <- conditional_prob(bootstrap_mles[[j]][[1]], bootstrap_mles[[j]][[2]], t_c, t_w)
    
    value_vector_temp <- ((0:num_at_risk_star) - num_at_risk_star*p_star)/sqrt(num_at_risk_star*p_star*(1-p_star))
    prob_vector_temp <- dbinom(0:num_at_risk_star, num_at_risk_star, p_hat)
    
    value_vector_temp_2 <- pbinom((0:num_at_risk_star),size = num_at_risk_star, prob = p_star)
    
    value_vector_2 <- c(value_vector_2, value_vector_temp_2)
    value_vector <- c(value_vector, value_vector_temp)
    prob_vector <- c(prob_vector, prob_vector_temp)
    
    PB_predictive <- PB_predictive + pbinom(0:(n-r_hat), n-r_hat, p_star)
    GPQ_predictive <- GPQ_predictive + pbinom(0:(n-r_hat), n-r_hat, p_star_star)
  }
  prob_vector <- prob_vector/B
  PB_predictive <- PB_predictive/B
  GPQ_predictive <- GPQ_predictive/B
  
  value_vector_order_2 <- order(value_vector_2)
  value_vector_sort_2 <- sort(value_vector_2)
  prob_vector_sort_2 <- prob_vector[value_vector_order_2]
  
  value_vector_order <- order(value_vector)
  value_vector_sort <- sort(value_vector)
  prob_vector_sort <- prob_vector[value_vector_order]
  
  cdf_value_2 <- cumsum(prob_vector_sort_2)
  cdf_value <- cumsum(prob_vector_sort)
  
  p_vector <- c(0.05, 0.1, 0.9, 0.95)
  
  boot_quantile_vector <- sapply(p_vector, function (x) {value_vector_sort[which(cdf_value > x)[1]]})
  boot_quantile_vector_2 <- sapply(p_vector, function (x) {value_vector_sort_2[which(cdf_value_2 > x)[1]]})
  
  quantile_vector <- boot_quantile_vector * sqrt( (n-r_hat)*p_hat*(1-p_hat) ) + (n-r_hat)*p_hat
  quantile_vector_2 <- qbinom(boot_quantile_vector_2, size = n - r_hat, prob = p_hat)
  
  quantile_vector[3:4] <- ceiling(quantile_vector[3:4])
  quantile_vector[1:2] <- floor(quantile_vector[1:2])
  
  quantile_vector_2[1:2] <- quantile_vector_2[1:2] - 1
  quantile_vector[quantile_vector < 0] <- 0
  quantile_vector_2[quantile_vector_2 < 0] <- 0
  # 
  # real_quantile <- sapply(quantile_vector, function (x) {
  #   out <- x
  #   if (out <= 0) {out <- 0}
  #   return(out)
  # })
  
  
  PB_quantile <- sapply(p_vector, function (x) {which(PB_predictive > x)[1]})
  if (is.na(PB_quantile[1])) {print(PB_predictive)}
  
  PB_quantile[1:2] <- PB_quantile[1:2] - 1
  PB_quantile[PB_quantile < 0] <- 0
  
  GPQ_quantile <- sapply(p_vector, function (x) {which(GPQ_predictive > x)[1]})
  if (is.na(GPQ_quantile[1])) {print(GPQ_predictive)}
  
  GPQ_predictive[1:2] <- GPQ_predictive[1:2] - 1
  GPQ_predictive[GPQ_predictive < 0] <- 0
  
  bootCP <- compute_coverage_prob(quantile_vector, n - r_hat, p_true)
  caliCP <- compute_coverage_prob(quantile_vector_2, n - r_hat, p_true)
  PBCP <- compute_coverage_prob(PB_quantile, n - r_hat, p_true)
  GPQCP <- compute_coverage_prob(GPQ_quantile, n - r_hat, p_true)
  
  return (list(BootP = quantile_vector, CaliP = quantile_vector_2, PB = PB_quantile,
               GPQ_PI = GPQ_quantile, BootCP = bootCP, CaliCP = caliCP, PBCP = PBCP, GPQCP = GPQCP) )
}) -> output
Sys.time() - t1
