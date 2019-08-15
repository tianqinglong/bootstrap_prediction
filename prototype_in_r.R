B <- 5000

r <- 15
pf1 <- 0.1
beta <- 2
eta <- 1

n <- r/pf1


censored_data <- simulate_data(r, pf1, beta, eta)

r_hat <- censored_data$Number_of_Failures

mles <- get_Weibull_mle_R(censored_data)

t1 <- Sys.time()
# simulate bootstrap samples
bootstrap_samples <- lapply(1:B , function(x) {simulate_data(r, pf1, mles[[1]], mles[[2]])})
bootstrap_mles <- lapply(bootstrap_samples, function(x) {get_Weibull_mle_R(x)})
Sys.time() - t1

pf_delta <- 0.1

t_first_censor <- qweibull(pf1, beta, eta)
t_second_censor <- qweibull(pf1 + pf_delta, beta, eta)

conditional_prob <- function(beta, eta, t_c, t_w){
  return( (pweibull(t_w, beta, eta) - pweibull(t_c, beta, eta))/
            (1-pweibull(t_c, beta, eta)) )
}

p_hat <- conditional_prob(mles$Shape_mle, mles$Scale_mle, t_first_censor, t_second_censor)

value_vector <- NULL
prob_vector <- NULL
for (i in 1:B)
{
  r_star <- bootstrap_samples[[i]][[1]]
  num_at_risk <- n - r_star
  
  p_star <- conditional_prob(bootstrap_mles[[i]][[1]], bootstrap_mles[[i]][[2]], t_first_censor,
                             t_second_censor)
  
  value_vector_temp <- ((0:num_at_risk) - num_at_risk*p_star)/sqrt(num_at_risk*p_star*(1-p_star))
  
  prob_vector_temp <- dbinom((0:num_at_risk), num_at_risk, p_hat)
  
  value_vector <- c(value_vector, value_vector_temp)
  prob_vector <- c(prob_vector, prob_vector_temp)
}

prob_vector <- prob_vector/B

# Get empirical distribution
value_vector_order <- order(value_vector)

value_vector_sort <- sort(value_vector)
prob_vector_sort <- prob_vector[value_vector_order]
cdf_value <- cumsum(prob_vector_sort)

p_vector <- c(0.05, 0.1, 0.9, 0.95)

boot_quantile_vector <- sapply(p_vector, function (x) {value_vector_sort[min(which(cdf_value > x))]})
quantile_vector <- boot_quantile_vector * sqrt( (n-r_hat)*p_hat*(1-p_hat) ) + (n-r_hat)*p_hat

quantile_vector[3:4] <- ceiling(quantile_vector[3:4])
quantile_vector[1:2] <- floor(quantile_vector[1:2])

real_quantile <- sapply(quantile_vector, function (x) {
  out <- x
  if (out < 0) {out <- 0}
  return(out)
}
)
