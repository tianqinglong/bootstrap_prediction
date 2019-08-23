library(survival)
library(rlist)

simulate_data <- function(r, Pt, shape, scale) {
  
  n <- round( r/Pt )
  y_vector <- numeric(n)
  
  x <- 0
  
  for (i in 0:(n-1)) 
    {
      u <- runif(1)
      y <- 1 - (1-x) * (1-u)^(1.0/(n-i))
      
      if(y > Pt) break
    
      y_vector[i+1] <- y
      x <- y
    }
  
  if (sum(y_vector > 0) < 2) {
    output <- simulate_data(r, Pt, shape, scale)
    return(output)
    }
    
    num_failure <- sum( y_vector > 0 )
    
    failure_times <- scale * ( -log( 1 - y_vector[ 1:num_failure ] ) ) ^ ( 1/shape )
    
    censor_time <- scale * ( -log( 1 - Pt ) ) ^ ( 1/shape )

    output <- list(Number_of_Failures = num_failure, Censor_Time = censor_time, Failure_Times = failure_times, Total_Number = n)
    
    return(output)
}
