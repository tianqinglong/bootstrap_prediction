dyn.load('mleFast.so')

findMLEFast <- function(censored_data)
{
  num_censor <- censored_data$Total_Number - censored_data$Number_of_Failures
  
  failure_time <- censored_data$Failure_Times
  
  censored_time <- rep(censored_data$Censor_Time, num_censor)
  
  mle <- .Call("findMLE", failure_time, censored_time)
  
  if (mle[1] < 0) { mle <- get_Weibull_mle_R(censored_data) }
  
  return(mle)
}
