get_Weibull_mle_R <- function(censor_data) {
  
  n_minus_r <- censor_data[[4]] - censor_data[[1]]
  sobj <- Surv(time = c( censor_data[[3]],
                         
                         rep(censor_data[[2]], n_minus_r)
                         
  ),
  event = c( rep(1, censor_data[[1]] ),
             
             rep(0, n_minus_r)
  ),
  type = 'right'
  )
  
  sfit <- survreg(sobj~1, dist = 'weibull')

  return( list( Shape_mle = 1/sfit$scale,
                
                Scale_mle = as.numeric( exp(sfit$coefficients) )))
}