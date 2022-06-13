strip_inla = function(object){
  
  # try to get the posterior of the ICC effect
  icc_posterior <-  try(get_iccposterior(object, n = 100))
  
  if (class(icc_posterior ) == "try-error") {
    icc_posterior  <- NULL
  }
  
  # save the WAIC score
  waic = object$waic
  
  list(icc_posterior = icc_posterior,
       waic = waic)
}

get_iccposterior = function(object, n = 100){
 
  hyper_sample = inla.hyperpar.sample(n = n,
                       result = object)
  
  binomial_error = pi^2 / 3
  
  if(ncol(hyper_sample) == 1){
    sigma = 1 / hyper_sample
    
    posterior = sigma / (sigma + 1 + binomial_error)
  }
  
  if(ncol(hyper_sample) == 2){
    sigma_1 = 1 / hyper_sample[,1]
    sigma_2 = 1 / hyper_sample[,2]
    
    posterior_1 = sigma_1 / (sigma_1 + sigma_2 + 1 + binomial_error)
    posterior_2 = sigma_2 / (sigma_1 + sigma_2 + 1 + binomial_error)
    
    posterior = cbind(posterior_1, posterior_2)
    colnames(posterior) = colnames(hyper_sample)
  }
  posterior
}
