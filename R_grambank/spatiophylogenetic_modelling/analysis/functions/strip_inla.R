object = phylo_only_model

strip_inla = function(object){
  
  # try to get the posterior of the ICC effect
  icc_posterior <-  try(get_iccposterior(object, n = 100))
  
  if ("try-error" %in% class(icc_posterior )) {
    cat(paste0("I failed in running strip_inla, most likely there's something wrong with the posteriors.\n"))
    icc_posterior  <- NULL
  }
  
  # save the WAIC score
  waic = object$waic

    if(is.null(waic)){
    cat(paste0("The waic score was NULL.\n")) 
    }
  
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
