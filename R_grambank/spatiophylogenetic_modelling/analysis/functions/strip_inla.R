strip_inla = function(object){
  
#  object <- trial_model
  # try to get the posterior of the ICC effect
  icc_posterior <-  try(get_hyper_sample(object, n = 100))
  
  if ("try-error" %in% class(icc_posterior )) {
    cat(paste0("I failed in running strip_inla, most likely there's something wrong with the posteriors.\n"))
    icc_posterior  <- NULL
  }
  
  # save the WAIC score
  waic = object$waic
  cpo =  object$cpo$cpo
  pit = object$cpo$pit
  cpo_failure = object$cpo$failure

  list(hyper_sample = hyper_sample,
       waic = waic, cpi = cpo, pit = pit, cpo_failure = cpo_failure)
}

get_hyper_sample = function(object, n = 100){
 
  hyper_sample = inla.hyperpar.sample(n = n,
                       result = object)
hyper_sample
}
