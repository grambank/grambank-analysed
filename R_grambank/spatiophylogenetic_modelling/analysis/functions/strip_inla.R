strip_inla = function(object){ 
  
#  object <- trial_model
  # try to get the posterior of the ICC effect
  hyper_sample <-  try(get_hyper_sample(object))
  
  # save the WAIC score
  waic = object$waic
  cpo =  object$cpo$cpo
  pit = object$cpo$pit
  cpo_failure = object$cpo$failure
  ## get marginal likelihood, note this will need to have a normalising constant added to it later
  ## this normalising constant only needs to be calculated once for each precision matrix
  ## and is 0.5 * log(det), where det is the determinant of the precision matrix. If the model uses two
  ## precision matrices both constants need to be added.
  ## Let me know if you want to do the calculation here. It can be done but strip_inla() will
  ## get a bit more complicated because the constant can't be extracted from the model object,
  ## we will have to pass in the precision matrices as well (or the precalculated constants,
  ## which would be better).
  mlik <- object$mlik

  list(hyper_sample = hyper_sample,
       waic = waic, cpi = cpo, pit = pit, cpo_failure = cpo_failure, mlik = mlik)
}

get_hyper_sample = function(object, n = 100){
 
  hyper_sample = inla.hyperpar.sample(n = n,
                       result = object)
hyper_sample
}
