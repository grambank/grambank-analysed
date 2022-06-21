#this script is to be applied to the output of INLA_all_models.R 
#This script is written by Sam Passmore, with minor adjustments by Hedvig Skirgård.

#The output could be either based on the real data or the simulated. Set the variable sim_or_real to indicate which is relevant.

sim_or_real <- "real"

if(sim_or_real == "sim"){
  fns <- list.files("output/spatiophylogenetic_modelling/simulated_output/", pattern = "*.qs", full.names = T)
}

if(sim_or_real == "real"){
  fns <- list.files("output/spatiophylogenetic_modelling/featurewise/", pattern = "*.qs", full.names = T) 
}

#This script takes the output of strip_inla() and calculates the appropriate posteriors. 

#The output of strip_inla() is calculated per effect in the INLA-models. We specify five models: phylo only, spatial only, autotyp-area only, dual (spatial + phylo) and trial (phylo, spatial and autoty-area). The strip_inla() function pulls out the waic, dic, cpo and pit scores as well as the hypersample. The hypersample is what we use here to calculate the posteriors. 

#The hypersample of the single models (i.e. phylo only, spatial only and autotyp-area) has only 1 column. The dual model has 2 and the trial 3. The names of the columns map onto the effects

#define the binominal error, which is uniform accross all models.
binomial_error = pi^2 / 3

#defining a function for pulling out the posterios appropriately
get_icc_posterior <- function(hyper_sample, ncol= NULL) {

  # commented out example for debugging:
  #  hyper_sample = hyper_sample_phylo_only
  
  #for the single models
  if(ncol == 1){
    
    sigma = 1 / hyper_sample
    
    posterior = sigma / (sigma + 1 + binomial_error) %>% 
    as.data.frame()
    
  }
  
  #for the dual models
  if(ncol == 2){
    sigma_1 = 1 / hyper_sample[,1]
    sigma_2 = 1 / hyper_sample[,2]
    
    posterior_1 = sigma_1 / (sigma_1 + sigma_2 + 1 + binomial_error)
    posterior_2 = sigma_2 / (sigma_1 + sigma_2 + 1 + binomial_error)
    
    posterior = cbind(posterior_1, posterior_2)
  }
  colnames(posterior) = colnames(hyper_sample)
  posterior
}


for(fn in fns){
#fn <- fns[112]

qs <- qs::qread(fn)

#phylo_only
hyper_sample_phylo_only_posterior <- get_icc_posterior(hyper_sample = qs[[1]][[1]], ncol = 1)

#spatial only
hyper_sample_spatial_only_posterior <- get_icc_posterior(hyper_sample = qs[[2]][[1]], ncol = 1)

#autotyp-area
hyper_sample_autotyp_area_only_posterior <- get_icc_posterior(hyper_sample = qs[[3]][[1]], ncol = 1)

#dual
hyper_sample_dual_posterior <- get_icc_posterior(hyper_sample = qs[[4]][[1]], ncol = 2)

#trial
hyper_sample_trial <- qs[[5]][1]
}
