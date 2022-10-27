#this script is to be applied to the output of INLA_all_models.R 
#This script is written by Sam Passmore, with minor adjustments by Hedvig Skirgård.

#The output could be either based on the real data or the simulated. Set the variable sim_or_real to indicate which is relevant. You can specify this from CLI if you want. If you don't specify, the default is real.

source("requirements.R")

args = commandArgs(trailingOnly = TRUE)
if(length(args) != 0){
  sim_or_real <- args[1] #you'd set this by going Rscript script.R "real"
} else { #if you're running this script chunkwise in Rstudio or similar instead of via command line, you'll read in the parameters this way:
  sim_or_real <- "real"
}

if(sim_or_real == "sim"){ #running the simulated data
  fns <- list.files("output/spatiophylogenetic_modelling/simulated_output/", pattern = "*.qs", full.names = T)
}

if(sim_or_real == "real"){ #running the on the real data
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
  #  hyper_sample = qs[[5]][[1]]
  
  
  #for the single models
  if(ncol == 1){
    
    sigma = 1 / hyper_sample
    
    posterior = sigma / (sigma + 1 + binomial_error) %>% 
    as.data.frame()
    colnames(posterior) = paste0(colnames(hyper_sample),"_in_single")
  }
  
  #for the dual models
  if(ncol == 2){
    sigma_1 = 1 / hyper_sample[,1]
    sigma_2 = 1 / hyper_sample[,2]
    
    posterior_1 = sigma_1 / (sigma_1 + sigma_2 + 1 + binomial_error)
    posterior_2 = sigma_2 / (sigma_1 + sigma_2 + 1 + binomial_error)
    
    posterior = cbind(posterior_1, posterior_2)
    colnames(posterior) = paste0(colnames(hyper_sample),"_in_dual")
    
  }
  if(ncol == 3){
    if(!is.null(hyper_sample)){

    sigma_1 = 1 / hyper_sample[,1]
    sigma_2 = 1 / hyper_sample[,2]
    sigma_3 = 1 / hyper_sample[,3]
    
    posterior_1 = sigma_1 / (sigma_1 + sigma_2 + sigma_3 + 1 + binomial_error)
    posterior_2 = sigma_2 / (sigma_1 + sigma_2 +sigma_3 + 1 + binomial_error)
    posterior_3 = sigma_3 / (sigma_1 + sigma_2 +sigma_3 + 1 + binomial_error)
    
    posterior = cbind(posterior_1, posterior_2, posterior_3)
  }else{
    posterior <- data.frame(matrix(ncol = 3, nrow = 100))
  }
    colnames(posterior) = c("Precision for spatial_id_in_trial", "Precision for phylo_id_in_trial",
                            "Precision for AUTOTYP_area_id_iid_model_in_trial")

  }
  posterior
}

#empty df to bind to

posteriors_df <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(posteriors_df) <- c( #"Precision for phylo_id_in_single"      ,
                              #"Precision for spatial_id_in_single"     ,          
                              #"Precision for AUTOTYP_area_id_iid_model_in_single", 
                              "Precision for spatial_id_in_dual"   ,              
                              "Precision for phylo_id_in_dual"          ,   
                              "Precision for spatial_id_in_trial"  ,              
                              "Precision for phylo_id_in_trial"               ,    
                              "Precision for AUTOTYP_area_id_iid_model_in_trial", "fn")

posteriors_df <- posteriors_df %>% 
  mutate_all(as.numeric)

posteriors_df$fn <- as.character() 

index <- 0
for(fn in fns){
#fn <- fns[49]

index <- index + 1
qs <- qs::qread(fn)
fn <- basename(fn) %>% str_replace_all(".qs", "")
cat(paste0("I'm on ", fn, ", i.e. index ", index, ".\n"))

#phylo_only
#hyper_sample_phylo_only_posterior <- get_icc_posterior(hyper_sample = qs[[1]][[1]], ncol = 1)

#spatial only
#hyper_sample_spatial_only_posterior <- get_icc_posterior(hyper_sample = qs[[2]][[1]], ncol = 1)

#autotyp-area
#hyper_sample_autotyp_area_only_posterior <- get_icc_posterior(hyper_sample = qs[[3]][[1]], ncol = 1)

#dual
hyper_sample_dual_posterior <- get_icc_posterior(hyper_sample = qs[[1]][[1]], ncol = 2)

#trial
hyper_sample_trial_posterior <- get_icc_posterior(hyper_sample = qs[[2]][[1]], ncol = 3)

posteriors_df_spec <- cbind(#hyper_sample_phylo_only_posterior, 
                       #hyper_sample_spatial_only_posterior, 
                       #hyper_sample_autotyp_area_only_posterior, 
                       hyper_sample_dual_posterior, 
                       hyper_sample_trial_posterior) %>% 
  as.data.frame() %>% 
  mutate(fn = fn)


posteriors_df <- posteriors_df %>% 
  full_join(posteriors_df_spec, by = c(#"Precision for phylo_id_in_single", 
                                       #"Precision for spatial_id_in_single", 
                                       #"Precision for AUTOTYP_area_id_iid_model_in_single",
                                       "Precision for spatial_id_in_dual", 
                                       "Precision for phylo_id_in_dual", 
                                       "Precision for spatial_id_in_trial", 
                                       "Precision for phylo_id_in_trial",
                                       "Precision for AUTOTYP_area_id_iid_model_in_trial", 
                                       "fn"))
}

posteriors_df %>% 
  write_tsv("output/spatiophylogenetic_modelling/featurewise/posteriors_df.tsv", na = "")