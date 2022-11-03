source("requirements.R")

fns <- c("output/spatiophylogenetic_modelling/featurewise_specific_reruns/GB129_kappa_2_sigma_1.15_pcprior0.1.qs", "output/spatiophylogenetic_modelling/featurewise_specific_reruns/GB197_kappa_2_sigma_1.15_pcprior0.1.qs")


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
  cat(paste0("I'm on ", fn, ", i.e. index ", index, " out of ", length(fns), ".\n"))
  
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

posteriors_df_reruns <- posteriors_df %>% 
  as.data.frame() %>% 
  mutate(Feature_ID = str_extract(fn, pattern = "[:alnum:]*")) %>% 
  mutate(Feature_ID = paste0(Feature_ID, "_rerun"))

posteriors_df <- read_tsv("output/spatiophylogenetic_modelling/featurewise/posteriors_df.tsv") %>% 
  full_join(posteriors_df_reruns) %>% 
  dplyr::select(Feature_ID, spatial = `Precision for spatial_id_in_dual`, phylogenetic = `Precision for phylo_id_in_dual`) %>% 
  group_by(Feature_ID) %>% 
  summarise(mean_phylogenetic = mean(phylogenetic),
            mean_spatial = mean(spatial),
            error_phylogenetic = sd(phylogenetic),
            error_spatial = sd(spatial)) %>% 
  filter(str_detect(Feature_ID, "GB129") | str_detect(Feature_ID, "GB197"))

