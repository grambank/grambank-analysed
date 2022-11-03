# This script tests the following:
# Dual model (spatial + phylo)
# Trial (spatial + phylo + AUTOTYP-area)

#this script re-run specific models. the features are listed below:
features <- c("GB197", "GB129")

OUTPUTDIR <- "output/spatiophylogenetic_modelling/featurewise_specific_reruns"
if(!dir.exists(OUTPUTDIR)){
  dir.create(OUTPUTDIR)
}

#load parameters
source("global_variables.R")
source("spatiophylogenetic_modelling/analysis/INLA_parameters.R")

pcprior_vec <- c(prior_ten_percent)
cat(paste0("I'm using just one pc prior. It is: \n", 
           pcprior_vec[[1]][2], "\n"
))
pcprior <- pcprior_vec[[1]]$param

source("requirements.R")

source("set_random_seed.R")

beep <- 0
if(beep == 1){
  h_load("beepr")
  beep(11) #starting sound
}


  sim_or_real <- "real"
  start <- 1
  end <- 40
  range <- start:end
  prec_matrices <- NULL
  pcprior_choice <- "default"

df_fn <- "output/GB_wide/GB_cropped_for_missing.tsv"

df <- readr::read_tsv(file =   df_fn,show_col_types = F)
  

#load functions
source('spatiophylogenetic_modelling/analysis/functions/varcov_spatial.R')
source('spatiophylogenetic_modelling/analysis/functions/strip_inla.R')
source("spatiophylogenetic_modelling/install_inla.R")


cat("\n###\nLoading covariance matrices...\n")

#using the defualt model paramters
precision_matrices_fn <- "output/spatiophylogenetic_modelling/processed_data/precision_matrices_kappa_2_sigma_1.15.RDS"

precision_matrices = readRDS(precision_matrices_fn)
phylo_prec_mat = precision_matrices$phylogenetic_precision
spatial_prec_mat = precision_matrices$spatial_precision

cat("\n###\nDone with covariance matrices.\n")

#reading in AUTOTYP-area
if (!file.exists("output/non_GB_datasets/glottolog_AUTOTYP_areas.tsv")) { s
  source("unusualness/processing/assigning_AUTOTYP_areas.R") }		
autotyp_area <- read.delim("output/non_GB_datasets/glottolog_AUTOTYP_areas.tsv", sep = "\t") %>%
  dplyr::select(Language_ID, AUTOTYP_area_id_iid_model = AUTOTYP_area)

#subsetting data to only entries which are also in the tree in case other slunk along
data <- lgs_in_analysis %>% 
  inner_join(df, by = "Language_ID") %>% 
  left_join(autotyp_area, by = "Language_ID")

x <- assert_that(all(data$Language_ID == lgs_in_analysis$Language_ID), msg = "Data doesn't match!")

## Since we are using a sparse phylogenetic matrix, we need to math taxa to the correct
## rows in the matrix
data$phylo_id = match(data$Language_ID, rownames(phylo_prec_mat))
data$spatial_id = match(data$Language_ID, rownames(spatial_prec_mat))
data$obs_id = 1:nrow(data)

  for(feature in features){
    #  feature <- features[1]
    cat(paste0("I'm on ", feature, " and the time is ", Sys.time(), ".\n"))
    cat(paste0("Precision matrix = ", basename(precision_matrices_fn), ".\n"))
    
    saved_file =   paste0(OUTPUTDIR, feature,"_",   substr(x = basename(precision_matrices_fn), 20, 37), "_pcprior",   pcprior[2], ".qs")
    
      # #### Dual Model ####
      # 
      # # The dual model contains the phylogenetic and spatial precision models using the same
      # # priors and scaled precision matrices. As well as a single effect for residuals (constrained to 1).
      
      formula <-  eval(substitute(expr = this_feature ~
                                    f(spatial_id,
                                      model = "generic0",
                                      Cmatrix = spatial_prec_mat,
                                      hyper = pcprior) +
                                    f(phylo_id,
                                      model = "generic0",
                                      Cmatrix = phylo_prec_mat,
                                      hyper = pcprior)  +
                                    f(obs_id, model = "iid", hyper = obs_hyper),
                                  list(this_feature=as.name(feature))))
      
      dual_model<- INLA::inla(formula = formula,
                              control.compute = list(waic = TRUE, cpo = TRUE, dic = TRUE), 
                              family = "binomial",
                              control.predictor=list(link=1),
                              data = data)
      
      dual_model %>% 
        qs::qsave(file = paste0(OUTPUTDIR, "/", feature, "_full_model_dual.qs"))
      
      
      # #### Trial Model ####
      # 
      # # The trial model contains the phylogenetic and spatial precision models using the same
      # # priors and scaled precision matrices. As well as a single effect for residuals (constrained to 1). In addition, there is also the categorical AUTOTYP-area included
      
      formula <-  eval(substitute(expr = this_feature ~
                                    f(spatial_id,
                                      model = "generic0",
                                      Cmatrix = spatial_prec_mat,
                                      hyper = pcprior) +
                                    f(phylo_id,
                                      model = "generic0",
                                      Cmatrix = phylo_prec_mat,
                                      hyper = pcprior)  +
                                    f(obs_id, model = "iid", hyper = obs_hyper) +
                                    f(AUTOTYP_area_id_iid_model, 
                                      hyper = pcprior,
                                      model = "iid"),
                                  list(this_feature=as.name(feature))))
      
      
      trial_model<- INLA::inla(formula = formula,
                               control.compute = list(waic = TRUE, cpo = TRUE, dic = TRUE), 
                               family = "binomial",
                               control.predictor=list(link=1),
                               data = data)
      
      trial_model %>% 
        qs::qsave(file = paste0(OUTPUTDIR, "/", feature, "_full_model_trial.qs"))
      
      ############
      #######strip interesting information from the INLA objects
      #########
      
      ## strip_inla function is located in another script. 
      #For each model we strip out the information we need to calculate heritability scores
      ## to save on storage space. 
      
      
      cat(paste0("Starting running strip inla on dual  of ", feature, " and the time is ", Sys.time(), ".\n"))
      
      dual_model_stripped = try(expr = {strip_inla(dual_model)})
      
      if (class(dual_model_stripped) == "try-error") {
        dual_model_stripped <- NULL
      }
      
      cat(paste0("Finished running strip inla on dual  of ", feature, " and the time is ", Sys.time(), ".\n"))
      
      cat(paste0("Starting running strip inla on trial model of  ", feature, " and the time is ", Sys.time(), ".\n"))
      
      trial_model_stripped <-try(expr = {strip_inla(trial_model)})
      
      if (class(trial_model_stripped) == "try-error") {
        trial_model_stripped <- NULL
      }
      
      cat(paste0("Finished running strip inla on trial model of  ", feature, " and the time is ", Sys.time(), ".\n"))
      
      
      
      #### Save output ####
      model_outputs = list(dual_model_stripped, 
                           trial_model_stripped)
      
      qs::qsave(model_outputs, 
                file = saved_file) }
    
    if(beep == 1){
      beep()
    }
   

