#This is a script for running binomial INLA over 113 binary Grambank features, with phylo and spatial effects.

#set this as 1 if you're just running this script on 50 lgs over 3 features to debug. Otherwise set to 0.
debug_run = 0

#should the scripts output an rds file for each output from INLA::inla() or not? Set to 0 if not, 1 if yes.
save_RDS_featurewise <- 0

source("global_variables.R")
source("spatiophylogenetic_modelling/analysis/INLA_parameters.R")

join_columns = c("2.5%", "50%", "97.5%", "Feature_ID", 
                 "effect", "waic", "model")

source("fun_def_h_load.R")

h_load(pkg = c("tidyverse", "reshape2","ape", "rlang", "assertthat" , "qs"))

#If the tree hasn't been prune yet - prune the tree :)
if (!file.exists("output/spatiophylogenetic_modelling/processed_data/jaeger_pruned.tree")) {
  source("spatiophylogenetic_modelling/processing/pruning_jagertree.R") 
}		

#library(INLA, quietly = T, warn.conflicts = F, verbose = F)

#make output dirs
if (!dir.exists("output/spatiophylogenetic_modelling/")) {
  dir.create("output/spatiophylogenetic_modelling/")
}

if (!dir.exists("output/spatiophylogenetic_modelling/featurewise/")) {
  dir.create("output/spatiophylogenetic_modelling/featurewise/")
}

if(debug_run == 1){
  OUTPUTDIR  <- file.path("output", "spatiophylogenetic_modelling","featurewise", "results_debug/")
} else{
  OUTPUTDIR <- file.path("output", "spatiophylogenetic_modelling","featurewise", "results/")
}

if (!dir.exists(  OUTPUTDIR )) {
  dir.create(  OUTPUTDIR )
}

if (!dir.exists(file.path(  OUTPUTDIR , "phylo_only/"))) {
  dir.create(file.path(  OUTPUTDIR , "phylo_only/"))
  dir.create(file.path(  OUTPUTDIR , "spatial_only/"))
  dir.create(file.path(  OUTPUTDIR , "autotyp_area_only/"))
  dir.create(file.path(  OUTPUTDIR , "dual_process_rdata/"))
  dir.create(file.path(  OUTPUTDIR , "trial_process_rdata/"))
  }		

cat("Starting INLA runs at", as.character(Sys.time()), ".\n")

#loading inputs
cat("\n###\nLoading covariance matrices...\n")
source("spatiophylogenetic_modelling/analysis/make_vcvs.R")

cat("\n###\nDone with covariance matrices.\n")

## Adding random effect ids
grambank_df = GB %>%
  left_join(tibble(Language_ID = rownames(phylo_prec_mat),
                  phy_id_generic = 1:nrow(phylo_prec_mat),
                  phy_id_iid_model = 1:nrow(phylo_prec_mat),
                  spatial_id_generic = 1:nrow(spatial_prec_mat),
                  spatial_id_iid_model = 1:nrow(spatial_prec_mat)), 
            by = "Language_ID") %>% 
  left_join(languages,  by = "Language_ID") %>% 
  rename(AUTOTYP_area_id_iid_model = AUTOTYP_area)

#################
###INLA LOOPS####
#################

#features to loop over
features <- GB %>% 
  dplyr::select(-Language_ID) %>% 
  colnames() 