#This is a script for running binomial INLA over 113 binary Grambank features, with phylo and spatial effects.

#set this as 1 if you're just running this script on 50 lgs over 3 features to debug. Otherwise set to 0.
debug_run = 1

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

# Phylogenetic matrix is right dims #comment out if debugging swiftly
#x <- assert_that(nrow(spatial_covar_mat) == n_overlap_imputed_and_EDGE_tree, 
#                 msg = "The phylogeny has changed and will not match the data")

#### Set up model priors ####

## Taken from Dinnage et al. (2020):
### As in the manuscript, we use a “PC” prior, which stand for “Penalizing Complexity”.
### This is a standard prior developed by the developers of INLA, which is “weakly informative”.
### It places slightly more of the prior probability density on values close to zero,
### but has a “long tail”, which allows the data to push the parameter away from zero if
### there is good evidence (i.e. the likelihood of the data is higher).
### The PC prior has two parameters, p1 and p2: p2 is the proportion of the prior probability
### density that falls above values greater than p1.
###
### Given we are modelling (somewhat) Gaussian data that we have standardised to a variance of 1,
### we would not expect any random factor to have a variance greater than one.
### So we will set our prior to only have about 10% of its prior probability density
### above 1. We will guestimate the total variance possibly explained by the
### phylogenetic effect based on it’s diagonal entries, which are all equal
### to 2.7373648. So to get to 1, the scaling factor would have to be about 0.36

pcprior = list(prec =list(prior="pc.prec", param = c(1, 0.1)))

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

#subsetting for debugging code swiftly
if(debug_run == 1) {
features <- features[3:7]
}