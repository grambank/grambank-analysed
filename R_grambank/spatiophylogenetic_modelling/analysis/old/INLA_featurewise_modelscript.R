# This script tests the following:
# Pagel's lambda
# Phylo only
# Spatial only
# Dual model

source("requirements.R")

#set this as 1 if you're just running this script on 50 lgs over 3 features to debug. Otherwise set to 0.
debug_run = 0

#should the scripts output an rds file for each output from INLA::inla() or not? Set to 0 if not, 1 if yes.
save_RDS_featurewise <- 0

#load funtions
source("global_variables.R")
source("spatiophylogenetic_modelling/analysis/INLA_parameters.R")

#dir setup
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

#load functions
source('spatiophylogenetic_modelling/analysis/functions/varcov_spatial.R')
source('spatiophylogenetic_modelling/analysis/functions/strip_inla.R')


#wrangle CLI input
args = commandArgs(trailingOnly = TRUE)

if(length(args) != 0){
  start = args[1]
  end <- args[2]
  range <- start:end
} else {
  start <- 1
  end <- 40
  range <- start:end
}


time <- as.character(Sys.time())

sink(file = paste0(OUTPUTDIR , "INLA_log_range_", start, "_", end,"_", time, ".txt" ), split = T)

#set up features to go to go through
#checking which ones have an qs file already and excluding them


#loading inputs
cat("\n###\nLoading covariance matrices...\n")

precision_matrices_fn <- "output/spatiophylogenetic_modelling/precision_matrices.RDS"
if(!(file.exists(precision_matrices_fn))){
  source("spatiophylogenetic_modelling/analysis/simulations/make_precisionmatrices.R")}

precision_matrices = readRDS(precision_matrices_fn)
phylo_prec_mat = precision_matrices$phylogenetic_precision
spatial_prec_mat = precision_matrices$spatial_precision

cat("\n###\nDone with covariance matrices.\n")

tree_fn <- "output/spatiophylogenetic_modelling/processed_data/EDGE_pruned_tree.tree"
if(!(file.exists(tree_fn))){
  source("spatiophylogenetic_modelling/processing/pruning_EDGE_tree.R")}
tree = read.tree(tree_fn)

tree_tips_df <- tree$tip.label %>% 
  as.data.frame() %>% 
  rename("Language_ID"= ".")

#reading in GB
GB_filename <- file.path("output", "GB_wide", "GB_cropped_for_missing.tsv")
if (!file.exists(GB_filename)) { 
  source("make_wide.R")
  source("make_wide_binarized.R")
  source("impute_missing_values.R")}		
GB <- tree_tips_df %>% 
  inner_join(  read.delim(GB_filename, sep = "\t"), by = "Language_ID")

## Since we are using a sparse phylogenetic matrix, we need to math taxa to the correct
## rows in the matrix
phylo_id = match(tree$tip.label, rownames(phylo_prec_mat))
GB$phylo_id = phylo_id

## Other effects are in the same order they appear in the dataset. 
GB$spatial_id = 1:nrow(spatial_prec_mat)
GB$obs_id = 1:nrow(spatial_prec_mat)

#features to loop over
features <- GB %>% 
  dplyr::select(-Language_ID) %>% 
  colnames() 

features <- features[range]

#### INLA Set up ####
## See the script make_precisionmatrices.R for details on how these were created
## Both matrices have their variance scaled to 1. 

index <- 0

for(feauture in features){
  #  feature <- features[1]
  cat(paste0("I'm on ", feature, " and the time is ", Sys.time(), ".\n"))
  
cat(paste0("done with set-up on ", feature, " and the time is ", Sys.time(), ".\n"))

#### Phylo Only model ####

## The phylogenetic model contains a additive effects for the phylogenetic component
## and a random effect for the residual variance which we have fixed to 1.
## We additionally calculate WAIC scores for this and each model for comparative purposes. 


formula <- eval(substitute(this_feature ~
                             f(phylo_id,
                               model = "generic0",
                               Cmatrix = phylo_prec_mat,
                               hyper = pcprior) +
                             f(obs_id, model = "iid", hyper = obs_hyper),
                            list(this_feature=as.name(feature))))

lambdaonly_model = inla(formula = formula,
                        family = "binomial",
                        control.compute = list(waic = TRUE), 
                        data = GB)

cat(paste0("Finished running phylo-only on ", feature, " and the time is ", Sys.time(), ".\n"))

#### Spatial Only model ####
## The spatial model is set up the same way as the phylogenetic model
## only containing the spatial precision matrix. 



formula <- eval(substitute(this_feature ~
                             f(spatial_id,
                               model = "generic0",
                               Cmatrix = spatial_prec_mat,
                               hyper = pcprior) +
                             f(obs_id, model = "iid", hyper = obs_hyper),
                           list(this_feature=as.name(feature))))

spatialonly_model = inla(formula = formula,
                         control.compute = list(waic = TRUE),
                         family = "binomial",
                         data = GB)

cat(paste0("Finished running spatial only on ", feature, " and the time is ", Sys.time(), ".\n"))

# #### Dual Model ####
# 
# # The dual model contains the phylogenetic and spatial precision models using the same
# # priors and scaled precision matrices. As well as a single effect for residuals (constrained to 1)

formula <-  eval(substitute(this_feature ~
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

dual_model = inla(formula = formula,
                  control.compute = list(waic = TRUE),
                  family = "binomial",
                  data = GB)

cat(paste0("Finished running dual on ", feature, " and the time is ", Sys.time(), ".\n"))

#### Save output ####
saved_file =   paste0(OUTPUTDIR, feature, ".qs")

phylogeny_only = strip_inla(lambdaonly_model)
cat(paste0("Finished running strip inla on phylo only  ", feature, " and the time is ", Sys.time(), ".\n"))

spatial_only = strip_inla(spatialonly_model)

cat(paste0("Finished running strip inla on spatial only  ", feature, " and the time is ", Sys.time(), ".\n"))

dual_model = strip_inla(dual_model)

cat(paste0("Finished running strip inla on dual  ", feature, " and the time is ", Sys.time(), ".\n"))

## For each model we strip out the information we need to calculate heritability scores
## to save on storage space. 
model_outputs = list(
  phylogeny_only,
  spatial_only,
  dual_model)

qs::qsave(model_outputs, 
          file = saved_file)

beep()

index <- index + 1

cat(paste0("I've finished ", feature, " and the time is ", Sys.time(), ".\n That means I'm ", round((index/length(range)), 2)*100, "% done.\n"))
rm(model_outputs, lambdaonly_model, spatialonly_model, dual_model)
}

sink(file = NULL)