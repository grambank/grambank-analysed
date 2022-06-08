# This script tests the following:
# Pagel's lambda
# Phylo only
# Spatial only
# Dual model

source("requirements.R")
beep <- 1

#set this as 1 if you're just running this script on 50 lgs over 3 features to debug. Otherwise set to 0.debug_run = 0

#wrangle CLI input
args = commandArgs(trailingOnly = TRUE)

if(length(args) != 0){
  sim_or_real <- args[1]
  start = args[2]
  end <- args[3]
  range <- start:end
} else { #if you're running this script chunkwise in Rstudio or similar instead of via command line, you'll read in the parameters this way:
  sim_or_real <- "sim"
  start <- 1
  end <- 40
  range <- start:end
}

if(sim_or_real == "sim"){
  df <- readr::read_tsv(file = "output/spatiophylogenetic_modelling/simulated_data/ simulated_data_df.tsv",show_col_types = F)
  OUTPUTDIR_top <- "output/spatiophylogenetic_modelling/simulated_output/"
  }
if(sim_or_real == "real"){
  df <- readr::read_tsv(file = "output/GB_wide/GB_cropped_for_missing.tsv",show_col_types = F)
  OUTPUTDIR <- "output/spatiophylogenetic_modelling/featurewise/"
}

#dir setup
#make output dirs
if (!dir.exists("output/spatiophylogenetic_modelling/")) {
  dir.create("output/spatiophylogenetic_modelling/")
}

if (!dir.exists(OUTPUTDIR)) {
  dir.create(OUTPUTDIR)
}

#load parameters
source("global_variables.R")
source("spatiophylogenetic_modelling/analysis/INLA_parameters.R")

#load functions
source('spatiophylogenetic_modelling/analysis/functions/varcov_spatial.R')
source('spatiophylogenetic_modelling/analysis/functions/strip_inla.R')

time <- as.character(Sys.time())

sink(file = paste0(OUTPUTDIR , "INLA_log_range_", start, "_", end,"_", time, ".txt" ), split = T)

#set up features to go to go through
#checking which ones have an qs file already and excluding them

#loading inputs

## See the script make_precisionmatrices.R for details on how these were created
## Both matrices have their variance scaled to 1. 

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

#reading in AUTOTYP-area
if (!file.exists("output/non_GB_datasets/glottolog_AUTOTYP_areas.tsv")) { source("unusualness/processing/assigning_AUTOTYP_areas.R") }		
autotyp_area <- read.delim("output/non_GB_datasets/glottolog_AUTOTYP_areas.tsv", sep = "\t") %>%
  dplyr::select(Language_ID, AUTOTYP_area_id_iid_model = AUTOTYP_area)

#subsetting data to only entries which are also in the tree in case other slunk along
data <- tree_tips_df %>% 
    inner_join(df, by = "Language_ID") %>% 
  left_join(autotyp_area, by = "Language_ID")

x <- assert_that(nrow(data) == Ntip(tree), msg = "Data doesn't match phylogeny!")

## Since we are using a sparse phylogenetic matrix, we need to math taxa to the correct
## rows in the matrix
phylo_id = match(tree$tip.label, rownames(phylo_prec_mat))
data$phylo_id = phylo_id

## Other effects are in the same order they appear in the dataset. 
data$spatial_id = 1:nrow(spatial_prec_mat)
data$obs_id = 1:nrow(spatial_prec_mat)

#features to loop over
features <- data %>% 
  dplyr::select(-Language_ID) %>% 
  colnames() 

features <- features[range]

index <- 0

for(feauture in features){
  #  feature <- features[1]
  cat(paste0("I'm on ", feature, " and the time is ", Sys.time(), ".\n"))
  
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
  
  phylo_only_model<- INLA::inla(formula = formula,
                          family = "binomial",
                          control.compute = list(waic = TRUE), 
                          control.predictor=list(link=1),
                          data = data)
  
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
  
  spatialonly_model<- INLA::inla(formula = formula,
                           control.compute = list(waic = TRUE),
                           family = "binomial",
                           control.predictor=list(link=1),
                           data = data)
  
  cat(paste0("Finished running spatial only on ", feature, " and the time is ", Sys.time(), ".\n"))
  
  #### AUTOTYP-area model ####
  
  formula = eval(substitute(this_feature ~ 
                              f(AUTOTYP_area_id_iid_model, 
                                hyper = pcprior,
                                model = "iid"),
                            list(this_feature=as.name(feature))))
  
  AUTOTYP_area_model <- INLA::inla(formula = formula,
                                  control.compute = list(waic = TRUE),
                                  family = "binomial",
                                  control.predictor=list(link=1),
                                  data = data)
  
 
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
                          control.compute = list(waic = TRUE),
                          family = "binomial",
                          control.predictor=list(link=1),
                          data = data)
  
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
                           control.compute = list(waic = TRUE),
                           family = "binomial",
                           control.predictor=list(link=1),
                           data = data)
  
  cat(paste0("Finished running trial on ", feature, " and the time is ", Sys.time(), ".\n"))
  
  #strip interesting information from the INLA objects
  phylogeny_only = strip_inla(phylo_only_model)
  cat(paste0("Finished running strip inla on phylo only  ", feature, " and the time is ", Sys.time(), ".\n"))
  
  spatial_only = strip_inla(spatialonly_model)
  
  cat(paste0("Finished running strip inla on spatial only of ", feature, " and the time is ", Sys.time(), ".\n"))

  AUTOTYP_area = strip_inla(  AUTOTYP_area_model)
  
  cat(paste0("Finished running strip inla on AUTOTYP_area_model of ", feature, " and the time is ", Sys.time(), ".\n"))
  
  dual_model = strip_inla(dual_model)
  
  cat(paste0("Finished running strip inla on dual  of ", feature, " and the time is ", Sys.time(), ".\n"))
  
  trial_model = strip_inla(trial_model)
  
  cat(paste0("Finished running strip inla on trial model of  ", feature, " and the time is ", Sys.time(), ".\n"))
  
  
  ## For each model we strip out the information we need to calculate heritability scores
  ## to save on storage space. 
  model_outputs = list(
    phylogeny_only,
    spatial_only,
    AUTOTYP_area,
    dual_model, 
    trial_model)
  
  #### Save output ####
  saved_file =   paste0(OUTPUTDIR, feature, ".qs")
  
  qs::qsave(model_outputs, 
            file = saved_file)
  
  if(beep == 1){
  beep()
    }
  
  index <- index + 1
  
  cat(paste0("I've finished ", feature, " and the time is ", Sys.time(), ".\n That means I'm ", round((index/length(range)), 2)*100, "% done.\n"))
  rm(model_outputs, phylo_only_model, spatialonly_model, dual_model)
}

sink(file = NULL)