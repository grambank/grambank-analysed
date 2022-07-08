# This script tests the following:
# Pagel's lambda
# Phylo only
# Spatial only
# Dual model

source("requirements.R")

OUTPUTDIR <- "output/spatiophylogenetic_modelling/simulated_output/"
if(!dir.exists(OUTPUTDIR)){
  dir.create(OUTPUTDIR, showWarnings = FALSE)
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

#set up files to go through
#checking which ones have an qs file already and excluding them
fns <- list.files("output/spatiophylogenetic_modelling/simulated_data/", full.names = T)
fns_base <- fns %>% basename()

done_fns <- list.files("output/spatiophylogenetic_modelling/simulated_output/", pattern = "*.qs") %>% 
  as_tibble() %>% 
  rename(fns_base = value) %>% 
  mutate(fns_base = str_replace_all(fns_base, ".qs", ""))

fns_df <- tibble(fns = fns, 
       fns_base = fns_base) %>% 
      mutate(fns_base = str_replace_all(fns_base, ".csv", "")) %>% 
      anti_join(done_fns, by = "fns_base")

fns <- fns_df$fns

#### INLA Set up ####
## See the script make_precisionmatrices.R for details on how these were created
## Both matrices have their variance scaled to 1. 

precision_matrices_fn <- "output/spatiophylogenetic_modelling/precision_matrices.RDS"
if(!(file.exists(precision_matrices_fn))){
  source("spatiophylogenetic_modelling/analysis/simulations/make_precisionmatrices.R")}

precision_matrices = readRDS(precision_matrices_fn)
phylo_prec_mat = precision_matrices$phylogenetic_precision
spatial_prec_mat = precision_matrices$spatial_precision

#fetching pc prior variable
source("spatiophylogenetic_modelling/analysis/INLA_parameters.R")

tree_fn <- "output/spatiophylogenetic_modelling/processed_data/EDGE_pruned_tree.tree"
if(!(file.exists(tree_fn))){
  source("spatiophylogenetic_modelling/processing/pruning_EDGE_tree.R")}
tree = read.tree(tree_fn)

index <- 0

for(i in range){
  #  i <- 1
  filename <- fns[i]
  cat(paste0("I'm on ", filename, " and the time is ", Sys.time(), ".\n"))
  
  ### Data ####
  data = read.csv(filename)

## Since we are using a sparse phylogenetic matrix, we need to math taxa to the correct
## rows in the matrix
phylo_id = match(tree$tip.label, rownames(phylo_prec_mat))
data$phylo_id = phylo_id

## Other effects are in the same order they appear in the dataset. 
data$spatial_id = 1:nrow(spatial_prec_mat)
data$obs_id = 1:nrow(spatial_prec_mat)

cat(paste0("done with set-up on ", filename, " and the time is ", Sys.time(), ".\n"))


#### Phylo Only model ####

## The phylogenetic model contains a additive effects for the phylogenetic component
## and a random effect for the residual variance which we have fixed to 1.
## We additionally calculate WAIC scores for this and each model for comparative purposes. 
lambdaonly_model = inla(formula = y ~
                          f(phylo_id,
                            model = "generic0",
                            Cmatrix = phylo_prec_mat,
                            hyper = pcprior) +
                          f(obs_id, model = "iid", hyper = obs_hyper),
                        family = "binomial",
                        control.compute = list(waic = TRUE), 
                        data = data)


cat(paste0("Finished running phylo-only on ", filename, " and the time is ", Sys.time(), ".\n"))

#### Spatial Only model ####
## The spatial model is set up the same way as the phylogenetic model
## only containing the spatial precision matrix. 
spatialonly_model = inla(formula = y ~
                           f(spatial_id,
                             model = "generic0",
                             Cmatrix = spatial_prec_mat,
                             hyper = pcprior) +
                           f(obs_id, model = "iid", hyper = obs_hyper),
                         control.compute = list(waic = TRUE),
                         family = "binomial",
                         data = data)

cat(paste0("Finished running spatial only on ", filename, " and the time is ", Sys.time(), ".\n"))

# #### Dual Model ####
# 
# # The dual model contains the phylogenetic and spatial precision models using the same
# # priors and scaled precision matrices. As well as a single effect for residuals (constrained to 1)
dual_model = inla(y ~
                    f(spatial_id,
                      model = "generic0",
                      Cmatrix = spatial_prec_mat,
                      hyper = pcprior) +
                    f(phylo_id,
                      model = "generic0",
                      Cmatrix = phylo_prec_mat,
                      hyper = pcprior)  +
                     f(obs_id, model = "iid", hyper = obs_hyper),
                  control.compute = list(waic = TRUE),
                  family = "binomial",
                  data = data)

cat(paste0("Finished running dual on ", filename, " and the time is ", Sys.time(), ".\n"))

#### Save output ####
saved_file = basename(filename) %>% 
  tools::file_path_sans_ext() %>% 
  paste0(OUTPUTDIR, ., ".qs")

phylogeny_only = strip_inla(lambdaonly_model)
cat(paste0("Finished running strip inla on phylo only  ", filename, " and the time is ", Sys.time(), ".\n"))

spatial_only = strip_inla(spatialonly_model)

cat(paste0("Finished running strip inla on spatial only  ", filename, " and the time is ", Sys.time(), ".\n"))

dual_model = strip_inla(dual_model)

cat(paste0("Finished running strip inla on dual  ", filename, " and the time is ", Sys.time(), ".\n"))

## For each model we strip out the information we need to calculate heritability scores
## to save on storage space. 
model_outputs = list(
  phylogeny_only,
  spatial_only,
  dual_model)

qs::qsave(model_outputs, 
          file = saved_file)

index <- index + 1

cat(paste0("I've finished ", filename, " and the time is ", Sys.time(), ".\n That means I'm ", round((index/length(range)), 2)*100, "% done.\n"))
rm(model_outputs, lambdaonly_model, spatialonly_model, dual_model)
}

sink(file = NULL)