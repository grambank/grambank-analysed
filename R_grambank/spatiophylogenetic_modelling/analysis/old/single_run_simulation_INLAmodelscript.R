#!/usr/bin/env Rscript

# This script tests the following:
# Pagel's lambda
# Phylo only
# Spatial only
# Dual model

suppressPackageStartupMessages({
  library(ape)
  library(geiger)
  library(dplyr)
  source('functions/varcov_spatial.R')
  source('functions/strip_inla.R')
})

library("INLA")
#if(cluster == 1){
#INLA::inla.binary.install(os = "Ubuntu-18.0", path = "../rlibs/INLA/bin/linux/")#}
INLA::inla.setOption(inla.mode="experimental")

args = commandArgs(trailingOnly = TRUE)

filename = args[1]

### Data ####
data = read.csv(filename)
tree = read.tree("data/EDGE_pruned_tree.tree")

#### INLA Set up ####
## See the script make_precisionmatrices.R for details on how these were created
## Both matrices have thier variance scaled to 1. 
precision_matrices = readRDS("output/precision_matrices.RDS")
phylo_prec_mat = precision_matrices$phylogenetic_precision
spatial_prec_mat = precision_matrices$spatial_precision


# we use a penalising complexity prior which are particulary suited to the
# analyses of additive models. 
# We should test the sensitivity of priors on the full data model
pcprior = list(prec = list(
  prior="pc.prec",
  param = c(1, 0.1)) # This prior suggests that the probability that variance for the random effect is greater than 1 is 10%
)


# We need to fix the residual variance to one, since it is not an identifiable quantity
# within a binomial model. 
obs_hyper <- list(prec = list(initial = log(1), fixed = TRUE))

## Since we are using a sparse phylogenetic matrix, we need to math taxa to the correct
## rows in the matrix
phylo_id = match(tree$tip.label, rownames(phylo_prec_mat))
data$phylo_id = phylo_id

## Other effects are in the same order they appear in the dataset. 
data$spatial_id = 1:nrow(spatial_prec_mat)
data$obs_id = 1:nrow(spatial_prec_mat)


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
# 
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

#### Save output ####
dir.create("output/simulated_output/", showWarnings = FALSE)

saved_file = basename(filename) %>% 
  tools::file_path_sans_ext() %>% 
  paste0("output/simulated_output/", ., ".qs")

## For each model we strip out the information we need to calculate heritability scores
## to save on storage space. 
model_outputs = list(
  phylogeny_only = strip_inla(lambdaonly_model),
  spatial_only = strip_inla(spatialonly_model),
  dual_model = strip_inla(dual_model)
)

qs::qsave(model_outputs, 
          file = saved_file)