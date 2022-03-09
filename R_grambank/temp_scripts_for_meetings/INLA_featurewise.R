source("requirements.R")

source("spatiophylogenetic_modelling/processing/pruning_jagertree.R")

# load variational covariance matrix taken from geoR::varcov_spatial
source('spatiophylogenetic_modelling/analysis/varcov_spatial.R')

# Check INLA is installed
if (!is_installed("INLA")) { cat("INLA wasn't installed, installing now.\n") 
  source(file.path("spatiophylogenetic_modelling", "install_inla.R")) } else {
    cat("Great, INLA was already installed, loading now.\n") }
suppressPackageStartupMessages(
  library(INLA, quietly = T, warn.conflicts = F, verbose = F)
)





#Running