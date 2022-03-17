## Binomial INLA test
#script written by Sam Passmore

source("requirements.R")
p_load(beepr, 
  phytools, 
  geiger, 
  ape, 
  caper, 
  assertthat)

# load variational covariance matrix function taken from geoR::varcov_spatial
source('spatiophylogenetic_modelling/analysis/varcov_spatial.R')

# Check that INLA is installed
if (!is_installed("INLA")) { cat("INLA wasn't installed, installing now.\n") 
  source(file.path("spatiophylogenetic_modelling", "install_inla.R")) } else {
    cat("Great, INLA was already installed, loading now.\n") }
suppressPackageStartupMessages(
  library(INLA, quietly = T, warn.conflicts = F, verbose = F)
)

# Set seed to make sure random data is the reproduced
set.seed(115588)

## Functions
## Making a spatial precision matrix
cov2precision = function(spatial_covar_mat){
  spatial_covar_mat = spatial_covar_mat / exp(determinant(spatial_covar_mat)$modulus[1] /
                                                nrow(spatial_covar_mat))
  spatial_prec_mat = solve(spatial_covar_mat)
  spatial_prec_mat
}

## Parameters
n = 1000 # I think grambank has about 1000 languages but I run on 100 to be fast
lambda = 0.6 # setting the amount of phylogenetic signal (0 - 1)

# simulate tree with n tips
random_tree = pbtree(n = n, tip.label = 1:n)
# simulate data discrete data with lambda value
y = rTraitDisc(
  rescale(random_tree, lambda, model = "lambda")
  )

table(y)

# Put into a dataframe
model_data = data.frame(y = y == "A", # Make y numeric binary (1's and 0's)
                        phy_id = names(y),# Id for phylo effects
                        spat_id = names(y), # id for spatial effects
                        ## Add Fake geographic data (should give no effect)
                        longitude = runif(n, min = -180, max = 180),
                        latitude = runif(n, min = -90, max = 90))
model_data$y = as.numeric(model_data$y)

# Check variable has the right signal using standard tests
pagels_lambda = fitDiscrete(random_tree, y, transform = "lambda")
pagels_lambda$opt$lambda # fitDiscrete not quite accurate

phylod_results = caper::phylo.d(data = model_data, 
                                 phy = random_tree, 
                                 binvar = y,
                                 names.col = phy_id)
phylod_results$DEstimate

#### Model with INLA
phylo_covar_mat <- ape::vcv(random_tree)
phylo_covar_mat <- phylo_covar_mat / max(phylo_covar_mat)
# The diagonal of phylo_covar_mat should inform our prior
phylo_covar_mat <- phylo_covar_mat / exp(determinant(phylo_covar_mat)$modulus[1] /
                                           nrow(phylo_covar_mat))
phylo_prec_mat <-solve(phylo_covar_mat)

## Use a standard prior
pcprior_phy = list(prec = list(
  prior="pc.prec", 
  param = c(1, 0.1))
  )


### Single Process
# GLM Null
no_predictors_glm = glm(y ~ 1,
                        family = "binomial",
                        data = model_data)
summary(no_predictors_glm)
exp(coef(no_predictors_glm))
sum(model_data$y == 1) / sum(model_data$y == 0)

## INLA Null model (should be the same as GLM)
no_predictors = inla(y ~ 1,
                     family = "binomial",
                     data = model_data)
# Intercept should be the ratio of 1's and 0's
exp(no_predictors$summary.fixed$mean)
sum(model_data$y == 1) / sum(model_data$y == 0)

## Grouped model
## Adding groups does nothing because there is only 1 value per group
## But just checking that nothing changes
group_model = inla(formula = y ~ f(phy_id, model = "iid"),
                   family = "binomial",
                   data= model_data)
exp(group_model$summary.fixed$mean)
sum(model_data$y == 1) / sum(model_data$y == 0)

## Adding a phylogentic effect
lambda_only = inla(formula = y ~
                     f(as.integer(phy_id), model = "generic0",
                       Cmatrix = phylo_prec_mat,
                       constr = TRUE, hyper = pcprior_phy),
                   data = model_data,
                   family = "binomial")
exp(lambda_only$summary.fixed$mean)
sum(model_data$y == 1) / sum(model_data$y == 0)
summary(lambda_only)

## What does the 50% estimate look like from the posterior?
inla.tmarginal(function(x) 1/sqrt(x),
  lambda_only$marginals.hyperpar$`Precision for as.integer(phy_id)`,
  method = "linear") %>%
  inla.qmarginal(c(0.5), .) 

## Spatial Model - should be bad at predicting our phylogenetic model
spatial_covar_mat = 
  geoR::varcov.spatial(model_data[,c("longitude", "latitude")], 
                 cov.pars = c(1, 1.15), kappa = 2)$varcov
dimnames(spatial_covar_mat) = list(model_data$spat_id, model_data$spat_id)

# Convert to spatial precision matrix
spatial_prec_mat = cov2precision(spatial_covar_mat)

spatial_only = inla(formula = y ~
                     f(as.integer(spat_id), model = "generic0",
                       Cmatrix = spatial_prec_mat,
                       constr = TRUE, hyper = pcprior_phy),
                   data = model_data,
                   family = "binomial")

# Check 50% marginal for spatial effect
inla.tmarginal(function(x) 1/sqrt(x),
               spatial_only$marginals.hyperpar$`Precision for as.integer(spat_id)`,
               method = "linear") %>%
  inla.qmarginal(c(0.5), .) 


## Dual Model
dual_model = inla(y ~ 
                    f(as.integer(phy_id), model = "generic0",
                      Cmatrix = phylo_prec_mat,
                      constr = TRUE, hyper = pcprior_phy) + 
                    f(as.integer(spat_id), model = "generic0",
                      Cmatrix = spatial_prec_mat,
                      constr = TRUE, hyper = pcprior_phy), 
                  data = model_data,
                  family = "binomial")

## Check 2.5% 50%, 97.5% range for phylogenetic and spatial effect

# Seems to track Pagel's lambda
inla.tmarginal(function(x) 1/sqrt(x),
               dual_model$marginals.hyperpar$`Precision for as.integer(phy_id)`,
               method = "linear") %>%
  inla.qmarginal(c(0.025, 0.5, 0.975), .) 

# Has a high estimate, but wide parameters suggest the model doesn't think this
# is an important parameter.
inla.tmarginal(function(x) 1/sqrt(x),
               dual_model$marginals.hyperpar$`Precision for as.integer(spat_id)`,
               method = "linear") %>%
  inla.qmarginal(c(0.025, 0.5, 0.975), .) 

