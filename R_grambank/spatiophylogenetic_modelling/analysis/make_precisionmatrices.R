### Precison matrix set-up for INLA analysis.

source("fun_def_h_load.R")
h_load(pkg = c("ape", "adephylo", "MCMCglmm", "assertthat", "stringr", "tidyverse"))

#prep necessary variables and functions
source('spatiophylogenetic_modelling/analysis/functions/varcov_spatial.R')
source("spatiophylogenetic_modelling/analysis/INLA_parameters.R")

#### Phylogenetic Precision ####
tree_fn <- "output/spatiophylogenetic_modelling/processed_data/EDGE_pruned_tree.tree"
if(!(file.exists(tree_fn))){
  source("spatiophylogenetic_modelling/processing/pruning_EDGE_tree.R")}
tree = read.tree(tree_fn)

#double check that subset to lgs in GB cropped dataset
tree <- ape::keep.tip(tree, lgs_in_analysis$Language_ID)

#rescale the branches
tree$edge.length = tree$edge.length / 1000

#check that all the tip labels in the tree match the GB and vice versa
x <- assertthat::assert_that(all(tree$tip.label %in% lgs_in_analysis$Language_ID))

# We want the phylogenetic matrix to have a variance of about 1 to make it comparable to
# the spatial matrix and to compare across trees.
# The easiest way to scale the phylogenetic covariance matrix is to scale the branch-lengths so the root to tip
# distance is 1 (for an ultrametric tree), or the maximum root to tip distance is one (for non-ultrametric)

# But I would recommend actually scaling the covariance matrix after it has been calculated by dividing its
# elements by its "typical" variance. So this is related to the idea of using the diagonals as I mentioned before,
# but applying it to scaling the covariance instead of to the prior.
# That way, you can base the prior on a typical variance of one, and use this scale for all variance factors.
# To do that you can apply  the same transformation on your spatial covariance matrix, to make them comparable,
# and use the same prior. The typical variance is calculated as exp(mean(log(diag(covar)))), where covar is your
# covariance matrix

## Calculate precision using typical variance rescaling
# Here, we calculate the precison matrix, including all nodes and tips
# By including the nodes, we create a sparse matrix, which results in significant
# time improvements within INLA. Note we don't want to scale the phylogeny
# because we are doing that ourselves in a moment
phy_inv_nodes = MCMCglmm::inverseA(tree,
                                   nodes = "ALL",
                                   scale = FALSE)$Ainv

# Next, we invert the precison matrix - creating the covariance matrix
# and standardize by the typical variance, to ensure variance is scaled to 1
phy_covar_nodes = solve(phy_inv_nodes)
typical_phylogenetic_variance = exp(mean(log(diag(phy_covar_nodes))))
phy_cov_std = phy_covar_nodes / typical_phylogenetic_variance

# show that variance is scaled to approximately 1
# this shows an average and upper bound of 1, but
# But a minimum value of ~0. This is probably the impact of internal nodes

mean <- summary(diag(phy_cov_std))[["Mean"]]
x <- all.equal(mean, 1, tolerance = 0.05)
cat(paste0("The mean of the phylo covariance matrix is ", round(mean, 4), ".\n"))

min <- summary(diag(phy_cov_std))[["Min."]]
x <- all.equal(min, 0, tolerance = 0.05)
cat(paste0("The min of the phylo covariance matrix is ", round(min, 4), ".\n"))

#set dim names. This includes internal nodes and tips. Internal nodes will have names like "Node54", tips will have glottocodes.
dimnames(phy_cov_std) = dimnames(phy_inv_nodes)

# If we look at the scale of variance for nodes and not nodes (i.e. tips)
# we see that the variance for tips is always approximately 1 (which is what we want to see)
node_idx = str_detect(rownames(phy_cov_std), "Node")
# Nodes
summary(diag(phy_cov_std)[node_idx])
# Tips
mean <- summary(diag(phy_cov_std)[!node_idx])[["Mean"]]
cat(paste0("The mean of the phylo covariance matrix for tips is ", round(mean, 4), ".\n"))

# Convert the typical variance standardized covariance matrix back to a precison matrix
# we will do this by transforming the branchlengths of the tree and then using inverseA
# again, otherwise numerically precision issues will mean our precision matrix won't
# be properly sparse

tree_scaled <- tree
tree_scaled$edge.length <- tree_scaled$edge.length / typical_phylogenetic_variance
phy_prec_mat_new <- MCMCglmm::inverseA(tree_scaled,
                                       nodes = "ALL",
                                       scale = FALSE)$Ainv

# double check the matrix is more or less the same
phy_prec_mat = solve(phy_cov_std)
dimnames(phy_prec_mat) = dimnames(phy_inv_nodes)

## matrices are identical up to tolerance
all.equal(as.matrix(phy_prec_mat), as.matrix(phy_prec_mat_new))

## But...
## technically not that sparse:
sum(phy_prec_mat != 0)
## dgeMatrix is a dense matrix format
class(phy_prec_mat)
## this is properly sparse
sum(phy_prec_mat_new != 0)
## Now we have a proper dgCMatrix, a sparse format
# class(phy_prec_mat_new)

## using this new phylo precision matrix should perform much better
phy_prec_mat <- phy_prec_mat_new

#### Spatial Precison ####
# Get the longitude and latitude data from the simulated datasets
locations_df = read.delim('output/non_GB_datasets/glottolog-cldf_wide_df.tsv', sep = "\t") %>%
  inner_join(lgs_in_analysis, by = "Language_ID") #subset to matches in tree and in cropped in GB.

#check that the locations df and tree tips match
x <- assertthat::assert_that(all(locations_df$Language_ID %in% tree$tip.label),  msg = "OH NO THE TREE AND LOCATIONS ARE DIFFERENT")

spatial_covar_mat = varcov.spatial(locations_df[,c("Longitude", "Latitude")],
                                   cov.pars = sigma,
                                   kappa = kappa)$varcov

## Repeat the typical variance standardisation from above
typical_variance_spatial = exp(mean(log(diag(spatial_covar_mat))))
spatial_cov_std = spatial_covar_mat / typical_variance_spatial
spatial_prec_mat = solve(spatial_cov_std)
dimnames(spatial_prec_mat) = list(locations_df$Language_ID, locations_df$Language_ID)

## Save the precision matrices to be used in each script
precision_matrices = list(
  phylogenetic_precision = phy_prec_mat,
  spatial_precision = spatial_prec_mat
)

saveRDS(precision_matrices, "output/spatiophylogenetic_modelling/processed_data/precision_matrices.RDS")
