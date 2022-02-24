source('requirements.R')

# load variaction covariance matrix taken from geoR::varcov_spatial 
source('spatiophylogenetic_modelling/analysis/varcov_spatial.R')

# Check INLA is installed
if (!is_installed("INLA")) { source(file.path("spatiophylogenetic_modelling", "install_inla.R")) } else {
  cat("Great, INLA was already installed.") }
suppressPackageStartupMessages(library(INLA, quietly = T, warn.conflicts = F, verbose = F))

cat("#### Detecting coder bias ####")

#### Inputs ####

# pca
pca_filename = 'PCA/PCA_language_values.tsv'
pca_components = read_tsv(pca_filename, col_types = cols()) %>% 
  mutate(PC1 =scale(PC1)) %>% 
  mutate(PC2 =scale(PC2)) %>% 
  mutate(PC3 =scale(PC3)) 

# language metadata
languages <- read_csv(GRAMBANK_LANGUAGES, col_types=LANGUAGES_COLSPEC) %>%		
  dplyr::select(Language_ID = Language_level_ID, Family_name, Coders, Longitude, Latitude) %>% 
  distinct(Language_ID, .keep_all = T) %>% 
  inner_join(pca_components, by = "Language_ID") %>% 
 dplyr::select(Language_ID, Coders, Longitude, Latitude, Family_name) 
  
# Check that there are 1as many languages as we expect
x <- assert_that(nrow(pca_components) == n_imputed, msg = "The total number of languages has changed. This might mean the data doesn't align with the phylogeny")

# trees
tree_filename = 'spatiophylogenetic_modelling/processed_data/jaeger_pruned.tree'
phylogenetic_tree = read.tree(tree_filename)

# Subset PCA and languages to Jaeger set
pca_components = pca_components[pca_components$Language_ID %in% phylogenetic_tree$tip.label,]
languages = languages[languages$Language_ID %in% pca_components$Language_ID,]

# prune tree to dataset
taxa = pca_components$Language_ID
phylogenetic_tree = keep.tip(phylogenetic_tree, tip = taxa)

#check that the dataframes and tree have the amount of entries we expect
x <- assert_that(nrow(pca_components) ==n_overlap_imputed_and_jaeger_tree , msg = "The number of languages in PCA have changed. It is no longer compatible with the phylogeny")
x <- assert_that(nrow(languages) == n_overlap_imputed_and_jaeger_tree , msg = "The number of languages in PCA have changed. It is no longer compatible with the phylogeny")
x <- assert_that(length(phylogenetic_tree$tip.label) == n_overlap_imputed_and_jaeger_tree , msg = "The number of tips in the phylogeny has changed. It will not be compatible with the dataset")

# All tips are in data
x <- assert_that(all(phylogenetic_tree$tip.label %in% pca_components$Language_ID), msg = "The data and phylogeny taxa do not match")

#### Spatial Jittering ####
## There are some number of languages that have identical spatial coordinates, which we cannot allow for the spatial analysis.
## I have jittered the coordinates that are identical
## Jittering moves locations randomly by about 0.3 and 1 degree in Longitude & Latitude
duplicate_coords = languages[duplicated(languages[,c("Longitude", "Latitude")]) | duplicated(languages[,c("Longitude", "Latitude")], fromLast = TRUE),"Language_ID"]
duplicate_rowid = languages$Language_ID %in% duplicate_coords
languages$Latitude[duplicate_rowid] = jitter(languages$Latitude[duplicate_rowid], factor = 1)
languages$Longitude[duplicate_rowid] = jitter(languages$Longitude[duplicate_rowid], factor = 1)

#### Phylogenetic covariance matrix ####
phylo_covar_mat <- ape::vcv(phylogenetic_tree)
phylo_covar_mat <- phylo_covar_mat / max(phylo_covar_mat)
# The diagonal of phylo_covar_mat should inform our prior
phylo_covar_mat <- phylo_covar_mat / exp(determinant(phylo_covar_mat)$modulus[1] /
                                           nrow(phylo_covar_mat))
phylo_prec_mat <-solve(phylo_covar_mat)

# Phylogenetic matrix is right dims
x <- assert_that(all(dim(phylo_prec_mat) == c(n_overlap_imputed_and_jaeger_tree , n_overlap_imputed_and_jaeger_tree )), msg = "The phylogeny has changed and will not match the data")

#### Spatial covariance matrix ####
## Ensure the order of languages matches the order within the phylogeny
languages = languages[order(match(languages$Language_ID, rownames(phylo_prec_mat))),]

spatial_covar_mat = varcov.spatial(languages[,c("Longitude", "Latitude")], cov.pars = sigma, kappa = kappa)$varcov 

dimnames(spatial_covar_mat) = list(languages$Language_ID, languages$Language_ID)

cov2precision = function(spatial_covar_mat){
  spatial_covar_mat = spatial_covar_mat / exp(determinant(spatial_covar_mat)$modulus[1] /
                                                nrow(spatial_covar_mat))
  spatial_prec_mat = solve(spatial_covar_mat)
  spatial_prec_mat
}

spatial_prec_mat = cov2precision(spatial_covar_mat)

#### Prior set-up & Random effect id ####

pcprior_phy = list(prec =list(prior="pc.prec", param = c(1, 0.1)))
pcprior_spa = list(prec =list(prior="pc.prec", param = c(1, 0.1)))

## Adding random effect ids
grambank_pca = pca_components %>%
  left_join(tibble(Language_ID = rownames(phylo_prec_mat),
                   phy_id = 1:nrow(phylo_prec_mat),
                   sp_id = 1:nrow(spatial_prec_mat)), by = "Language_ID") %>% 
  left_join(., languages, by = "Language_ID")

## Testing for coder influence would be easier if we have one coder / language
## Here I assume that the first coder had the most impact. 
grambank_pca$first_coder = str_extract(grambank_pca$Coders, "[^;]*")

#### Spatial & Phylo Model: Coder as fixed effect ####
cat("Building Spatiophylogenetic models: PC1\n")
PC1_FE = inla(PC1 ~ first_coder + 
                                f(phy_id, model = "generic0", Cmatrix = phylo_prec_mat, constr = TRUE, hyper = pcprior_phy) + 
                                f(sp_id, model = "generic0", Cmatrix = spatial_prec_mat, constr = TRUE, hyper = pcprior_spa), 
                              control.compute = list(waic=TRUE, dic = TRUE, config = TRUE),
                              data = grambank_pca)

cat("Building Spatiophylogenetic models: PC2\n")
PC2_FE = inla(PC2 ~ first_coder +
                                f(phy_id, model = "generic0", Cmatrix = phylo_prec_mat, constr = TRUE, hyper = pcprior_phy) + 
                                f(sp_id, model = "generic0", Cmatrix = spatial_prec_mat, constr = TRUE, hyper = pcprior_spa), 
                              control.compute = list(waic=TRUE, dic = TRUE),
                              data = grambank_pca)

cat("Building Spatiophylogenetic models: PC3\n")
PC3_FE = inla(PC3 ~ first_coder +
                                f(phy_id, model = "generic0", Cmatrix = phylo_prec_mat, constr = TRUE, hyper = pcprior_phy) + 
                                f(sp_id, model = "generic0", Cmatrix = spatial_prec_mat, constr = TRUE, hyper = pcprior_spa), 
                              control.compute = list(waic=TRUE, dic = TRUE),
                              data = grambank_pca)

## Coders w/ significant fixed effects:
idx = (grambank_pca$first_coder %in% c("SR", "MY", "TM", "GB", "AR", "HPG", "HW", "MJ", "RHE", "TWE", "TWI"))
table(grambank_pca$first_coder[idx])
length(table(grambank_pca$Family_name[grambank_pca$first_coder == "TWI"]))

#### Spatial & Phylo Model: Coder as random effect ####
cat("Building Spatiophylogenetic models: PC1\n")
PC1_RE = inla(PC1 ~ f(first_coder, model = "iid") + 
                                f(phy_id, model = "generic0", Cmatrix = phylo_prec_mat, constr = TRUE, hyper = pcprior_phy) + 
                                f(sp_id, model = "generic0", Cmatrix = spatial_prec_mat, constr = TRUE, hyper = pcprior_spa), 
                              control.compute = list(waic=TRUE, dic = TRUE, config = TRUE),
                              data = grambank_pca)

cat("Building Spatiophylogenetic models: PC2\n")
PC2_RE = inla(PC2 ~  f(first_coder, model = "iid") + 
                                f(phy_id, model = "generic0", Cmatrix = phylo_prec_mat, constr = TRUE, hyper = pcprior_phy) + 
                                f(sp_id, model = "generic0", Cmatrix = spatial_prec_mat, constr = TRUE, hyper = pcprior_spa), 
                              control.compute = list(waic=TRUE, dic = TRUE),
                              data = grambank_pca)

cat("Building Spatiophylogenetic models: PC3\n")
PC3_RE = inla(PC3 ~ f(first_coder, model = "iid") +
                                f(phy_id, model = "generic0", Cmatrix = phylo_prec_mat, constr = TRUE, hyper = pcprior_phy) + 
                                f(sp_id, model = "generic0", Cmatrix = spatial_prec_mat, constr = TRUE, hyper = pcprior_spa), 
                              control.compute = list(waic=TRUE, dic = TRUE),
                              data = grambank_pca)


inla.tmarginal(function(x) 1/x, 
               PC1_RE$marginals.hyperpar$`Precision for phy_id`, 
               method = "linear") %>%
  inla.qmarginal(c(0.025, 0.5, 0.975), .)

inla.tmarginal(function(x) 1/x, 
               PC1_RE$marginals.hyperpar$`Precision for sp_id`, 
               method = "linear") %>%
  inla.qmarginal(c(0.025, 0.5, 0.975), .)

## Coder effects

inla.tmarginal(function(x) 1/x, 
               PC1_RE$marginals.hyperpar$`Precision for first_coder`, 
               method = "linear") %>%
  inla.qmarginal(c(0.025, 0.5, 0.975), .)

inla.tmarginal(function(x) 1/x, 
               PC2_RE$marginals.hyperpar$`Precision for first_coder`, 
               method = "linear") %>%
  inla.qmarginal(c(0.025, 0.5, 0.975), .)

inla.tmarginal(function(x) 1/x, 
               PC3_RE$marginals.hyperpar$`Precision for first_coder`, 
               method = "linear") %>%
  inla.qmarginal(c(0.025, 0.5, 0.975), .)
