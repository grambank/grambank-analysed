source('requirements.R')
#script written by Sam Passmore

# load variation covariance matrix taken from geoR::varcov_spatial 
source('spatiophylogenetic_modelling/analysis/varcov_spatial.R')

# Check INLA is installed
if (!is_installed("INLA")) { source(file.path("spatiophylogenetic_modelling", "install_inla.R")) } else {
  cat("Great, INLA was already installed.\n") }
suppressPackageStartupMessages(library(INLA, quietly = T, warn.conflicts = F, verbose = F))

cat("#### Testing Spatial parameters via model comparison ####")

#### Inputs ####
# language metadata
languages <- read_csv(GRAMBANK_LANGUAGES, col_types=LANGUAGES_COLSPEC) %>%		
  dplyr::select(Language_ID = Language_level_ID, Family_name, Name, Longitude, Latitude, Macroarea) %>% 
  distinct(Language_ID, .keep_all = T)

# pca
pca_filename = 'PCA/PCA_language_values.tsv'
pca_components = read_tsv(pca_filename, col_types = cols()) %>% 
  mutate(PC1 =scale(PC1)) %>% 
  mutate(PC2 =scale(PC2)) %>% 
  mutate(PC3 =scale(PC3)) 
  

#check that there are as many observations as we expect
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

#check that there are as many observations as we expect
x <- assert_that(nrow(pca_components) == n_overlap_imputed_and_jaeger_tree, msg = "The number of languages in PCA have changed. It is no longer compatible with the phylogeny")
x <- assert_that(nrow(languages) == n_overlap_imputed_and_jaeger_tree, msg = "The number of languages in PCA have changed. It is no longer compatible with the phylogeny")

#check that there are as many observations as we expect
x <- assert_that(length(phylogenetic_tree$tip.label) == n_overlap_imputed_and_jaeger_tree, msg = "The number of tips in the phylogeny has changed. It will not be compatible with the dataset")

# All tips are in data
x <- assert_that(all(phylogenetic_tree$tip.label %in% pca_components$Language_ID), msg = "The data and phylogeny taxa do not match")

#### Parameters ####
pca_basename = tools::file_path_sans_ext(pca_filename) %>% basename()
tree_basename = tools::file_path_sans_ext(tree_filename) %>% basename()
out_name = paste(pca_basename, tree_basename, sep = "_")

#### Spatial Jittering ####
## There are some number of languages that have identical spatial coordinates, which we cannot allow for the spatial analysis.
## I have jittered the coordinates that are identical
## Jittering moves locations randomly by about 0.3 and 1 degree in Longitude & Latitude
duplicate_coords = languages[duplicated(languages[,c("Longitude", "Latitude")]) | duplicated(languages[,c("Longitude", "Latitude")], fromLast = TRUE),"Language_ID"]
duplicate_rowid = languages$Language_ID %in% duplicate_coords
languages$Latitude[duplicate_rowid] = jitter(languages$Latitude[duplicate_rowid], factor = 1)
languages$Longitude[duplicate_rowid] = jitter(languages$Longitude[duplicate_rowid], factor = 1)

#### Spatial covariance matrix ####
## Ensure the order of languages matches the order within the phylogeny
languages = languages[order(match(languages$Language_ID, phylogenetic_tree$tip.label)),]

cov2precision = function(spatial_covar_mat){
  spatial_covar_mat = spatial_covar_mat / 
    exp(determinant(spatial_covar_mat)$modulus[1] / nrow(spatial_covar_mat))
  spatial_prec_mat = solve(spatial_covar_mat)
  spatial_prec_mat
}

cat("These are the sigma and kappa values we'll be going through.\n")
spatial_mat_list = map2(kappa_vec, sigma_vec, function(k, s){
  cat(paste0("Kappa: ", k, "; Sigma: ", s, ".\n"))
  spatial_covar_mat = varcov.spatial(languages[,c("Longitude", "Latitude")], cov.pars = c(1, s), kappa = k)$varcov
  dimnames(spatial_covar_mat) = list(languages$Language_ID, languages$Language_ID)
  cov2precision(spatial_covar_mat)
})

pcprior_spa = list(prec =list(prior="pc.prec", param = c(1, 0.1)))

## Adding random effect ids
grambank_pca = pca_components %>%
  left_join(tibble(Language_ID = languages$Language_ID,
                   sp_id = 1:nrow(spatial_mat_list[[1]])), by = "Language_ID")


cat(paste0("Building Spatial models: sigma = ", 
             sigma_vec[1], 
             "; kappa = ", kappa_vec[1],".\n"))
spatial_s.115 = inla(PC1 ~ 
                     f(sp_id, model = "generic0", Cmatrix = spatial_mat_list[[1]], 
                       constr = TRUE, hyper = pcprior_spa), 
                   control.compute = list(waic=TRUE, dic = TRUE),
                   data = grambank_pca)

cat(paste0("Building Spatial models: sigma = ", 
             sigma_vec[2], 
             "; kappa = ", kappa_vec[2],".\n"))
spatial_s.3 = inla(PC1 ~ 
                     f(sp_id, model = "generic0", Cmatrix = spatial_mat_list[[2]], 
                       constr = TRUE, hyper = pcprior_spa), 
                   control.compute = list(waic=TRUE, dic = TRUE),
                   data = grambank_pca)

cat(paste0("Building Spatial models: sigma = ", 
             sigma_vec[3], 
             "; kappa = ", kappa_vec[3],".\n"))
spatial_s.6 = inla(PC1 ~ 
                     f(sp_id, model = "generic0", Cmatrix = spatial_mat_list[[3]], 
                       constr = TRUE, hyper = pcprior_spa), 
                   control.compute = list(waic=TRUE, dic = TRUE),
                   data = grambank_pca)

cat(paste0("Building Spatial models: sigma = ", 
             sigma_vec[4], 
             "; kappa = ", kappa_vec[4] ,".\n"))

spatial_s.9  <- inla(PC1 ~ 
                        f(sp_id, model = "generic0", Cmatrix = spatial_mat_list[[4]], 
                          constr = TRUE, hyper = pcprior_spa), 
                      control.compute = list(waic=TRUE, dic = TRUE),
                      data = grambank_pca,
                      verbose = FALSE)
 
null_model = inla(PC1 ~ 1,
                  control.compute = list(waic=TRUE, dic = TRUE),
                  data = grambank_pca)

lapply(list(spatial_s.115, spatial_s.3, spatial_s.6, spatial_s.9, null_model),
       function(x) x$neffp)

lapply(list(spatial_s.115, spatial_s.3, spatial_s.6, spatial_s.9, null_model),
       function(x) x$waic$waic)
