#This is a script for running binomial INLA over 113 binary Grambank features, with phylo and spatial effects.

#set this as 1 if you're just running this script on 50 lgs over 3 features to debug. Otherwise set to 0.
source("fun_def_h_load.R")
h_load(pkg = c("ape", "tidyverse", "assertthat"))

# load variational covariance matrix function taken from geoR::varcov_spatial
source('spatiophylogenetic_modelling/analysis/varcov_spatial.R')

#make output dirs
if (!dir.exists("output/spatiophylogenetic_modelling/")) {
  dir.create("output/spatiophylogenetic_modelling/")
}
 OUTPUTDIR_vcv  <- file.path("output/spatiophylogenetic_modelling/")

#### Functions ####

cov2precision = function(spatial_covar_mat){
  spatial_covar_mat = spatial_covar_mat / exp(determinant(spatial_covar_mat)$modulus[1] /
                                                nrow(spatial_covar_mat))
  spatial_prec_mat = solve(spatial_covar_mat)
  spatial_prec_mat
}

#### Main Analyses ####

#### Inputs ####
 #reading in GB
 GB_filename <- file.path("output", "GB_wide", "GB_cropped_for_missing.tsv")
 if (!file.exists(GB_filename)) { 
   source("make_wide.R")
   source("make_wide_binarized.R")
   source("impute_missing_values.R")}		
 GB <- read.delim(GB_filename, sep = "\t")

random_cols <- sample(1:length(GB), size = 3)

if(debug_run == 1) { #if we're running to debug, pick 50 random languages and 3 random features
   GB <- GB %>% 
     sample_n(50) %>% 
     .[,c(1, random_cols)]
}

 # language metadata
if (!file.exists("output/non_GB_datasets/glottolog_AUTOTYP_areas.tsv")) { source("unusualness/processing/assigning_AUTOTYP_areas.R") }		
autotyp_area <- read.delim("output/non_GB_datasets/glottolog_AUTOTYP_areas.tsv", sep = "\t") %>%
  dplyr::select(Language_ID, AUTOTYP_area)

glottolog_df_fn = "output/non_GB_datasets/glottolog-cldf_wide_df.tsv"
if (!file.exists(glottolog_df_fn)) { source("make_glottolog-cldf_table.R") }		

languages <- read.delim(glottolog_df_fn, sep = "\t") %>%		
  dplyr::select(Language_ID, Family_ID, Name, Longitude, Latitude, Macroarea) %>% 
  distinct(Language_ID, .keep_all = T) %>% 
  inner_join(dplyr::select(GB, "Language_ID"), by = "Language_ID") %>% 
  mutate(Longitude = round(Longitude, 3)) %>% # let's cut down the precision of the lat/long to make the process go quicker. See stack exchange thread where they say "The third decimal place is worth up to 110 m: it can identify a large agricultural field or institutional campus." https://gis.stackexchange.com/questions/8650/measuring-accuracy-of-latitude-and-longitude
  mutate(Latitude = round(Latitude, 3)) %>% 
  left_join(autotyp_area, by = "Language_ID")

# trees
tree_filename = 'output/spatiophylogenetic_modelling/processed_data/EDGE_pruned_tree.tree'
if (!file.exists(tree_filename)) { source("spatiophylogenetic_modelling/processing/pruning_EDGE_tree.R") }		
phylogenetic_tree = read.tree(tree_filename)

# Subset GB and languages to overlap with tree
GB  <- GB[GB$Language_ID %in% phylogenetic_tree$tip.label,]
languages = languages[languages$Language_ID %in% GB $Language_ID,]

# prune tree to dataset
taxa = GB$Language_ID
phylogenetic_tree = keep.tip(phylogenetic_tree, tip = taxa)

#### Parameters ####
#### Spatial Jittering ####
## There are some number of languages that have identical spatial coordinates, which we cannot allow for the spatial analysis.
## I have jittered the coordinates that are identical
languages  <- languages %>% 
  mutate(longlat = paste0(Longitude, " - ", Latitude)) %>% 
  mutate(longlat = ifelse(longlat == "NA - NA", NA, longlat)) %>% 
  mutate(dup_longlat = (duplicated(longlat) + duplicated(longlat, fromLast = TRUE))) %>% 
  mutate(dup_longlat = ifelse(is.na(Longitude), NA, dup_longlat )) %>% 
  mutate(Longitude = ifelse(dup_longlat >= 1, jitter(Longitude, factor = 1), Longitude)) %>% 
  mutate(Latitude = ifelse(dup_longlat >= 1, jitter(Latitude, factor = 1), Latitude))

#### Phylogenetic covariance matrix ####
cat("Calculating the phylogenetic variance covariance matrix.\n")

phylo_covar_mat <- ape::vcv(phylogenetic_tree)
phylo_covar_mat <- phylo_covar_mat / max(phylo_covar_mat)
# The diagonal of phylo_covar_mat should inform our prior
phylo_prec_mat = cov2precision(phylo_covar_mat)

phylo_covar_mat %>% 
  as.data.frame() %>% 
  rownames_to_column("Language_ID") %>% 
  write_tsv(file = paste0(OUTPUTDIR_vcv, "/phylo_covar_mat.tsv"))


phylo_prec_mat  %>% 
  as.data.frame() %>% 
  rownames_to_column("Language_ID") %>% 
  write_tsv(file = paste0(OUTPUTDIR_vcv, "/phylo_prec_mat.tsv"))

#### Spatial covariance matrix ####

cat("Calculating the spatial variance covariance matrix.\n")
## Ensure the order of languages matches the order within the phylogeny
languages = languages[order(match(languages$Language_ID, rownames(phylo_prec_mat))),]

spatial_covar_mat = varcov.spatial(languages[,c("Longitude", "Latitude")], 
                                   cov.pars = sigma, kappa = kappa)$varcov
dimnames(spatial_covar_mat) = list(languages$Language_ID, languages$Language_ID)

spatial_prec_mat = cov2precision(spatial_covar_mat)

# Do Phylo and Spatial rownames match?
x = assert_that(all(rownames(phylo_prec_mat) == rownames(spatial_covar_mat)),
                msg = "Spatial and Phylo matrices do not align")


spatial_prec_mat %>% 
  as.data.frame() %>% 
  rownames_to_column("Language_ID") %>% 
  write_tsv(file = paste0(OUTPUTDIR_vcv, "/spatial_prec_mat.tsv"))

spatial_covar_mat %>% 
  as.data.frame() %>% 
  rownames_to_column("Language_ID") %>% 
  write_tsv(file = paste0(OUTPUTDIR_vcv, "/spatial_covar_mat.tsv"))