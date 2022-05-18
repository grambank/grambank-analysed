#This is a script for running binomial INLA over 113 binary Grambank features, with phylo and spatial effects.

#set this as 1 if you're just running this script on 50 lgs over 3 features to debug. Otherwise set to 0.
debug_run = 1

#should the scripts output an rds file for each output from INLA::inla() or not? Set to 0 if not, 1 if yes.
save_RDS_featurewise <- 1

kappa = 2
sigma = c(1, 1.15)

join_columns = c("2.5%", "50%", "97.5%", "Feature_ID", 
                 "effect", "waic", "model")

source("fun_def_h_load.R")

h_load(pkg = c("tidyverse", "reshape2","ape", "rlang", "assertthat" , "qs"))

#If the tree hasn't been prune yet - prune the tree :)
if (!file.exists("output/spatiophylogenetic_modelling/processed_data/jaeger_pruned.tree")) {
  source("spatiophylogenetic_modelling/processing/pruning_jagertree.R") 
}		

# load variational covariance matrix function taken from geoR::varcov_spatial
source('spatiophylogenetic_modelling/analysis/varcov_spatial.R')

# Check that INLA is installed

if (!("INLA" %in% rownames(installed.packages()))) { 
  cat("INLA wasn't installed, installing now.\n") 
  source(file.path("spatiophylogenetic_modelling", "install_inla.R")) 
}

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

if (!dir.exists(file.path(  OUTPUTDIR , "phylo_only/"))) {
  dir.create(file.path(  OUTPUTDIR , "phylo_only/"))
  dir.create(file.path(  OUTPUTDIR , "spatial_only/"))
  dir.create(file.path(  OUTPUTDIR , "autotyp_area_only/"))
  dir.create(file.path(  OUTPUTDIR , "dual_process_rdata/"))
  dir.create(file.path(  OUTPUTDIR , "trial_process_rdata/"))
  }		



cat("Starting INLA runs at", as.character(Sys.time()), ".\n")

#### Functions ####

cov2precision = function(spatial_covar_mat){
  spatial_covar_mat = spatial_covar_mat / exp(determinant(spatial_covar_mat)$modulus[1] /
                                                nrow(spatial_covar_mat))
  spatial_prec_mat = solve(spatial_covar_mat)
  spatial_prec_mat
}

#### Main Analyses ####

cat("#### Building Jaeger tree models ####\n")

#reading in GB
GB_imputed_filename <- file.path("output", "GB_wide", "GB_wide_imputed_binarized.tsv")
if (!file.exists(GB_imputed_filename)) { 
  source("make_wide.R")
  source("make_wide_binarized.R")
  source("impute_missing_values.R")}		
GB_imputed <- read.delim(GB_imputed_filename, sep = "\t")

#subset GB to test code for debugging
if(debug_run == 1){
GB_imputed <- GB_imputed[1:150,]
}

#### Inputs ####
# language metadata
if (!file.exists("output/non_GB_datasets/glottolog_AUTOTYP_areas.tsv")) { source("unusualness/processing/assigning_AUTOTYP_areas.R") }		
autotyp_area <- read.delim("output/non_GB_datasets/glottolog_AUTOTYP_areas.tsv", sep = "\t") %>%
  dplyr::select(Language_ID, AUTOTYP_area)

glottolog_df_fn = "output/non_GB_datasets/glottolog-cldf_wide_df.tsv"
if (!file.exists(glottolog_df_fn)) { source("make_glottolog-cldf_table.R") }		

languages <- read.delim(glottolog_df_fn, sep = "\t") %>%		
  dplyr::select(Language_ID, Family_ID, Name, Longitude, Latitude, Macroarea) %>% 
  distinct(Language_ID, .keep_all = T) %>% 
  inner_join(dplyr::select(GB_imputed, "Language_ID"), by = "Language_ID") %>% 
  mutate(Longitude = round(Longitude, 3)) %>% # let's cut down the precision of the lat/long to make the process go quicker. See stack exchange thread where they say "The third decimal place is worth up to 110 m: it can identify a large agricultural field or institutional campus." https://gis.stackexchange.com/questions/8650/measuring-accuracy-of-latitude-and-longitude
  mutate(Latitude = round(Latitude, 3)) %>% 
  left_join(autotyp_area, by = "Language_ID")

# trees
tree_filename = 'output/spatiophylogenetic_modelling/processed_data/jaeger_pruned.tree'
if (!file.exists(tree_filename)) { source("spatiophylogenetic_modelling/processing/pruning_jagertree.R") }		
phylogenetic_tree = read.tree(tree_filename)

# Subset GB and languages to Jaeger set
GB_imputed  <- GB_imputed[GB_imputed$Language_ID %in% phylogenetic_tree$tip.label,]
languages = languages[languages$Language_ID %in% GB_imputed $Language_ID,]

# prune tree to dataset
taxa = GB_imputed$Language_ID
phylogenetic_tree = keep.tip(phylogenetic_tree, tip = taxa)

#### Parameters ####
#### Spatial Jittering ####
## There are some number of languages that have identical spatial coordinates, which we cannot allow for the spatial analysis.
## I have jittered the coordinates that are identical
## Jittering moves locations randomly by about 0.3 and 1 degree in Longitude & Latitude
duplicate_coords = languages[duplicated(languages[,c("Longitude", "Latitude")]) | 
                               duplicated(languages[,c("Longitude", "Latitude")], 
                                          fromLast = TRUE),"Language_ID"]
duplicate_rowid = languages$Language_ID %in% duplicate_coords
languages$Latitude[duplicate_rowid] = jitter(languages$Latitude[duplicate_rowid], 
                                             factor = 1)
languages$Longitude[duplicate_rowid] = jitter(languages$Longitude[duplicate_rowid], 
                                              factor = 1)

#### Phylogenetic covariance matrix ####
cat("Calculating the phylogenetic variance covariance matrix.\n")

phylo_covar_mat <- ape::vcv(phylogenetic_tree)
phylo_covar_mat <- phylo_covar_mat / max(phylo_covar_mat)
# The diagonal of phylo_covar_mat should inform our prior
phylo_prec_mat = cov2precision(phylo_covar_mat)

# Phylogenetic matrix is right dims #comment out if debugging swiftly
if(debug_run != 1){
x <- assert_that(all(dim(phylo_prec_mat) == c(n_overlap_imputed_and_jaeger_tree,
                n_overlap_imputed_and_jaeger_tree)), 
                msg = "The phylogeny has changed and will not match the data")
}
#### Spatial covariance matrix ####

cat("Calculating the spatial variance covariance matrix.\n")
## Ensure the order of languages matches the order within the phylogeny
languages = languages[order(match(languages$Language_ID, rownames(phylo_prec_mat))),]

spatial_covar_mat = varcov.spatial(languages[,c("Longitude", "Latitude")], 
                                   cov.pars = sigma, kappa = kappa)$varcov
dimnames(spatial_covar_mat) = list(languages$Language_ID, languages$Language_ID)

spatial_prec_mat = cov2precision(spatial_covar_mat)

# Do Phylo and Spatial rownames match?
if(debug_run != 1){
x = assert_that(all(rownames(phylo_prec_mat) == rownames(spatial_covar_mat)),
                msg = "Spatial and Phylo matrices do not align")
}

#### Set up model priors ####

## Taken from Dinnage et al. (2020):
### As in the manuscript, we use a “PC” prior, which stand for “Penalizing Complexity”.
### This is a standard prior developed by the developers of INLA, which is “weakly informative”.
### It places slightly more of the prior probability density on values close to zero,
### but has a “long tail”, which allows the data to push the parameter away from zero if
### there is good evidence (i.e. the likelihood of the data is higher).
### The PC prior has two parameters, p1 and p2: p2 is the proportion of the prior probability
### density that falls above values greater than p1.
###
### Given we are modelling (somewhat) Gaussian data that we have standardised to a variance of 1,
### we would not expect any random factor to have a variance greater than one.
### So we will set our prior to only have about 10% of its prior probability density
### above 1. We will guestimate the total variance possibly explained by the
### phylogenetic effect based on it’s diagonal entries, which are all equal
### to 2.7373648. So to get to 1, the scaling factor would have to be about 0.36

pcprior = list(prec =list(prior="pc.prec", param = c(1, 0.1)))

## Note that sigma and kappa values are set in requirements.R

## Adding random effect ids
grambank_df = GB_imputed %>%
  left_join(tibble(Language_ID = rownames(phylo_prec_mat),
                  phy_id_generic = 1:nrow(phylo_prec_mat),
                  phy_id_iid_model = 1:nrow(phylo_prec_mat),
                  spatial_id_generic = 1:nrow(spatial_prec_mat),
                  spatial_id_iid_model = 1:nrow(spatial_prec_mat)), 
            by = "Language_ID") %>% 
  left_join(languages,  by = "Language_ID") %>% 
  rename(AUTOTYP_area_id_iid_model = AUTOTYP_area)

#################
###INLA LOOPS####
#################

#features to loop over
features <- GB_imputed %>% 
  dplyr::select(-Language_ID) %>% 
  colnames() 

#subsetting for debugging code swiftly
if(debug_run == 1) {
features <- features[3:7]
}