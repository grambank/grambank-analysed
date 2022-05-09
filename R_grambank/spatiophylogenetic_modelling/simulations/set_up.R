source("requirements.R")

source("spatiophylogenetic_modelling/install_inla.R")

# load variational covariance matrix function taken from geoR::varcov_spatial
source("spatiophylogenetic_modelling/analysis/varcov_spatial.R")

## functions
cov2precision = function(spatial_covar_mat){
  spatial_covar_mat = spatial_covar_mat / 
    exp(determinant(spatial_covar_mat)$modulus[1] /
          nrow(spatial_covar_mat))
  spatial_prec_mat = solve(spatial_covar_mat)
  spatial_prec_mat
}

get_lambda_inla = function(fit, effect){
  hyper_summary = fit$summary.hyperpar
  eff_1 = hyper_summary[effect, "mean"]
  eff_2 = sum(hyper_summary[,"mean"])
  binomial_constant = pi^2/3
  (1 / eff_1) / 
    ((1/eff_1) + (1/eff_2) + (1/binomial_constant))
}

# grambank metadata

#the cldf languages table only contains meta-data on the languoids, not necessarily their language-level parents. That's why we read in glottolog in full and then match that to the list of languages that are going into the other analysis, which is the imputed set.
glottolog_df <- read.delim(file = "output/non_GB_datasets/glottolog-cldf_wide_df.tsv", sep = "\t")

GB_imputed <- read.delim(file = "output/GB_wide/GB_wide_imputed_binarized.tsv", sep = "\t")

grambank_metadata = GB_imputed %>%		
  dplyr::select(Language_ID) %>% 
  left_join(glottolog_df, by = "Language_ID") %>% 
  dplyr::select(Language_ID,
                 Name, 
                Longitude, 
                Latitude) %>% 
  mutate(Longitude = round(Longitude, 3)) %>% # let's cut down the precision of the lat/long to make the process go quicker. See stack exchange thread where they say "The third decimal place is worth up to 110 m: it can identify a large agricultural field or institutional campus." https://gis.stackexchange.com/questions/8650/measuring-accuracy-of-latitude-and-longitude
  mutate(Latitude = round(Latitude, 3))

# jager tree
tree_fn <- "output/spatiophylogenetic_modelling/processed_data/jaeger_pruned.tree"
if(!file.exists(tree_fn)) {
source("spatiophylogenetic_modelling/processing/pruning_jagertree.R")}
tree <- read.tree(tree_fn)

## Subset grambank to all of those in Jager tree
keep_languages = grambank_metadata$Language_ID %in% tree$tip.label
grambank_metadata = grambank_metadata[keep_languages,]

# Sort Grambank to match tree tips
grambank_metadata = 
  grambank_metadata[match(tree$tip.labe,
                          grambank_metadata$Language_ID),]

# Use real longitude and latitude from Grambank, but let's jitter them a bit first in case they're at the same exact location
duplicate_coords = grambank_metadata[duplicated(grambank_metadata[,c("Longitude", "Latitude")]) | duplicated(grambank_metadata[,c("Longitude", "Latitude")], fromLast = TRUE),"Language_ID"]
duplicate_rowid = grambank_metadata$Language_ID %in% duplicate_coords
grambank_metadata$Latitude[duplicate_rowid] = jitter(grambank_metadata$Latitude[duplicate_rowid], factor = 1)
grambank_metadata$Longitude[duplicate_rowid] = jitter(grambank_metadata$Longitude[duplicate_rowid], factor = 1)

longitude = grambank_metadata$Longitude
latitude = grambank_metadata$Latitude

## Make matrices
#### Phylogenetic matrix
phylo_covar_mat <- ape::vcv(tree)
phylo_prec_mat = cov2precision(phylo_covar_mat)

x = assert_that(nrow(grambank_metadata) == nrow(phylo_prec_mat))

#priors
pcprior_phy = list(prec = list(
  prior="pc.prec",
  param = c(1, 0.1)) # probability that lambda is 0.1 is 10%
)

GB_freq <- GB_imputed %>% 
  column_to_rownames("Language_ID") %>% 
  as.matrix() %>% 
  as.numeric() %>% 
  mean()

cat("The percentage of 1:s in GB is ", round(GB_freq, 2) * 100, "%.\n")

OUTPUTDIR <- "output/spatiophylogenetic_modelling/simulation/"
if(!dir.exists(OUTPUTDIR)){
  dir.create(OUTPUTDIR)
}