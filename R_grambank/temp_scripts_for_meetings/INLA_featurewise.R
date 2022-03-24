source("requirements.R")
p_load(beepr)

#If the tree hasn't been prune yet - prune the tree :)
if (!file.exists("spatiophylogenetic_modelling/processed_data/jaeger_pruned.tree")) { source("spatiophylogenetic_modelling/processing/pruning_jagertree.R") }		

#make output dirs

if (!dir.exists("spatiophylogenetic_modelling/results/")) {
dir.create("spatiophylogenetic_modelling/results/")
  }
if (!dir.exists("spatiophylogenetic_modelling/results/phylo_only/")) {
dir.create("spatiophylogenetic_modelling/results/phylo_only/")
dir.create("spatiophylogenetic_modelling/results/spatial_only/")
dir.create("spatiophylogenetic_modelling/results/autotyp_area_only/")
dir.create("spatiophylogenetic_modelling/results/dual_process_rdata/")
   }		

# load variational covariance matrix function taken from geoR::varcov_spatial
source('spatiophylogenetic_modelling/analysis/varcov_spatial.R')

# Check that INLA is installed
if (!is_installed("INLA")) { cat("INLA wasn't installed, installing now.\n") 
  source(file.path("spatiophylogenetic_modelling", "install_inla.R")) } else {
    cat("Great, INLA was already installed, loading now.\n") }
suppressPackageStartupMessages(
  library(INLA, quietly = T, warn.conflicts = F, verbose = F)
)

OUTPUTDIR <- file.path("spatiophylogenetic_modelling", "results")
if (!dir.exists(OUTPUTDIR)) { dir.create(OUTPUTDIR) }		

cat("#### Building Jaeger tree models ####\n")

#reading in GB
GB_imputed_filename <- file.path("GB_wide", "GB_wide_imputed_binarized.tsv")
GB_imputed <- read_tsv(GB_imputed_filename, col_types= cols())

#### Inputs ####
# language metadata
if (!file.exists("non_GB_datasets/glottolog_AUTOTYP_areas.tsv")) { source("unusualness/processing/assigning_AUTOTYP_areas.R") }		
autotyp_area <- read.delim("non_GB_datasets/glottolog_AUTOTYP_areas.tsv", sep = "\t") %>%
  dplyr::select(Language_ID, AUTOTYP_area)

languages <- read.delim("non_GB_datasets/glottolog-cldf_wide_df.tsv") %>%		
  dplyr::select(Language_ID, Family_ID, Name, Longitude, Latitude, Macroarea) %>% 
  distinct(Language_ID, .keep_all = T) %>% 
  inner_join(dplyr::select(GB_imputed, "Language_ID"), by = "Language_ID") %>% 
  mutate(Longitude = round(Longitude, 3)) %>% # let's cut down the precision of the lat/long to make the process go quicker. See stack exchange thread where they say "The third decimal place is worth up to 110 m: it can identify a large agricultural field or institutional campus." https://gis.stackexchange.com/questions/8650/measuring-accuracy-of-latitude-and-longitude
  mutate(Latitude = round(Latitude, 3)) %>% 
  left_join(autotyp_area, by = "Language_ID")

# trees
tree_filename = 'spatiophylogenetic_modelling/processed_data/jaeger_pruned.tree'
phylogenetic_tree = read.tree(tree_filename)

# Subset PCA and languages to Jaeger set
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
duplicate_coords = languages[duplicated(languages[,c("Longitude", "Latitude")]) | duplicated(languages[,c("Longitude", "Latitude")], fromLast = TRUE),"Language_ID"]
duplicate_rowid = languages$Language_ID %in% duplicate_coords
languages$Latitude[duplicate_rowid] = jitter(languages$Latitude[duplicate_rowid], factor = 1)
languages$Longitude[duplicate_rowid] = jitter(languages$Longitude[duplicate_rowid], factor = 1)

#### Phylogenetic covariance matrix ####
cat("Calculating the phylogenetic variance covariance matrix.\n")

phylo_covar_mat <- ape::vcv(phylogenetic_tree)
phylo_covar_mat <- phylo_covar_mat / max(phylo_covar_mat)
# The diagonal of phylo_covar_mat should inform our prior
phylo_covar_mat <- phylo_covar_mat / exp(determinant(phylo_covar_mat)$modulus[1] /
                                           nrow(phylo_covar_mat))
phylo_prec_mat <-solve(phylo_covar_mat)

# Phylogenetic matrix is right dims
x <- assert_that(all(dim(phylo_prec_mat) == c(n_overlap_imputed_and_jaeger_tree, n_overlap_imputed_and_jaeger_tree)), msg = "The phylogeny has changed and will not match the data")

#### Spatial covariance matrix ####

cat("Calculating the spatial variance covariance matrix.\n")
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

# Do Phylo and Spatial rownames match?
if(all(rownames(phylo_prec_mat) == rownames(spatial_covar_mat))){
  spatial_covar_mat = spatial_covar_mat / exp(determinant(spatial_covar_mat)$modulus[1] /
                                                nrow(spatial_covar_mat))
  spatial_prec_mat = solve(spatial_covar_mat)
  #spatial_prec_mat[1:5, 1:5]
} else {
  stop("Spatial and Phylo matrices do not align")
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

prior_phylo = list(prec =list(prior="pc.prec", param = c(1, 0.1)))
prior_spatial = list(prec =list(prior="pc.prec", param = c(1, 0.1)))
prior_autotyp_area <- list(prec = list(prior = "pc.prec", param = c(1, 0.1)))

## Note that sigma and kappa are set in requirements.R

## Adding random effect ids
grambank_df = GB_imputed %>%
  left_join(tibble(Language_ID = rownames(phylo_prec_mat),
                   phy_id = 1:nrow(phylo_prec_mat),
                   sp_id = 1:nrow(spatial_prec_mat)), by = "Language_ID") %>% 
  left_join(languages,  by = "Language_ID")

#################
###INLA LOOPS####
#################
#features to loop over

features <- GB_imputed %>% 
  dplyr::select(-Language_ID) %>% 
  colnames()

cat("#### Phylogenetic only model ####\n")

#make empty df to bind to
df_phylo_only <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(df_phylo_only) <- c("2.5%","50%", "97.5%", "Feature_ID", "effect", "waic", "marginals.hyperpar.phy_id") 
df_phylo_only$`2.5%` <- as.numeric(df_phylo_only$`2.5%`)
df_phylo_only$`50%` <- as.numeric(df_phylo_only$`50%`)
df_phylo_only$`97.5%` <- as.numeric(df_phylo_only$`97.5%`)
df_phylo_only$Feature_ID <- as.character(df_phylo_only$Feature_ID)
df_phylo_only$effect <- as.character(df_phylo_only$effect)
df_phylo_only$waic <- as.numeric(df_phylo_only$waic)
df_phylo_only$marginals.hyperpar.phy_id <- as.list(df_phylo_only$marginals.hyperpar.phy_id)

index <- 0

for(feature in features){
  
  #feature <- features[1]
  
  cat(paste0("# Running the phylo-only model on feature ", feature, ". That means I'm ", round(index/length(features) * 100, 2), "% done.\n"))
  index <- index + 1 
  
  output <- eval(substitute(inla(formula = this_feature ~
                                   f(phy_id, model = "generic0", Cmatrix = phylo_prec_mat,
                                     constr = TRUE, hyper = prior_phylo),
                                 control.compute = list(waic=TRUE, dic = TRUE, mlik = FALSE, config = TRUE),
                                 control.predictor = list(compute = TRUE),
                                 data = grambank_df,family = "binomial"),
                            list(this_feature=as.name(feature))))
  
  output %>% 
    saveRDS(file = paste0("spatiophylogenetic_modelling/results/phylo_only/phylo_only_", feature, ".rdata"))
    
  phylo_effect = inla.tmarginal(function(x) 1/sqrt(x),
                                output$marginals.hyperpar$`Precision for phy_id`,
                                method = "linear") %>%
    inla.qmarginal(c(0.025, 0.5, 0.975), .)
  
  df <- phylo_effect %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
    mutate(Feature_ID = feature) %>% 
    mutate(effect = "phylo_only") %>% 
    mutate(waic = output$waic$waic)  %>% 
    mutate(marginals.hyperpar.phy_id = output$marginals.hyperpar[1])

    df_phylo_only <- df_phylo_only  %>% 
    full_join(df, by = c("2.5%", "50%", "97.5%", "Feature_ID", "effect", "waic", "marginals.hyperpar.phy_id"))
  
}
cat("All done with the phylo only model, 100% done!")

df_phylo_only %>% write_tsv("spatiophylogenetic_modelling/results/df_phylo_only.tsv")
df_phylo_only %>% saveRDS("spatiophylogenetic_modelling/results/df_phylo_only.Rdata")

###

cat("#### Spatial only Model ####\n")


#make empty df to bind to
df_spatial_only <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(df_spatial_only) <- c("2.5%","50%", "97.5%", "Feature_ID", "effect", "waic", "marginals.hyperpar.sp_id") 
df_spatial_only$`2.5%` <- as.numeric(df_spatial_only$`2.5%`)
df_spatial_only$`50%` <- as.numeric(df_spatial_only$`50%`)
df_spatial_only$`97.5%` <- as.numeric(df_spatial_only$`97.5%`)
df_spatial_only$Feature_ID <- as.character(df_spatial_only$Feature_ID)
df_spatial_only$effect <- as.character(df_spatial_only$effect)
df_spatial_only$waic <- as.numeric(df_spatial_only$waic)
df_spatial_only$marginals.hyperpar.sp_id <- as.list(df_spatial_only$marginals.hyperpar.sp_id)


index <- 0

cat("#### spatial only model ####\n")
for(feature in features){
  
  #feature <- features[1]
  
  cat(paste0("# Running the spatial-only model on feature ", feature, ". That means I'm ", round(index/length(features) * 100, 2), "% done.\n"))
  index <- index + 1 
  
  output <- eval(substitute(inla(formula = this_feature ~
                                   f(sp_id, model = "generic0", Cmatrix = spatial_prec_mat,
                                     constr = TRUE, hyper = prior_spatial),
                                 control.compute = list(waic=TRUE, dic = TRUE),
                                 data = grambank_df,family = "binomial"),
                            list(this_feature=as.name(feature))))
  
  output %>% 
    saveRDS(file = paste0("spatiophylogenetic_modelling/results/spatial_only/spatial_only_", feature, ".rdata"))
  
  spatial_effect = inla.tmarginal(function(x) 1/sqrt(x),
                                  output$marginals.hyperpar$`Precision for sp_id`,
                                  method = "linear") %>%
    inla.qmarginal(c(0.025, 0.5, 0.975), .)
  
  
  df <- spatial_effect %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
    mutate(Feature_ID = feature) %>% 
    mutate(effect = "spatial_only") %>% 
    mutate(waic = output$waic$waic) %>% 
    mutate(marginals.hyperpar.sp_id = output$marginals.hyperpar[1])
  
  df_spatial_only <- df_spatial_only  %>%
    full_join(df, by = c("2.5%", "50%", "97.5%", "Feature_ID", "effect", "waic", "marginals.hyperpar.sp_id"))
  
}

df_spatial_only %>% write_tsv("spatiophylogenetic_modelling/results/df_spatial_only.tsv")
df_spatial_only %>% saveRDS("spatiophylogenetic_modelling/results/df_spatial_only.Rdata")

cat("All done with the spatial only model, 100% done!")

###
cat("#### autotyp_area only model ####\n")

#make empty df to bind to
df_autotyp_area_only <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(df_autotyp_area_only) <- c("2.5%","50%", "97.5%", "Feature_ID", "effect", "waic", "marginals.hyperpar.autotyp_area") 
df_autotyp_area_only$`2.5%` <- as.numeric(df_autotyp_area_only$`2.5%`)
df_autotyp_area_only$`50%` <- as.numeric(df_autotyp_area_only$`50%`)
df_autotyp_area_only$`97.5%` <- as.numeric(df_autotyp_area_only$`97.5%`)
df_autotyp_area_only$Feature_ID <- as.character(df_autotyp_area_only$Feature_ID)
df_autotyp_area_only$effect <- as.character(df_autotyp_area_only$effect)
df_autotyp_area_only$waic <- as.numeric(df_autotyp_area_only$waic)
df_autotyp_area_only$marginals.hyperpar.autotyp_area <- as.list(df_autotyp_area_only$marginals.hyperpar.autotyp_area)

index <- 0

for(feature in features){
  
  #feature <- features[1]
  
  cat(paste0("# Running the autotyp_area-only model on feature ", feature, ". That means I'm ", round(index/length(features) * 100, 2), "% done.\n"))
  index <- index + 1 
  
  output <- eval(substitute(inla(formula = this_feature ~ f(AUTOTYP_area, hyper = prior_autotyp_area),
                                 control.compute = list(waic=TRUE, dic = TRUE, mlik = FALSE, config = TRUE),
                                 control.predictor = list(compute = TRUE),
                                 data = grambank_df,family = "binomial"),
                            list(this_feature=as.name(feature))))
  
  output %>% 
    saveRDS(file = paste0("spatiophylogenetic_modelling/results/autotyp_area_only/autotyp_area_only_", feature, ".rdata"))
  
  


      autotyp_area_effect = inla.tmarginal(function(x) 1/sqrt(x),
                                output$marginals.hyperpar$`Precision for AUTOTYP_area`,
                                method = "linear") %>%
    inla.qmarginal(c(0.025, 0.5, 0.975), .)
  
  df <- autotyp_area_effect %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
    mutate(Feature_ID = feature) %>% 
    mutate(effect = "autotyp_area_only") %>% 
    mutate(waic = output$waic$waic)  %>% 
    mutate(marginals.hyperpar.autotyp_area = output$marginals.hyperpar[1])
  
  df_autotyp_area_only <- df_autotyp_area_only  %>% 
    full_join(df, by = c("2.5%", "50%", "97.5%", "Feature_ID", "effect", "waic", "marginals.hyperpar.autotyp_area"))
  
}
cat("All done with the autotyp_area only model, 100% done!")

df_autotyp_area_only %>% write_tsv("spatiophylogenetic_modelling/results/df_autotyp_area_only.tsv")
df_autotyp_area_only %>% saveRDS("spatiophylogenetic_modelling/results/df_autotyp_area_only.Rdata")


df_autotyp_area_only %>% View()

cat("#### Spatial & Phylo Model ####\n")

index <- 0

#something was going awry with calculating the marginal effects of a specific features, GB051, so the for loop above has been split in twain: one which saves the entire output of inla() as an rdata object in a directory and one that calculates the effects and renders the same kind of df as above. This way running the for loop with inla() can still happen and we can debug the particulars after. 

for(feature in features){
  
  #feature <- features[16]
  
  cat(paste0("# Running the spatial-phylo (double-process) model on feature ", feature, ". That means I'm ", round(index/length(features) * 100, 2), "% done.\n"))
  index <- index + 1 
  
  output <- eval(substitute(inla(formula = this_feature ~
                                   f(phy_id, model = "generic0", Cmatrix = phylo_prec_mat, constr = TRUE, hyper = prior_phylo) + 
                                   f(sp_id, model = "generic0", Cmatrix = spatial_prec_mat, constr = TRUE, hyper = prior_spatial) ,
                                 control.compute = list(waic=TRUE),
                                 data = grambank_df, family = "binomial"),
                            list(this_feature=as.name(feature))))
  
  output %>% 
    saveRDS(file = paste0("spatiophylogenetic_modelling/results/dual_process_rdata/spatial_phylo_", feature, ".rdata"))
}


spatial_phylo_rdata_fns <- list.files("spatiophylogenetic_modelling/results/dual_process_rdata/", full.names = T, pattern = ".*rdata")

#second for loop for the dual process model because of previously discussing debugging workflow the for loop is split in twain.

#make empty df to bind to
df_spatial_phylo <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(df_spatial_phylo) <- c("2.5%","50%", "97.5%", "Feature_ID", "effect", "waic") 
df_spatial_phylo$`2.5%` <- as.numeric(df_spatial_phylo$`2.5%`)
df_spatial_phylo$`50%` <- as.numeric(df_spatial_phylo$`50%`)
df_spatial_phylo$Feature_ID <- as.character(df_spatial_phylo$Feature_ID)
df_spatial_phylo$effect <- as.character(df_spatial_phylo$effect)
df_spatial_phylo$waic <- as.numeric(df_spatial_phylo$waic)
df_spatial_phylo$marginals.hyperpar.sp_id <- as.list(df_spatial_phylo$marginals.hyperpar.sp_id)
df_spatial_phylo$marginals.hyperpar.phy_id <- as.list(df_spatial_phylo$marginals.hyperpar.phy_id)

for(fn in spatial_phylo_rdata_fns){

#  fn <- spatial_phylo_rdata_fns[1]
# fn <- "spatiophylogenetic_modelling/results/dual_process_rdata/spatial_phylo_GB051.rdata"

output <- readRDS(fn)

feature <- fn %>% str_extract("GB[0-9]*[a|b]?")

cat(paste0("I'm processing the inla output for feature ", feature, ".\n" ))
    
phylo_effect = inla.tmarginal(function(x) 1/sqrt(x), 
                              output$marginals.hyperpar$`Precision for phy_id`, 
                              method = "linear") %>%
  inla.qmarginal(c(0.025, 0.5, 0.975), .)

spatial_effect = inla.tmarginal(function(x) 1/sqrt(x), 
                                output$marginals.hyperpar$`Precision for sp_id`, 
                                method = "linear") %>%
  inla.qmarginal(c(0.025, 0.5, 0.975), .)

df_phylo <- phylo_effect %>% 
  as.data.frame() %>% 
  t() %>% 
  as.data.frame() %>% 
  rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
  mutate(Feature_ID = feature) %>% 
  mutate(effect = "phylo_in_double") %>% 
  mutate(waic = output$waic$waic) %>% 
  mutate(marginals.hyperpar.phy_id = output$marginals.hyperpar[1])

df_space <- spatial_effect %>% 
  as.data.frame() %>% 
  t() %>% 
  as.data.frame() %>% 
  rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
  mutate(Feature_ID = feature) %>% 
  mutate(effect = "spatial_in_double") %>% 
  mutate(waic = output$waic$waic) %>% 
  mutate(marginals.hyperpar.sp_id = output$marginals.hyperpar[2])

df_spatial_phylo <- df_spatial_phylo  %>% 
  full_join(df_space, by = c("2.5%", "50%", "97.5%", "Feature_ID", "effect", "waic", "marginals.hyperpar.sp_id")) %>% 
  full_join(df_phylo, by = c("2.5%", "50%", "97.5%", "Feature_ID", "effect", "waic", "marginals.hyperpar.phy_id"))

}

df_spatial_phylo %>% write_tsv("spatiophylogenetic_modelling/results/df_spatial_phylo.tsv")
df_spatial_phylo %>% saveRDS("spatiophylogenetic_modelling/results/df_spatial_phylo.Rdata")
