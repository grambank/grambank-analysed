source('requirements.R')

source("spatiophylogenetic_modelling/analysis/INLA_parameters.R")
pcprior <- prior_ten_percent

OUTPUTDIR <- "output/spatiophylogenetic_modelling/INLA_spec_results_summaries/"
if(!dir.exists(OUTPUTDIR)){
  dir.create(OUTPUTDIR)
}

beep = 0

#description of features
GB_id_desc <- readr::read_tsv("output/GB_wide/parameters_binary.tsv", show_col_types = F) %>%
  dplyr::select(Feature_ID = ID, Grambank_ID_desc, Name)

source("spatiophylogenetic_modelling/analysis/INLA_parameters.R")

cat("\n###\nLoading covariance matrices...\n")

precision_matrices_fn <- "output/spatiophylogenetic_modelling/processed_data/precision_matrices_kappa_2_sigma_1.15.RDS"
if(!(file.exists(precision_matrices_fn))){
  source("spatiophylogenetic_modelling/analysis/make_precisionmatrices.R")}

precision_matrices = readRDS(precision_matrices_fn)
phylo_prec_mat = precision_matrices$phylogenetic_precision
spatial_prec_mat = precision_matrices$spatial_precision

color_vector <-c("#593d9cff", "#f68f46ff", "#d6d4d4")

#### Format Posterior Data ####
## Read in model posteriors
model_output_files = list.files(
  #path = "output/spatiophylogenetic_modelling/featurewise/",
  path = "/Users/skirgard/Nextcloud/Git_output/grambank-analysed/output_old/spatiophylogenetic_modelling/featurewise/",
  pattern = ".*kappa_2_sigma_1.15_pcprior0.1.*.qs",
  full.names = TRUE)

model_output = lapply(model_output_files, qread)

names(model_output) = basename(model_output_files) %>%
  tools::file_path_sans_ext(.)

## Extract the posterior distrubtions for each dual_process model from the hypersample
dual_posterior = lapply(model_output, function(m) {
  dd = m[[1]]$hyper_sample
  binomial_error = pi^2 / 3
  # Calculate h^2
  posterior = (1 / dd) / (rowSums(1 / dd) + 1 + binomial_error)
  
  posterior
})
# Convert this to a dataframe
dual_posterior = map_df(dual_posterior, ~as.data.frame(.x), .id="id")

dual_hyperpar = lapply(model_output, function(m) {
  dd = apply(m[[1]]$hyper_sample, 2, function(x) median(log(x)))
  dd
})

dual_hyperpar <- do.call(rbind, dual_hyperpar) %>%
  as.data.frame()
dual_hyperpar$Feature_ID_model <- rownames(dual_hyperpar)

colnames(dual_posterior) = c(
  "Feature_ID_model",
  "spatial",
  "phylogenetic"
)

dual_posterior$Feature_ID <- str_extract(dual_posterior$Feature_ID_model, "GB[0-9]{3}a?b?")
dual_hyperpar$Feature_ID <- str_extract(dual_hyperpar$Feature_ID_model, "GB[0-9]{3}a?b?")

# join feature metadata to posterior
dual_posterior = left_join(dual_posterior, GB_id_desc, by ="Feature_ID") 

# Summarise the posterior distributions
dual_summary = dual_posterior %>%
  group_by(Feature_ID) %>%
  summarise(mean_phylogenetic = mean(phylogenetic),
            mean_spatial = mean(spatial),
            error_phylogenetic = sd(phylogenetic),
            error_spatial = sd(spatial),
            #            domain = first(Main_domain),
            cor = list(cor(cbind(phylogenetic, spatial))),
            median_phylogenetic = median(phylogenetic),
            median_spatial = median(spatial)) %>%
  arrange(desc(mean_phylogenetic)) %>%
  left_join(GB_id_desc, by = "Feature_ID")

glottolog_df <- read_tsv("output/non_GB_datasets/glottolog-cldf_wide_df.tsv", show_col_types = F) %>%
  mutate(Family_ID = ifelse(is.na(Family_ID), "Isolate", Family_ID))

#reading in GB
GB_fn <- "output/GB_wide/GB_cropped_for_missing.tsv"
if(!file.exists(GB_fn)){
  cat(paste0("Making GB data wide, binarising, cropping etc..\n"))
  source("make_wide.R")
  source("make_wide_binarized.R")
  source("impute_missing_values.R")
}
GB_df <- readr::read_tsv(file =   GB_fn,show_col_types = F) %>%
  inner_join(lgs_in_analysis, by = "Language_ID")

#reading in tree
tree_fn <- "output/spatiophylogenetic_modelling/processed_data/EDGE_pruned_tree.tree"
if(!(file.exists(tree_fn))){
  source("spatiophylogenetic_modelling/processing/pruning_EDGE_tree.R")}
tree = read.tree(tree_fn)

#double check that subset to lgs in GB cropped dataset
tree <- ape::keep.tip(tree, lgs_in_analysis$Language_ID)

