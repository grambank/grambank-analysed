source('requirements.R')

#### Format Posterior Data ####
## Feature Metadata
feature_groupings <- read_tsv("output/GB_wide/table_theo_scores_supp.tsv", show_col_types = F)

## Identify results by spatial decay and pcprior
filename_suffix = "_kappa_2_sigma_1.15_pcprior0.1"

## Read in model posteriors
posteriors_df <- read_tsv("output/spatiophylogenetic_modelling/featurewise/posteriors_df.tsv", show_col_types = F) %>% 
  filter(str_detect(fn, filename_suffix)) %>% 
  dplyr::select(Feature_ID, spatial = `Precision for spatial_id_in_dual`, phylogenetic = `Precision for phylo_id_in_dual`)


