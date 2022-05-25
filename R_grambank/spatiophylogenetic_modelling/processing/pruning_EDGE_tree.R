source("requirements.R")

#reading in glottolog
glottolog_fn <- "output/non_GB_datasets/glottolog-cldf_wide_df.tsv"
if(!file.exists(glottolog_fn)) {
  source("make_glottolog-cldf_table.R")
}
glottolog_df <- read_tsv(glottolog_fn,col_types = cols()) %>% 
  dplyr::select(Language_ID, Language_level_ID, level, aes) %>% 
  mutate(Language_level_ID = ifelse(is.na(Language_level_ID), Language_ID, Language_level_ID)) #making language-level entities their own parent, so that we can use this column for aggregation easier.

#reading in Grambank. we're reading in the one that is cropped for missing languages and merged for dialects, since that what we're using in the analysis
GB_fn <- "output/GB_wide/GB_cropped_for_missing.tsv"
if(!file.exists(GB_fn)){
  source("make_wide.R")
  source("make_wide_binarized.R")
  source("impute_missing_values.R")
}
GB_languages <- read_tsv(GB_fn,col_types = cols()) %>% 
  dplyr::select(Language_ID) #this column is already aggregated for dialects in make_wide.R

#reading in tree
EDGE_tree <- ape::read.nexus("spatiophylogenetic_modelling/phylogenies/EDGE6635-merged-relabelled.tree")

#subsetting the tips to those in Grambank and such that there is only one per language_level_id, i.e. merging dialects and the like.
to_keep <- EDGE_tree$tip.label %>% 
  as.data.frame() %>% 
  rename(tip.label = ".") %>% 
  separate(col = tip.label , into = c("Language_ID", "Name_EDGE"), remove = F, sep = 8) %>% 
  left_join(glottolog_df, by = "Language_ID") %>% 
  inner_join(GB_languages, by = "Language_ID") %>% 
  group_by(Language_level_ID) %>% 
  sample_n(1)

#actually rpuning the tree itself
pruned_tree <- ape::keep.tip(EDGE_tree, to_keep$tip.label)

#renaming tip labels to just glottocodes
pruned_tree$tip.label <- pruned_tree$tip.label %>% 
  as.data.frame() %>% 
  rename(tip.label = ".") %>% 
  separate(col = tip.label , into = c("Language_ID", "Name_EDGE"), sep = 8) %>% 
  dplyr::select(Language_ID) %>% 
  .[,1]

pruned_tree %>% 
  write.tree("output/spatiophylogenetic_modelling/processed_data/EDGE_pruned_tree.tree")