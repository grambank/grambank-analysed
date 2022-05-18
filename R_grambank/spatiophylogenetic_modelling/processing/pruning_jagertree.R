source("requirements.R")

#script written by Sam Passmore and Hedvig Skirgård

cat("Pruning the global tree from Jäger to the Grambank dataset.\n")

OUTPUT_DIR_top <- "output/spatiophylogenetic_modelling/"
if (!dir.exists(OUTPUT_DIR_top)) {
  dir.create(OUTPUT_DIR_top)
}

OUTPUT_DIR <- "output/spatiophylogenetic_modelling/processed_data/"
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR)
}

#### Inputs ####
#Glottolog-cldf table to have a match for all dialects to their language parent. Note that the particular dialects may differ from the dialects in GB which is why we cann't use the language table from the grambank-cldf relase
glottolog_df <- read_tsv("output/non_GB_datasets/glottolog-cldf_wide_df.tsv",col_types = cols()) %>% 
  dplyr::select(Language_ID, Language_level_ID, level) %>% 
  mutate(Language_level_ID = ifelse(is.na(Language_level_ID), Language_ID, Language_level_ID)) #making language-level entities their own parent, so that we can use this column for aggregation easier.

taxa_pairing <- read.csv('spatiophylogenetic_modelling/phylogenies/taxa.csv') %>% 
  dplyr::rename(Language_ID = glottocode) %>% 
  left_join(glottolog_df, by = "Language_ID")

jaeger_tree <- read.tree('spatiophylogenetic_modelling/phylogenies/world.tre')

GB_languages <- read_tsv("output/GB_wide/GB_wide_imputed_binarized.tsv",col_types = cols()) %>% 
  dplyr::select(Language_ID) #this column is already aggregated for dialects in make_wide.R

#### Tree - Data pairs ####
coverage = sum(GB_languages$Language_ID %in% taxa_pairing$Language_level_ID) / nrow(GB_languages) 
cat("Tips in the global Jaeger tree  can be matched to", round(coverage, 2) * 100, "% of the langauges in Grambank.\n")

#### Subset to Grambank langauges ####
in_tree <- jaeger_tree$tip.label %>% 
  as.data.frame() %>% 
  dplyr::rename(taxon_full = ".") %>% #making a data frame with a column representing all the tip labels, in the right order
  mutate(taxon = str_extract(taxon_full, "[^.]+$")) %>% #extracting the part of the tip label that matches to the taxa file
  left_join(taxa_pairing, by = "taxon") %>% #pairing with the taxa file
  dplyr::select(taxon_full, Language_ID = Language_level_ID) %>% #only keeping the necessary cols for pruning and matching with GB
  inner_join(GB_languages, by = "Language_ID") %>% #prune to only tips which are also in GB
  group_by(Language_ID) %>%
  sample_n(size = 1)  #if there is more than one tip with the same language ID (duplicates or dialects), randomly remove all but one

jaeger_pruned = keep.tip(jaeger_tree, in_tree$taxon_full)

#### Relabel taxa in tree ####
#taxons are matched to the glottocode of themselves or their parent which has the level "language" in glottolog

jaeger_pruned$tip.label <- jaeger_pruned$tip.label %>% 
  as.data.frame() %>% 
  dplyr::rename(taxon = ".") %>% 
  mutate(taxon = str_extract(taxon, "[^.]+$")) %>% 
  left_join(taxa_pairing, by = "taxon") %>%
  dplyr::select(Language_level_ID) %>% 
  .[,1]

write.tree(jaeger_pruned, "output/spatiophylogenetic_modelling/processed_data/jaeger_pruned.tree")

removed <- GB_languages %>% 
  anti_join(taxa_pairing %>% dplyr::select(Language_ID = Language_level_ID), by = "Language_ID")

write.csv(removed, 'output/spatiophylogenetic_modelling/processed_data/jaeger_removed.csv')

cat("Jager tree created with", length(jaeger_pruned$tip.label), "tips.\n")
cat(nrow(removed), "languages were not paired.\n")
cat("Languages that are in Grambank, but not in the Jaeger tree can be found in the file output/spatiophylogenetic_modelling/processed_data/jaeger_removed.csv.\n")
