source("requirements.R")

h_load("adephylo")

#setting up output dir
OUTPUTDIR <- file.path("output/PCA")
if (!dir.exists(OUTPUTDIR)) { dir.create(OUTPUTDIR) }

cat("Running PHYLO Principal Component Analysis on the cropped, dialect-merged, binarised and imputed dataset.\n")

#reading in Gb
GB_imputed <- read_tsv(file.path("output", "GB_wide", "GB_wide_imputed_binarized.tsv"), col_types= cols()) %>%
  column_to_rownames("Language_ID") %>%
  as.matrix()

#reading in pruned and wrangled tree
tree <- ape::read.tree("output/spatiophylogenetic_modelling/processed_data/EDGE_pruned_tree.tree")

#subsetting gb to records also in tree
GB_imputed <- GB_imputed[tree$tip.label,]

ppca_obj <- phytools::phyl.pca(tree = tree, Y = GB_imputed,method ="BM", mode="cov")


  ppca_obj$S %>% 
  as.data.frame() %>% 
  rownames_to_column("Language_ID") %>% 
  write_tsv("output/PCA/ppca_language_values.tsv")

t(summary(ppca_obj)$importance) %>%
    as.data.frame() %>%
    rownames_to_column("PC") %>%
    mutate(`Proportion of Variance` = round(`Proportion of Variance`, digits = 2)) %>%
    dplyr::select(`Proportion of Variance`, everything()) %>% 
  write_tsv("output/PCA/ppca_importance.tsv")

  
ppca_obj$L %>% 
  as.data.frame() %>% 
  rownames_to_column("Parameter_ID") %>% 
  write_tsv("output/PCA/ppca_rotations.tsv")
