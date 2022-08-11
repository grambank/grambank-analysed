source("requirements.R")		

tree <- ape::read.tree("output/spatiophylogenetic_modelling/processed_data/EDGE_pruned_tree.tree")
tree_tip_df <- tree$tip.label %>% 
  as.data.frame() %>% 
  rename(Language_ID = ".")

theo_scores_fn <- "output/PCA/theo_scores.tsv"

if(!file.exists(theo_scores_fn)){
  source("theo_scores_compare_PCA_loadings.R")
}

theo_scores_df <- read_tsv(theo_scores_fn, show_col_types = F) %>% 
  inner_join(tree_tip_df) %>% 
  as.data.frame()

rownames(theo_scores_df) <- theo_scores_df$Language_ID

comp_data <- caper::comparative.data(phy = tree, data = theo_scores_df, names.col = "Language_ID")

test <- caper::pgls(PC1 ~ Fusion, data = comp_data)

