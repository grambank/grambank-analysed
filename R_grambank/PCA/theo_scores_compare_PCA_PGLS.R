source("requirements.R")		

#reading in data
tree <- ape::read.tree("output/spatiophylogenetic_modelling/processed_data/EDGE_pruned_tree.tree")
tree_tip_df <- tree$tip.label %>% 
  as.data.frame() %>% 
  rename(Language_ID = ".")

theo_scores_fn <- "output/PCA/theo_scores.tsv"

if(!file.exists(theo_scores_fn)){
  source("theo_scores_compare_PCA_loadings.R")
}

theo_scores_df <- read_tsv(theo_scores_fn, show_col_types = F) %>%
  inner_join(tree_tip_df,by = "Language_ID") %>% 
  column_to_rownames("Language_ID") %>% 
  as.data.frame() 

#divide all theo scores by their standard deviation
theo_scores_df <- apply(theo_scores_df,2, function(x) x/sd(x)) %>% as.data.frame()

theo_scores_df <- theo_scores_df %>% 
  rownames_to_column("Language_ID")

#setting up the comparative data object for caper
comp_data <- caper::comparative.data(phy = tree, data = theo_scores_df, names.col = "Language_ID")

#testing combinations

theo_scores_cols <- c("word order"   ,       "Flexivity"      ,     "Gender/\nnoun class", "locus of\nmarking"  ,
                      "Fusion"        ,      "Informativity" )

PCS <- c("PC1_scaled" ,         "PC2_scaled"     ,     "PC3_scaled"         )

results_df <- matrix(nrow = 0, ncol = 4) %>% 
  as.data.frame()

colnames(results_df) <- c(
  "theo_score", 
    "PC", "coef" , 
  "p_value"  
)

results_df <- results_df %>% 
  mutate_all(as.numeric) %>% 
  mutate(PC = as.character(PC))  %>% 
  mutate(theo_score = as.character(theo_score)) 

for(theo_score in theo_scores_cols){
#  theo_score <- theo_scores_cols[1]
  cat("i'm on ", theo_score, ".\n")
  output_PC1 <- eval(substitute(caper::pgls(formula = PC1 ~ this_var, data = comp_data), list(this_var=as.name(theo_score)))) %>% summary()
  output_PC2 <- eval(substitute(caper::pgls(formula = PC2 ~ this_var, data = comp_data), list(this_var=as.name(theo_score))))  %>% summary()
  output_PC3 <- eval(substitute(caper::pgls(formula = PC3 ~ this_var, data = comp_data), list(this_var=as.name(theo_score))))  %>% summary()

PC1_df <- data.frame(PC = "PC1",
  theo_score = theo_score,
  coef = output_PC1$coefficients[2,1] ,
  p_value = output_PC1$coefficients[2,4] )

PC2_df <- data.frame(PC = "PC2",
                     theo_score = theo_score,
  coef = output_PC2$coefficients[2,1] ,
  p_value = output_PC2$coefficients[2,4] )
  
PC3_df <- data.frame(PC = "PC3",
                     theo_score = theo_score,
                     coef = output_PC3$coefficients[2,1] ,
                     p_value = output_PC3$coefficients[2,4] )

results_df <- results_df  %>% 
  full_join(PC1_df, by = c("theo_score", "PC", "coef", "p_value")) %>% 
  full_join(PC2_df, by = c("theo_score", "PC", "coef", "p_value")) %>% 
  full_join(PC3_df, by = c("theo_score", "PC", "coef", "p_value"))
}

results_df <- results_df %>% 
  mutate(sig = ifelse(p_value < 0.05, "sig", "non_sig")) 

results_df  %>% 
  mutate(p_value = round(p_value, 5)) %>% 
  mutate(coef = round(coef, 5)) %>% 
  mutate(theo_score = str_replace_all(theo_score, "\n", "")) %>%
  arrange(PC) %>% 
    dplyr::select(PC, `Theoretical score` = theo_score, coef, `p-value (of t)` =p_value) %>% 
  write_tsv("output/PCA/PGLS_theo_score_correlations.tsv")
