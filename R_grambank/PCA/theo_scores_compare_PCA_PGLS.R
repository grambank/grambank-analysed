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
  inner_join(tree_tip_df,by = "Language_ID") 

PCA_df <- read_tsv("output/PCA/PCA_language_values.tsv", show_col_types = F) %>% 
  dplyr::select(Language_ID, PC1, PC2, PC3)

theo_scores_df <- theo_scores_df %>% 
  left_join(PCA_df, by = "Language_ID") %>% 
  mutate_if(is.numeric, scale) %>% 
  as.data.frame()

#setting up the comparative data object for caper
comp_data <- caper::comparative.data(phy = tree, data = theo_scores_df, names.col = "Language_ID")

#testing combinations

theo_scores_cols <- c("Word_Order"   ,       "Flexivity"      ,     "Gender_or_Noun_Class", "Locus_of_Marking"  ,
                      "Fusion"        ,      "Informativity" )

results_df <- matrix(nrow = 0, ncol = 5) %>% 
  as.data.frame()

colnames(results_df) <- c(
  "theo_score", 
    "PC", "coef" , 
  "p_value"  , "t_value"
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
  p_value = output_PC1$coefficients[2,4] ,
  t_value =  output_PC1$coefficients[2,3])

PC2_df <- data.frame(PC = "PC2",
                     theo_score = theo_score,
  coef = output_PC2$coefficients[2,1] ,
  p_value = output_PC2$coefficients[2,4] ,
  t_value =  output_PC1$coefficients[2,3])
  
PC3_df <- data.frame(PC = "PC3",
                     theo_score = theo_score,
                     coef = output_PC3$coefficients[2,1] ,
                     p_value = output_PC3$coefficients[2,4] ,
                     t_value =  output_PC1$coefficients[2,3])

results_df <- results_df  %>% 
  full_join(PC1_df, by = c("theo_score", "PC", "coef", "p_value", "t_value")) %>% 
  full_join(PC2_df, by = c("theo_score", "PC", "coef", "p_value", "t_value")) %>% 
  full_join(PC3_df, by = c("theo_score", "PC", "coef", "p_value", "t_value"))
}

results_df <- results_df %>% 
  mutate(sig = ifelse(p_value < 0.05, "sig", "non_sig")) 

results_df  %>% 
  mutate(p_value = round(p_value, 5)) %>% 
  mutate(coef = round(coef, 5)) %>% 
  mutate(t_value = round(t_value, 5)) %>% 
    mutate(theo_score = str_replace_all(theo_score, "\n", "")) %>% 
  arrange(PC) %>% 
    dplyr::select(PC, `Theoretical score` = theo_score, coef, `t-value` = t_value, `p-value (of t)` =p_value) %>% 
  write_tsv("output/PCA/PGLS_theo_score_correlations.tsv")

