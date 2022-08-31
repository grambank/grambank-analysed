source("requirements.R")		

pPCA_importance <- read_tsv("output/PCA/ppca_importance.tsv", show_col_types = F) %>% 
  distinct(PC, `Proportion of Variance`) %>% 
  mutate(PC = paste0("p", PC))

PCA_importance <- read_tsv("output/PCA/PCA_rotations.tsv", show_col_types = F) %>% 
  distinct(PC, `Proportion of Variance`)