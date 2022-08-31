source("requirements.R")		

pPCA_importance <- read_tsv("output/PCA/ppca_importance.tsv", show_col_types = F) %>% 
  distinct(PC, `Proportion of Variance`) %>% 
  mutate(kind = "pPCA") %>% 
  .[1:19,]

pPCA_importance$PC <- fct_reorder(pPCA_importance$PC, 1-pPCA_importance$`Proportion of Variance`)

PCA_importance <- read_tsv("output/PCA/PCA_rotations.tsv", show_col_types = F) %>% 
  distinct(PC, `Proportion of Variance`) %>% 
  mutate(kind = "PCA") %>% 
  .[1:19,]

PCA_importance$PC <- fct_reorder(PCA_importance$PC, 1-PCA_importance$`Proportion of Variance`)

joined_df <- full_join(pPCA_importance, PCA_importance)

joined_df  %>% 
ggplot(aes(x = PC, y = `Proportion of Variance`, color = kind)) +
  geom_point() +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1))

ggsave(filename = "output/PCA/PCA_pPCA_variance.png")

  