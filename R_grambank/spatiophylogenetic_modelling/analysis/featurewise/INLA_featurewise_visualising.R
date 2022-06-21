source("requirements.R")

parameters_binary <- read_csv("feature_grouping_for_analysis.csv", show_col_types = F)

posteriors_df <- read_tsv("output/spatiophylogenetic_modelling/featurewise/posteriors_df.tsv", show_col_types = F)

df_for_plot <- posteriors_df %>%
  dplyr::select(Feature_ID, `Precision for phylo_id_in_dual`, `Precision for spatial_id_in_dual`) 

df_for_plot %>% 
  group_by(Feature_ID) %>% 
  summarise(mean_spatial_in_dual = mean(`Precision for spatial_id_in_dual`), 
            mean_phylo_in_dual = mean(`Precision for phylo_id_in_dual`), 
  ) %>% 
  left_join(parameters_binary) %>% 
  ggplot(aes(x = mean_phylo_in_dual, y = mean_spatial_in_dual, color = Main_domain)) +
  geom_point() +
  theme_classic() +
  xlim(c(0,1))+
  ylim(c(0,1)) +
  coord_fixed()

ggsave(filename = "output/spatiophylogenetic_modelling/featurewise/dual_effect_scatterplot.png")