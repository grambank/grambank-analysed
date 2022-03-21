source("requirements.R")
p_load(beepr)

parameters <- read_csv("../../grambank_grambank/grambank/docs/feature_groupings/feature_grouping_for_analysis.csv", show_col_types = F, col_types = cols())

df_phylo_only <- readRDS("spatiophylogenetic_modelling/results/df_phylo_only.Rdata") %>%  
  left_join(parameters, by = "Feature_ID") %>% 
  mutate(model = "phylo_only")

df_phylo_only$Feature_ID <- fct_reorder(df_phylo_only$Feature_ID, df_phylo_only$waic)

df_phylo_only %>% 
  ggplot() +
  geom_point(aes(x = Feature_ID, y = waic, col = Main_domain                 )) +
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 30), 
        axis.text.y = element_text(angle = 30),
        axis.ticks = element_blank(),
        axis.line = element_blank()) 

df_phylo_only %>% 
  ggplot() +
  geom_boxplot(aes(x = Main_domain, y = waic, col = Main_domain)) +
theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 30), 
        axis.text.y = element_text(angle = 30),
        axis.ticks = element_blank(),
        axis.line = element_blank()) 

df_spatial_only <- readRDS("spatiophylogenetic_modelling/results/df_spatial_only.Rdata") %>%
  left_join(parameters, by = "Feature_ID") %>% 
  mutate(model = "spatial_only")

df_dual <- readRDS("spatiophylogenetic_modelling/results/df_spatial_phylo.Rdata") %>%
  left_join(parameters, by = "Feature_ID") %>%
  mutate(model = "dual") %>%

joined <- full_join(df_phylo_only, df_spatial_only) %>% 
  full_join(df_dual) 

joined %>% 
  distinct(model, waic, Feature_ID) %>% 
  group_by(Feature_ID) %>% 
  slice(which.min(waic)) %>% View()




