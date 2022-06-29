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

#visualising model fits
model_scores_df_fn <- "output/spatiophylogenetic_modelling/featurewise/model_scores.tsv"
if(!file.exists(model_scores_df_fn)){
  source("spatiophylogenetic_modelling/analysis/functions/extract_model_fit_scores.R")
}
model_scores_df <- read_tsv(model_scores_df_fn, na = "", show_col_types = F)

model_scores_df %>% 
  reshape2::melt(id.vars = "Feature_ID") %>% 
  filter(str_detect(variable, "waic")) %>% 
  group_by(Feature_ID) %>% 
  slice_min(order_by = value)  %>% 
  group_by(variable) %>% 
  summarise(n = n() )

colnames(model_scores_df) <- colnames(model_scores_df) %>% str_replace_all("_", "\n ")

model_scores_df %>%
  reshape2::melt(id.vars = "Feature\n ID") %>% 
  filter(str_detect(variable, "cpo")) %>% 
  ggplot(aes(x = variable, y = value, color = variable))+
  geom_boxplot() +
  geom_jitter() +
  theme_classic() +  
  theme(legend.position = "none", 
        axis.title.x = element_blank()) +
  ylab("CPO")

ggsave("output/spatiophylogenetic_modelling/featurewise/boxplot_waic_scorse_featurewise.png")


model_scores_df %>%
  reshape2::melt(id.vars = "Feature\n ID") %>% 
  filter(str_detect(variable, "waic")) %>% 
  group_by(`Feature\n ID`) %>% 
  slice_min(value) %>% 
  group_by(variable) %>% 
  summarise(n = n())

png("output/spatiophylogenetic_modelling/featurewise/SLOM_phylo_only_model_fits.png")
model_scores_df %>% 
  dplyr::select("phylogeny\n only\n waic", "phylogeny\n only\n cpo","phylogeny\n only\n pit",  "phylogeny\n only\n mlik\n gaussian") %>% 
  pairs.panels( 
    method = "pearson", # correlation method
    hist.col = "#a3afd1",# "#a9d1a3","",""),
    density = TRUE,  # show density plots
    ellipses = F, # show correlation ellipses
    cex.labels= 1,
    #           smoother= T,
    cor=T,
    lm=T,
    ci = T, cex.cor = 0.9,stars = T
  )
x <- dev.off()

