source("requirements.R")
p_load(beepr)

parameters <- read_csv("feature_grouping_for_analysis.csv", show_col_types = F, col_types = cols()) 

INPUTDIR <- "output/spatiophylogenetic_modelling/results_debug/"
OUTPUTDIR <- "output/spatiophylogenetic_modelling/figures/"
if (!dir.exists(  OUTPUTDIR )) {
  dir.create(  OUTPUTDIR )
}

df_phylo_only <- readRDS(file.path(INPUTDIR, "df_phylo_only.Rdata")) 
df_spatial_only <- readRDS(file.path(INPUTDIR, "df_spatial_only.Rdata")) 
df_dual <- readRDS(file.path(INPUTDIR, "/df_spatial_phylo.Rdata"))

joined <- full_join(df_phylo_only, df_spatial_only) %>% 
  full_join(df_dual) %>% 
  left_join(parameters, by = "Feature_ID")

lowest_waic_df <- joined %>% 
  distinct(model, waic, Feature_ID) %>% 
  group_by(Feature_ID) %>% 
  slice(which.min(waic))

#shaping my df to look like the one in Sam's spmodel plots.
plot_df_hed <- joined %>% 
  dplyr::select(Feature_ID, effect, model, "2.5%", "50%", "97.5%") %>% 
  reshape2::melt(.id.vars = c(Feature_ID ,  effect, model)) 

col_vector <- c("purple4", "turquoise3")


plot_df_summ_hed_subset <- plot_df_hed %>% 
  reshape2::dcast(Feature_ID+ effect+ model ~ variable, value.var =    "value" ) 

plot <- plot_df_summ_hed_subset %>% 
  ggplot(aes(y = Feature_ID, x = `50%`, fill = effect, col = effect), stat = "identity") +
  theme_classic() +
#  scale_colour_manual(values = col_vector) + 
#  scale_fill_manual(values = col_vector) +
  geom_errorbar(width = 0.2, 
                aes(xmin = `2.5%`, 
                    xmax = `97.5%`),
                position = position_dodge(width = 0.4)) + 
  geom_point(shape=21, size=2, position = position_dodge(width = 0.4)) + 
  theme(legend.title = element_blank(), legend.position="bottom") +
  theme(axis.text.y = element_text(size = 14))

plot(plot)  

png( "output/spatiophylogenetic_modelling/figures/featurewise_dual_process_effect_debug_normal.png", width = 10, height = 25, units = "in", res = 100)
plot(plot)
x <- dev.off()

simmo_plot_draft_1 <- plot_df_summ_hed_subset %>% 
  reshape2::dcast(Feature_ID ~ effect, value.var = "Mean") %>% 
  left_join(parameters, by = "Feature_ID") %>% 
  ggplot() +
  geom_point(aes(x = phylo, y = spatial, col = Main_domain)) +
  ylim(c(0, 1)) +
  xlim(c(0, 2.7)) +
  theme_classic() +
  theme_set(theme_gray(base_size = 18))

png( "spatiophylogenetic_modelling/figures/featurewise_dual_process_effects_phylo_vs_spatial.png", width = 12, height = 12, units = "in", res = 100)
plot(simmo_plot_draft_1)
x <- dev.off()

