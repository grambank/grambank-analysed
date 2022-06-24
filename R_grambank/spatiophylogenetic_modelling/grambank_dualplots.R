## visualisation
source('requirements.R')

library(qs)
library(purrr)
library(dplyr)
library(readr)
library(ggplot2)

## Colour blind colour palette
colour_blind  <- c("#009E73", "#0072B2", "#F0E442", "#D55E00")

# Make data for plot
parameters_binary <- read_csv("feature_grouping_for_analysis.csv", show_col_types = F)

model_output_files = list.files(path = "output/spatiophylogenetic_modelling/featurewise/", 
                                pattern = "*.qs", 
                                full.names = TRUE)

model_output = lapply(model_output_files, qread)
names(model_output) = basename(model_output_files) %>% 
  tools::file_path_sans_ext(.)

## Extract the posterior estimates for each dual_process model
dual_posterior = lapply(model_output, function(m) {
  dd = m[[4]]$hyper_sample
  binomial_error = pi^2 / 3
  posterior = (1 / dd) / (rowSums(1 / dd) + 1 + binomial_error)
  
  posterior
})
dual_posterior = map_df(dual_posterior, ~as.data.frame(.x), .id="id")

colnames(dual_posterior) = c(
  "Feature_ID",
  "spatial",
  "phylogenetic"
)

dual_posterior = left_join(dual_posterior, parameters_binary, by ="Feature_ID")

dual_summary = dual_posterior %>% 
  group_by(Feature_ID) %>% 
  summarise(mean_phylogenetic = mean(phylogenetic),
            mean_spatial = mean(spatial),
            error_phylogenetic = sd(phylogenetic),
            error_spatial = sd(spatial),
            domain = first(Main_domain))

## Russell Gray suggestion
center_plot = 
  ggplot() + 
  geom_point(data = dual_summary,
             aes(x = mean_phylogenetic, 
                 y = mean_spatial, 
                 col = domain),
             size = 0.25) + 
  geom_ellipse(data = dual_summary, 
    aes(x0 = mean_phylogenetic, 
                   y0 = mean_spatial,
                   a = error_phylogenetic/2, 
                   b = error_spatial/2, 
                   angle = 0,
                   fill = domain),
               alpha = 0.3,
               color = NA) + 
  ylim(c(0, 1)) + xlim(c(0, 1)) + 
  theme_minimal(base_size = 6) + 
  xlab("Variance explained by Phylogeny") + 
  ylab("Variance explained by Geography") + 
  scale_colour_manual(values = colour_blind) + 
  scale_fill_manual(values = colour_blind) + 
  theme(legend.position = "bottom",
        legend.title = element_blank())
  

dual_posterior2 = dual_posterior %>% 
  dplyr::filter(!Feature_ID %in% c("GB198", "GB098", "GB096", "GB111", "GB116"))
center_plot2 = 
  ggplot() + 
  geom_point(data = dual_summary,
             aes(x = mean_phylogenetic, 
                 y = mean_spatial, 
                 col = domain),
             size = 0.25)  + 
  ylim(c(0, 1)) + xlim(c(0, 1)) + 
  theme_minimal(base_size = 6) + 
  xlab("Variance explained by Phylogeny") + 
  ylab("Variance explained by Geography") + 
  scale_colour_manual(values = colour_blind) + 
  scale_fill_manual(values = colour_blind) + 
  theme(legend.position = "bottom",
        legend.title = element_blank())

center_plot_wlabels = 
  ggplot() + 
  geom_text(data = dual_summary,
             aes(x = mean_phylogenetic, 
                 y = mean_spatial, 
                 col = domain,
                 label = Feature_ID),
             size = 2) + 
  geom_ellipse(data = dual_summary, 
               aes(x0 = mean_phylogenetic, 
                   y0 = mean_spatial,
                   a = error_phylogenetic/2, 
                   b = error_spatial/2, 
                   angle = 0,
                   fill = domain),
               alpha = 0.3,
               color = NA) + 
  ylim(c(0, 1)) + xlim(c(0, 1)) + 
  theme_minimal(base_size = 6) + 
  xlab("Variance explained by Phylogeny") + 
  ylab("Variance explained by Geography") + 
  theme(legend.position = "bottom",
        legend.title = element_blank())


density_spatial_domain = 
  ggplot(data = dual_posterior,
         aes(y = spatial, 
             fill = Main_domain,
             col = Main_domain)) + 
  geom_density(alpha = 0.3, size = 0.25) + theme_void() + 
    theme(legend.position = "none")

density_spatial_feature = 
  ggplot(data = dual_posterior2,
         aes(y = spatial, 
             group = Feature_ID,
             fill = Main_domain,
             col = Main_domain)) + 
  geom_density(alpha = 0.3, size = 0.25) + theme_void() + 
  scale_colour_manual(values = colour_blind) + 
  scale_fill_manual(values = colour_blind) + 
  theme(legend.position = "none")

density_phylogeny_domain = 
  ggplot(data = dual_posterior,
         aes(x = phylogenetic, 
             fill = Main_domain,
             col = Main_domain)) + 
  geom_density(alpha = 0.3, size = 0.25) + theme_void() + 
  theme(legend.position = "none")

density_phylogeny_feature = 
  ggplot(data = dual_posterior2,
         aes(x = phylogenetic, 
             group = Feature_ID,
             fill = Main_domain,
             col = Main_domain)) + 
  geom_density(alpha = 0.3, size = 0.25) + theme_void() + 
  scale_colour_manual(values = colour_blind) + 
  scale_fill_manual(values = colour_blind) + 
  theme(legend.position = "none")

design = "
  14
  23
  23
"

figure = 
  (density_phylogeny / center_plot) + density_spatial + 
  plot_layout(design = design, widths = c(1, 0.2),
              heights = c(0.2, 1))

figure_2 = (density_phylogeny_feature / center_plot2) + density_spatial_feature + 
  plot_layout(design = design, widths = c(1, 0.2),
              heights = c(0.2, 1))


ggsave(plot = figure,
       filename = "output/spatiophylogenetic_modelling/spatiophylogenetic_figure.jpg",
       width = 210 / 2,
       height = 180 / 2,
       units = "mm")

ggsave(plot = figure_2,
       filename = "output/spatiophylogenetic_modelling/spatiophylogenetic_figure2.jpg",
       width = 210 / 2,
       height = 180 / 2,
       units = "mm")

