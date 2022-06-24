## This script visualises the output from a dual-process model for each Grambank feature
source('requirements.R')

## Change this vector to change the colour palette of the plots
colour_blind  <- c("#009E73", "#0072B2", "#F0E442", "#D55E00")

#### Format Posterior Data ####
## Feature Metadata
parameters_binary <- read_csv("feature_grouping_for_analysis.csv", show_col_types = F)

## Read in model posteriors
model_output_files = list.files(path = "output/spatiophylogenetic_modelling/featurewise/", 
                                pattern = "*.qs", 
                                full.names = TRUE)

model_output = lapply(model_output_files, qread)
names(model_output) = basename(model_output_files) %>% 
  tools::file_path_sans_ext(.)

## Extract the posterior distrubtions for each dual_process model from the hypersample
dual_posterior = lapply(model_output, function(m) {
  dd = m[[4]]$hyper_sample
  binomial_error = pi^2 / 3
  # Calculate h^2
  posterior = (1 / dd) / (rowSums(1 / dd) + 1 + binomial_error)
  
  posterior
})
# Convert this to a dataframe
dual_posterior = map_df(dual_posterior, ~as.data.frame(.x), .id="id")

colnames(dual_posterior) = c(
  "Feature_ID",
  "spatial",
  "phylogenetic"
)

# join feature metadata to posterior
dual_posterior = left_join(dual_posterior, parameters_binary, by ="Feature_ID")

# Summarise the posterior distributions
dual_summary = dual_posterior %>% 
  group_by(Feature_ID) %>% 
  summarise(mean_phylogenetic = mean(phylogenetic),
            mean_spatial = mean(spatial),
            error_phylogenetic = sd(phylogenetic),
            error_spatial = sd(spatial),
            domain = first(Main_domain))

#### Make Plot ####
center_plot = 
  ggplot() + 
  geom_point(data = dual_summary,
             aes(x = mean_phylogenetic, 
                 y = mean_spatial, 
                 col = domain),
             size = 0.25) + 
  geom_ellipse(data = dual_summary, # remove this geom to remove the error ellipses 
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
  coord_equal() + 
  theme(legend.position = "bottom",
        legend.title = element_blank())
  

dual_posterior2 = dual_posterior %>% 
  dplyr::filter(!Feature_ID %in% c("GB198", "GB098", "GB096", "GB111", "GB116"))

## Creating marginal density plots 

## Density plots by Domain
density_spatial_domain = 
  ggplot(data = dual_posterior,
         aes(y = spatial, 
             fill = Main_domain,
             col = Main_domain)) + 
  geom_density(alpha = 0.3, size = 0.25) + theme_void() + 
    theme(legend.position = "none")

density_phylogeny_domain = 
  ggplot(data = dual_posterior,
         aes(x = phylogenetic, 
             fill = Main_domain,
             col = Main_domain)) + 
  geom_density(alpha = 0.3, size = 0.25) + theme_void() + 
  theme(legend.position = "none")

## Desnity plots by feature
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

# Figure design grid
design = "
  14
  23
  23
"

# Build figure with domain densities
figure_domain = 
  (density_phylogeny_domain / center_plot) + density_spatial_domain + 
  plot_layout(design = design, widths = c(1, 0.2),
              heights = c(0.2, 1))

# Build figure with feature densities 
figure_feature = (density_phylogeny_feature / center_plot) + density_spatial_feature + 
  plot_layout(design = design, widths = c(1, 0.2),
              heights = c(0.2, 1))


ggsave(plot = figure_domain,
       filename = "output/spatiophylogenetic_modelling/spatiophylogenetic_figure_domaindensity.jpg",
       width = 210 / 2,
       height = 210 / 2,
       units = "mm")

ggsave(plot = figure_feature,
       filename = "output/spatiophylogenetic_modelling/spatiophylogenetic_figure_featuredensity.jpg",
       width = 210 / 2,
       height = 210 / 2,
       units = "mm")

