## Prior varying figure

source('requirements.R')

model_files = list.files('output/spatiophylogenetic_modelling/featurewise/', 
                         pattern = "GB[0-9]{3}_kappa_.*_pcprior0.1.*.qs", 
                         full.names = TRUE)

model_output_list = lapply(model_files, function(m){
  model = qread(m)
  dd = model[[1]]$hyper_sample
  binomial_error = pi^2 / 3
  # Calculate h^2
  (1 / dd) / (rowSums(1 / dd) + 1 + binomial_error)
})
names(model_output_list) = basename(model_files)

model_output = map_df(model_output_list, ~as.data.frame(.x), .id="id")

# parameters
model_output$feature = str_extract(model_output$id, "GB[0-9]{3}")
model_output$settings = str_extract(model_output$id, "kappa_\\d+([.,]\\d+)?_sigma_\\d+([.,]\\d+)?")

colnames(model_output) = str_replace_all(colnames(model_output), 
                                         pattern = " ", 
                                         replacement = ".")

# summary table
model_summary = 
  model_output %>% 
  group_by(feature, settings) %>% 
  summarise(Spatial_estimate = mean(Precision.for.spatial_id),
            Phylogenetic_estimate = mean(Precision.for.phylo_id),
            )

model_long = pivot_longer(model_summary, cols = c("Spatial_estimate", "Phylogenetic_estimate"))

col_vector <- c("purple4", "turquoise3")

p =   ggplot(model_long, aes(x = settings, y = value, group = feature, color = name)) + 
  geom_point() + 
  geom_line() + 
  ylim(c(0, 1)) + 
  ylab("Spatiophylogenetic parameter estimates") + 
  xlab("Matern spatial parameters") + 
  scale_x_discrete(labels = c('kappa = 2; sigma = 1.15',
                              'kappa = 2; sigma = 2',
                              'kappa = 2.5; sigma = 3')) + 
  theme_classic() +
  scale_color_manual(values = col_vector) +
  theme(legend.position = 'None',
        axis.text.x = element_text(angle = 45, hjust=1)) + 
  facet_wrap(~name)

ggsave(plot = p, 
       filename = "output/spatiophylogenetic_modelling/spatialparameter_effects.png", 
       width = 150, 
       height = 100, units = "mm")
