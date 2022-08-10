## Prior varying figure

source('requirements.R')

model_files = list.files('output/spatiophylogenetic_modelling/featurewise/', 
                         pattern = "GB[0-9]{3}_kappa_.*.qs", 
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
  summarise(spatial_estimate = mean(Precision.for.spatial_id),
            spatial_sd = sd(Precision.for.spatial_id),
            phylogenetic_estimate = mean(Precision.for.phylo_id),
            phylogenetic_sd = sd(Precision.for.phylo_id))


ggplot(model_summary, aes(x = settings, y = spatial_estimate, col = feature, group = feature)) + 
  geom_point() + 
  geom_line() + 
  theme(legend.position = 'none')
