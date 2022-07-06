## This script visualises the output from a trial-process model for each Grambank feature
source('requirements.R')

## Change this vector to change the colour palette of the plots
col_vector <- c("#039e37", "purple4",  "#c23c3c", "turquoise3")

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

## Extract the posterior distrubtions for each trial_process model from the hypersample
trial_posterior = lapply(model_output, function(m) {
  dd = m[[5]]$hyper_sample
  binomial_error = pi^2 / 3
  # Calculate h^2
  posterior = (1 / dd) / (rowSums(1 / dd) + 1 + binomial_error)
  
  posterior
})
# Convert this to a dataframe
trial_posterior = map_df(trial_posterior, ~as.data.frame(.x), .id="id")

colnames(trial_posterior) = c(
  "Feature_ID",
  "spatial",
  "phylogenetic", "autotyp_area"
)


# join feature metadata to posterior
trial_posterior = left_join(trial_posterior, parameters_binary, by ="Feature_ID")

# Summarise the posterior distributions
trial_summary = trial_posterior %>%
  group_by(Feature_ID) %>%
  summarise(mean_phylogenetic = mean(phylogenetic),
            mean_spatial = mean(spatial),
            mean_autotyp_area = mean(autotyp_area),
            error_phylogenetic = sd(phylogenetic),
            error_spatial = sd(spatial),
            errro_autotyp_area = sd(autotyp_area),
            domain = first(Main_domain),
            Nichols_1995_prediction = first(Nichols_1995_prediction),
            cor = list(cor(cbind(phylogenetic, spatial))))


trial_summary[,2:4] %>% 
  psych::pairs.panels()
