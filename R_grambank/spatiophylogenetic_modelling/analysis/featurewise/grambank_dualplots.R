## This script visualises the output from a dual-process model for each Grambank feature
source('requirements.R')

## Change this vector to change the colour palette of the plots
col_vector <- c("#039e37", "purple4",  "#c23c3c", "turquoise3")

#make outputdir
if(!dir.exists("output/spatiophylogenetic_modelling/effect_plots")){
  dir.create("output/spatiophylogenetic_modelling/effect_plots")
}

#### Format Posterior Data ####
## Feature Metadata
feature_groupings <- read_csv("feature_grouping_for_analysis.csv", show_col_types = F)

## Read in model posteriors
model_output_files = list.files(path = "output/spatiophylogenetic_modelling/featurewise/",
                                pattern = "*.qs",
                                full.names = TRUE)
model_output_files <- model_output_files[!str_detect(model_output_files, "kapp")]

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
dual_posterior = left_join(dual_posterior, feature_groupings, by ="Feature_ID")

# Summarise the posterior distributions
dual_summary = dual_posterior %>%
  group_by(Feature_ID) %>%
  summarise(mean_phylogenetic = mean(phylogenetic),
            mean_spatial = mean(spatial),
            error_phylogenetic = sd(phylogenetic),
            error_spatial = sd(spatial),
            domain = first(Main_domain),
            Nichols_1995_prediction = first(Nichols_1995_prediction),
            cor = list(cor(cbind(phylogenetic, spatial))))

dual_summary %>% 
  dplyr::select(Feature_ID, `Phylogenetic effect (mean)` = mean_phylogenetic, `Phylogenetic effect (Standard Deviation)`=error_phylogenetic, `Spatial effect (mean)`= mean_spatial, `Spatial effect (Standard Deviation)` = error_spatial) %>% 
  mutate_if(is.numeric, round,3) %>% 
  arrange(desc(`Phylogenetic effect (mean)`)) %>% 
  write_tsv("output/spatiophylogenetic_modelling/featurewise/dual_summary_effects.tsv", na = "")

trinf <- data.frame(x=c(0,1,1),y=c(1,0,1))

## make upper triangle work on log plot
trinf_sf <- st_linestring(as.matrix(rbind(trinf, trinf[1, ]))) %>%
  st_segmentize(0.01)

trinf <- trinf_sf %>%
  st_coordinates() %>%
  as.data.frame() %>%
  dplyr::select(x = X, y = Y)

#make col for plot labels with a, b, c, d
dual_summary <- dual_summary %>% 
  mutate(letter_plot_label = str_replace(domain, "clause", "a")) %>% 
  mutate(letter_plot_label = str_replace(letter_plot_label, "nominal domain", "b")) %>% 
  mutate(letter_plot_label = str_replace(letter_plot_label, "pronoun", "c")) %>% 
  mutate(letter_plot_label = str_replace(letter_plot_label, "verbal domain", "d"))

## ellipses
ellipses <- dual_summary %>%
  rowwise() %>%
  summarise(Feature_ID = Feature_ID,
            letter_plot_label =letter_plot_label,
            domain = domain,
            Nichols_1995_prediction = Nichols_1995_prediction,
            ellipse = ellipse::ellipse(cor, centre = c(mean_phylogenetic, mean_spatial),
                                       scale = c(error_phylogenetic, error_spatial),
                                       level = 0.68)) %>%
  mutate(x = ellipse[ , 1],
         y = ellipse[ , 2]) %>%
  mutate(x = ifelse(x <= 0.001, 0.001, x),
         y = ifelse(y <= 0.001, 0.001, y)) %>%
  mutate(x = ifelse(x > 1, 1, x),
         y = ifelse(y > 1, 1, y)) %>%
  ungroup()

#h_load("lemon")

plot_function <- function(label = c("letter_plot_label", "domain", "Nichols_1995_prediction"), facet, fn = spatiophylogenetic_figure_panels_ellipses){
  #label <- "letter_plot_label"
  
  dual_summary <- dual_summary %>% 
    dplyr::select(Feature_ID, all_of(label), mean_phylogenetic, mean_spatial) %>% 
    rename(label = 2) %>% 
    filter(!is.na(label))
  
  ellipses <- ellipses %>% 
    dplyr::select(Feature_ID, all_of(label), x, y) %>% 
    rename(label = 2) %>% 
    filter(!is.na(label))
  
center_plot =   ggplot(data = dual_summary,
                       aes(x = mean_phylogenetic,
                           y = mean_spatial)) +
  geom_polygon(aes(x, y,
                   fill = label,
                   group = Feature_ID),
               data = ellipses,
               alpha = 0.1,
               color = NA) +
  geom_point(aes(col = label, fill = label),
             size = 0.5) +
  theme_classic(base_size = 10) +
  xlab("Variance explained by Phylogeny") +
  ylab("Variance explained by Geography") +
  scale_colour_manual(values = col_vector) +
  scale_fill_manual(values = col_vector) +
  scale_x_continuous(expand=c(0,0), breaks=c(0.25, 0.5, 0.75, 1), labels = scales::percent_format(scale = 100)) +
  scale_y_continuous(expand=c(0,0), breaks=c(0, 0.25, 0.5, 0.7), limits = c(0, 0.7), labels = scales::percent_format(scale = 100)) +
  coord_equal() +
  geom_path(aes(x, y), data = data.frame(x = seq(0, 1, length.out = 100),
                                         y = seq(1, 0, length.out = 100)),
            linetype = "dashed", color = "#e6e6e6", size = 1.5) +
  geom_polygon(aes(x=x, y=y), data=trinf, fill="#F9F9F9") +
  theme(legend.title = element_blank(),
        panel.spacing = unit(2, "lines"), 
        strip.text = element_text(size = 10)
        ) 

if(facet == T){
center_plot  <- center_plot  +
#  lemon::facet_rep_wrap(~label,nrow = 2, repeat.tick.labels = T) + 
    facet_wrap(~label,nrow = 2) + 
    theme(legend.position = "None",
        strip.background = element_blank())
}


if(facet == F ){
  center_plot  <- center_plot  +
    geom_point(aes(shape = label)) 
  }


if(str_detect(fn, "ichols")){
  center_plot  <- center_plot  +
    ggtitle(label = "Nichols (1995) predictions")
  }

plot(center_plot)

ggsave(plot = center_plot,
       filename = paste0("output/spatiophylogenetic_modelling/effect_plots/", fn, ".jpg"),
       width = 230 / 2,
       height = 210 / 2,
       units = "mm")

}

plot_function(label = "letter_plot_label", facet = T, fn = "spatiophylogenetic_figure_panels_ellipses")

plot_function(label = "domain", facet = T, fn = "spatiophylogenetic_figure_panels_ellipses_domain")

plot_function(label = "Nichols_1995_prediction", facet = F, fn = "spatiophylogenetic_figure_panels_ellipses_nichols_prediction")

#tables for supplementary

parameters_binary <- data.table::fread(file.path("output", "GB_wide", "parameters_binary.tsv") ,
                                    encoding = 'UTF-8', 
                                    quote = "\"", header = TRUE, 
                                    sep = "\t")  %>% 
  rename(Feature_ID = ID)

feature_groupings %>% 
  full_join(parameters_binary, by = "Feature_ID") %>% 
  dplyr::select(Feature_ID, boundness, Flexivity, "Gender/noun class","locus of marking", "word order", "informativity", Main_domain, Nichols_1995_label, Nichols_1995_prediction) %>% 
write_tsv(file ="output/GB_wide/table_theo_scores_supp.tsv", na = "")
