## This script visualises the output from a dual-process model for each Grambank feature
source('requirements.R')

## Change this vector to change the colour palette of the plots
col_vector <- c("#c23c3c","orange", "purple4", "turquoise3")
shape_vector <- c(21, 22, 23, 24)

#make outputdir
if(!dir.exists("output/spatiophylogenetic_modelling/effect_plots")){
  dir.create("output/spatiophylogenetic_modelling/effect_plots")
}

#### Format Posterior Data ####
## Feature Metadata
feature_groupings <- read_csv("feature_grouping_for_analysis.csv", show_col_types = F)

## Identify results by spatial decay
# filename_suffix = "_kappa_2.5_sigma_3..qs"
# filename_suffix = "_kappa_2_sigma_2.RD.qs"
filename_suffix = "kappa_2_sigma_1.15_pcprior0.1"

## Read in model posteriors
posteriors_df <- read_tsv("output/spatiophylogenetic_modelling/featurewise/posteriors_df.tsv", show_col_types = F) %>% 
  filter(str_detect(fn, filename_suffix)) %>% 
  dplyr::select(Feature_ID, spatial = `Precision for spatial_id_in_dual`, phylogenetic = `Precision for phylo_id_in_dual`)

# join feature metadata to posterior
dual_posterior = left_join(posteriors_df, feature_groupings, by ="Feature_ID")

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

## make upper triangle work on log plot
trinf <- data.frame(x=c(0,1,1),y=c(1,0,1))

trinf_sf <- st_linestring(as.matrix(rbind(trinf, trinf[1, ]))) %>%
  st_segmentize(0.01)

trinf <- trinf_sf %>%
  st_coordinates() %>%
  as.data.frame() %>%
  dplyr::select(x = X, y = Y)

#make col for plot labels with a, b, c, d
dual_summary <- dual_summary %>% 
  mutate(letter_plot_label = str_replace(domain, "clause", "a\\) clause")) %>% 
  mutate(letter_plot_label = str_replace(letter_plot_label, "nominal domain", "b\\) nominal domain")) %>% 
  mutate(letter_plot_label = str_replace(letter_plot_label, "pronoun", "c\\) pronoun")) %>% 
  mutate(letter_plot_label = str_replace(letter_plot_label, "verbal domain", "d\\) verbal domain"))

dual_summary_summary <- dual_posterior %>%
  summarise(mean_phylogenetic = mean(phylogenetic),
            mean_spatial = mean(spatial),
            error_phylogenetic = sd(phylogenetic),
            error_spatial = sd(spatial)) %>% 
  t()

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
#  mutate(x = ifelse(x <= 0.001, 0.001, x),
#         y = ifelse(y <= 0.001, 0.001, y)) %>%
#  mutate(x = ifelse(x > 1, 1, x),
#         y = ifelse(y > 1, 1, y)) %>%
  ungroup()


plot_function <- function(label = c("letter_plot_label", "domain", "Nichols_1995_prediction"), facet, fn = spatiophylogenetic_figure_panels_){
  #label <- "Nichols_1995_prediction"
  
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
  geom_point(aes(col = label, fill = label, shape = label),
             size = 2, alpha = 0.6) +
  scale_shape_manual(values = shape_vector) +
  theme_classic(base_size = 10) +
  xlab("Variance explained by Phylogeny") +
  ylab("Variance explained by Geography") +
  scale_colour_manual(values = col_vector) +
  scale_fill_manual(values = col_vector) +
  scale_x_continuous(expand=c(0,0), breaks=c(0.25, 0.5, 0.75, 1), labels = scales::percent_format(scale = 100)) +
  scale_y_continuous(expand=c(0,0), breaks=c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent_format(scale = 100)) +
  coord_equal() +
  geom_path(aes(x, y), data = data.frame(x = seq(0, 1, length.out = 100),
                                         y = seq(1, 0, length.out = 100)),
            linetype = "dashed", color = "#e6e6e6", linewidth = 1.5) +
  geom_polygon(aes(x=x, y=y), data=trinf, fill="#F9F9F9") +
  theme(legend.title = element_blank(),
        panel.spacing = unit(2, "lines"), 
        strip.text = element_text(size = 10),
        plot.margin = ggplot2::margin(t = 5,  # Top margin
                             r = 12,  # Right margin
                             b = 5,  # Bottom margin
                             l = 5)) 


##IFS 
if(facet == T){
center_plot  <- center_plot  +
  lemon::facet_rep_wrap(~label,nrow = 2, repeat.tick.labels = T) + 
#    facet_wrap(~label,nrow = 2) + 
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
       filename = paste0(fn, ".tiff"),
       width = 230 / 2,
       height = 210 / 2,
       units = "mm", dpi = 600)

}

plot_function(label = "letter_plot_label", 
              facet = T, 
              fn = paste0("output/spatiophylogenetic_modelling/effect_plots/spatiophylogenetic_figure_panels_", filename_suffix)
)

plot_function(label = "domain", 
              facet = T, 
              fn = paste0("output/spatiophylogenetic_modelling/effect_plots/spatiophylogenetic_figure_panels_domain_", filename_suffix)
)

plot_function(label = "Nichols_1995_prediction", 
              facet = F, 
              fn = paste0("output/spatiophylogenetic_modelling/effect_plots/spatiophylogenetic_figure_panels_nichols_prediction_", filename_suffix)
)

#tables for supplementary

parameters_binary <- data.table::fread(file.path("output", "GB_wide", "parameters_binary.tsv") ,
                                    encoding = 'UTF-8', 
                                    quote = "\"", header = TRUE, 
                                    sep = "\t")  %>% 
  rename(Feature_ID = ID)

feature_groupings %>% 
  full_join(parameters_binary, by = "Feature_ID") %>% 
  dplyr::select(Feature_ID, Fusion = boundness, Flexivity, "Gender/noun class","locus of marking", "word order", "informativity", Main_domain, Nichols_1995_label, Nichols_1995_prediction) %>% 
  mutate(Fusion = ifelse(Fusion ==0, NA, Fusion)) %>% #we decided to not include 0 for fusion into the metric
write_tsv(file ="output/GB_wide/table_theo_scores_supp.tsv", na = "")
