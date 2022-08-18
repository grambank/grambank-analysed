source("spatiophylogenetic_modelling/analysis/featurewise/INLA_results_most_prep.R")

#making map plots to show the feature that has the most spatial effect in the model

five_most_spatial_features <- dual_summary %>%
  top_n(n = 5, wt = mean_spatial) %>%
  arrange(desc(mean_spatial)) %>%
  dplyr::select(Feature_ID) %>%
  as.matrix() %>%
  as.vector()

#read in meta data
Language_meta_data_fn <- "output/non_GB_datasets/glottolog-cldf_wide_df.tsv"
if(!file.exists(Language_meta_data_fn)){
  source("make_glottolog-cldf_table.R")
}
Language_meta_data  <- read_tsv(Language_meta_data_fn, show_col_types = F) %>%
  dplyr::select(Language_ID, Longitude, Latitude)


df_for_maps <- GB_df %>%
  inner_join(Language_meta_data, by = "Language_ID") %>%
  mutate(Longitude = if_else(Longitude <= -25, Longitude + 360, Longitude)) #shift for pacific

#basemap
#rendering a worldmap that is pacific centered
world <- map_data('world', wrap=c(-25,335), ylim=c(-56,80), margin=T)
lakes <- map_data("lakes", wrap=c(-25,335), col="white", border="gray", ylim=c(-55,65), margin=T)

basemap <- ggplot(df_for_maps) +
  geom_polygon(data=world, aes(x=long, y=lat, group=group),
               colour="gray87",
               fill="gray87", size = 0.5) +
  geom_polygon(data=lakes, aes(x=long, y=lat, group=group),
               colour="gray87",
               fill="white", size = 0.3) +
  theme(
    # all of these lines are just removing default things like grid lines, axes etc
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  coord_map(projection = "vandergrinten", ylim=c(-56,67)) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

index <- 0

for(feature in five_most_spatial_features){
  
  #  feature <- five_most_spatial_features[1]
  
  index <- index + 1
  #generating plot title
  plot_title <- GB_id_desc %>%
    filter(Feature_ID == feature) %>%
    dplyr::select(Name) %>%
    as.matrix() %>%
    as.vector()
  
  plot_title <- paste0(feature, " ", plot_title)
  plot_title <- gsub('(.{1,60})(\\s|$)', '\\1\n', plot_title)
  
  #filename
  filename <- paste(OUTPUTDIR, "/most_signal_spatial_",index, "_" , feature, ".png", sep = "")
  
  feature_df <-df_for_maps %>%
    filter(!is.na("Longitude")) %>%
    dplyr::select(Language_ID, all_of(feature), Longitude, Latitude) %>%
    rename(Feature = 2) %>%
    mutate(point.color = as.character(Feature)) %>%
    mutate(point.color = str_replace_all(point.color, "1", color_vector[1])) %>%
    mutate(point.color = str_replace_all(point.color, "0", color_vector[2])) %>%
    mutate(point.color = ifelse(is.na(point.color), color_vector[3], point.color)) %>%
    arrange(Feature)
  
  basemap +
    geom_jitter(data = feature_df, mapping = aes(x = Longitude, y = Latitude),  color = feature_df$point.color, alpha = 0.6, width = 3) +
    ggtitle(plot_title) +
    theme(plot.title = element_text(size = 18, face = "bold"))
  
  ggsave(filename = filename, width = 9.3, height = 9.2, units = "in")
  
}
