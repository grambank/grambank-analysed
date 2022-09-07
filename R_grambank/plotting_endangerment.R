source("requirements.R")

tree <- ape::read.tree("output/spatiophylogenetic_modelling/processed_data/EDGE_pruned_tree.tree")

bromham_fn <- "output/non_GB_datasets/bromham_supp_data.tsv"
if(!file.exists(bromham_fn)){
  source("dists/lang_lists_endangerment.R")
}
bromham_df <- read_tsv(file = bromham_fn, show_col_types = F)

#subsetting to overlapping sets
bromham_df_tree <- tree$tip.label %>% 
  as.data.frame() %>% 
  rename(Language_ID = ".") %>% 
  inner_join(bromham_df, by = "Language_ID")

tree <- ape::keep.tip(tree, tip = bromham_df_tree$Language_ID)

bromham_df_80 <- bromham_df_tree %>% 
  filter(period == "80") %>% 
  column_to_rownames("Language_ID") %>% 
  dplyr::select(level)

h_load("ggtree")

circ <- ggtree(tree, layout = "circular")

gheatmap(circ, bromham_df_80, width = 0.15, colnames = F)


#map

Language_meta_data_fn <- "output/non_GB_datasets/glottolog-cldf_wide_df.tsv"
if(!file.exists(Language_meta_data_fn)){
  source("make_glottolog-cldf_table.R")
}
Language_meta_data  <- read_tsv(Language_meta_data_fn, show_col_types = F) %>% 
  dplyr::select(Language_ID, Longitude, Latitude) %>% 
  mutate(Longitude = if_else(Longitude <= -25, Longitude + 360, Longitude))

world <- map_data('world', wrap=c(-25,335), ylim=c(-56,80), margin=T)
lakes <- map_data("lakes", wrap=c(-25,335), col="white", border="gray", ylim=c(-55,65), margin=T)

bromham_df_80_for_map <- bromham_df %>% 
  left_join(Language_meta_data) %>% 
  filter(period == "80") 

basemap <- ggplot(bromham_df_80_for_map) +
  geom_polygon(data=world, aes(x=long, y=lat, group=group),
               colour="gray87",
               fill="gray87", size = 0.5) +
  geom_polygon(data=lakes, aes(x=long, y=lat, group=group),
               colour="gray87",
               fill="white", size = 0.3) +
  viridis::scale_color_viridis( name='% of data', 
                                breaks = 0.25*0:4, labels = percent(0.25*0:4), direction = -1) +
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
  coord_map(projection = "vandergrinten", ylim=c(-56,67)) 

basemap +
  geom_jitter(mapping = aes(x = Longitude, y = Latitude, color = level), alpha = 0.6, width = 2)



