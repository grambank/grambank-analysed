source("fun_def_h_load.R")
h_load(c("tidyverse", "viridis"))

#read in grambank data
df_fn <- "output/GB_wide/GB_wide_strict.tsv"
if(!file.exists(df_fn)){
  cat(paste0("Making GB data wide, binarising, cropping etc..\n"))
  source("make_wide.R")}
gb <- read_tsv(df_fn, show_col_types = F) %>% 
  dplyr::select(Language_ID, na_prop)
  
#read in meta data
Language_meta_data_fn <- "output/non_GB_datasets/glottolog-cldf_wide_df.tsv"
if(!file.exists(Language_meta_data_fn)){
  source("make_glottolog-cldf_table.R")
}
Language_meta_data  <- read_tsv(Language_meta_data_fn, show_col_types = F) %>% 
dplyr::select(Language_ID, Language_level_ID, Longitude, Latitude, level) %>% 
  mutate(Language_level_ID = ifelse(is.na(Language_level_ID), Language_ID, Language_level_ID)) %>% 
  dplyr::select(-Language_ID) %>% 
  filter(level != "dialect") %>% 
  filter(level != "family") %>% 
  rename(Language_ID = Language_level_ID) %>% 
  distinct(Language_ID, .keep_all = T) 

#combine df
df <- gb %>% 
  inner_join(Language_meta_data, by = "Language_ID")

#rendering a worldmap that is pacific centered
world <- map_data('world', wrap=c(-25,335), ylim=c(-56,80), margin=T)
lakes <- map_data("lakes", wrap=c(-25,335), col="white", border="gray", ylim=c(-55,65), margin=T)

#shifting the longlat of the dataframe to match the pacific centered map
df_long_shifted <- df %>%
  mutate(Longitude = if_else(Longitude <= -25, Longitude + 360, Longitude))

df_long_shifted$prop <- 1 - df_long_shifted$na_prop 
legend_breaks <- c(min(df_long_shifted$prop), 0.8, 0.85, 0.9, 0.95, 1)

#basemap
basemap <- ggplot(df_long_shifted) +
  geom_polygon(data=world, aes(x=long, y=lat, group=group),
               colour="gray87",
               fill="gray87", size = 0.5) +
  geom_polygon(data=lakes, aes(x=long, y=lat, group=group),
               colour="gray87",
               fill="white", size = 0.3) +
  viridis::scale_color_viridis( name='% of data', 
                                breaks = legend_breaks, labels = percent(legend_breaks), begin = 0.5, end = 1) +
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
  expand_limits(x = df_long_shifted$Longitude, y = df_long_shifted$Latitude) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

basemap +
  geom_jitter(mapping = aes(x = Longitude, y = Latitude, color = prop), alpha = 0.6, width = 2)

ggsave(filename = "output/coverage_plots/worldmap_coverage.png", width = 10, height = 6)
