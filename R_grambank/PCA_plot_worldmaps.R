#PCA analysis over the imputed data. Here we use the imputed data from random forests to run PCA.
source("requirements.R")

#script written by Hedvig Skirgård and Simon Greenhill

OUTPUTDIR <- file.path('.', 'PCA')

Language_meta_data <-  read_csv(GRAMBANK_LANGUAGES, col_types=LANGUAGES_COLSPEC) %>%		
  dplyr::select(Language_ID = Language_level_ID, Family_name, Name, Longitude, Latitude, Macroarea) %>% 
  distinct(Language_ID, .keep_all = T) %>% 
  mutate(Family_name = ifelse(is.na(Family_name), "Isolate", Family_name))

#reading in the dataframe with PCA values for each language
GB_PCA_df <- read_tsv(file.path(OUTPUTDIR, 'PCA_language_values.tsv'), col_types = cols()) %>%
  left_join(Language_meta_data, by = "Language_ID" ) %>%
  dplyr::select(Language_ID, Longitude, Latitude, everything())


PCA_prop_variance_df <- read_tsv(file.path(OUTPUTDIR, 'PCA_rotations.tsv'), col_types = cols()) %>% 
  distinct(PC,  `Proportion of Variance`) %>% 
  mutate(Title = paste0(PC, " Proportion of Variance: ", `Proportion of Variance`*100, "%" ))  
  
#tided dataframe based on the prcomp object, every row is a component with the propotion of variance information
PCA_summary_importance <- read_tsv(file.path(OUTPUTDIR, 'PCA_rotations.tsv'), col_types = cols()) %>% 
  mutate(`Proportion of Variance` = round(`Proportion of Variance`, digits = 2)) %>%
  dplyr::select(`Proportion of Variance`, everything())

#worldmaps
#rendering a worldmap that is pacific centered
world <- map_data('world', wrap=c(-25,335), ylim=c(-56,80), margin=T)

lakes <- map_data("lakes", wrap=c(-25,335), col="white", border="gray", ylim=c(-55,65), margin=T)

#shifting the longlat of the dataframe to match the pacific centered map
GB_PCA_df_long_shifted <- GB_PCA_df %>%
  mutate(Longitude = if_else(Longitude <= -25, Longitude + 360, Longitude))

#Basemap
basemap <- ggplot(GB_PCA_df_long_shifted) +
  geom_polygon(data=world, aes(x=long, y=lat, group=group),
               colour="gray87",
               fill="gray87", size = 0.5) +
  geom_polygon(data=lakes, aes(x=long, y=lat, group=group),
               colour="gray87",
               fill="white", size = 0.3) +
  scale_fill_viridis() +
  theme(
        legend.position = "none",
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
  expand_limits(x = GB_PCA_df_long_shifted$Longitude, y = GB_PCA_df_long_shifted$Latitude)

## plot PC 1 - 4 on map
filename <- file.path(OUTPUTDIR, "worldmap_PC1.tiff")
cat(paste("writing", filename, "\n"))
pc1 <- basemap +
    geom_point(aes(x=Longitude, y=Latitude, color=PC1), alpha = 0.6) +
    scale_color_viridis() +
    ggtitle(paste("a:", PCA_prop_variance_df[1,3])) +
  theme(title = element_text(face = "bold"))
ggsave(filename, pc1, height = 5, width = 8)


filename <- file.path(OUTPUTDIR, "worldmap_PC2.tiff")
cat(paste("writing", filename, "\n"))
pc2 <- basemap +
    geom_point(aes(x=Longitude, y=Latitude, color=PC2), alpha = 0.6) +
    scale_color_viridis() +
    ggtitle(paste("b:", PCA_prop_variance_df[2,3]))+
  theme(title = element_text(face = "bold"))
ggsave(filename, pc2, height = 5, width = 8)


filename <- file.path(OUTPUTDIR, "worldmap_PC3.tiff")
cat(paste("writing", filename, "\n"))
pc3 <- basemap +
    geom_point(aes(x=Longitude, y=Latitude, color=PC3), alpha = 0.6) +
    scale_color_viridis() +
    ggtitle(paste("c:", PCA_prop_variance_df[3,3]))+
  theme(title = element_text(face = "bold"))
ggsave(filename, pc3, height = 5, width = 8)

filename <- file.path(OUTPUTDIR, "worldmap_all_three_PC.tiff")
cat(paste("writing", filename, "\n"))
map4 <- pc1 / pc2 / pc3
ggsave(filename, map4, height = 16, width = 8)


###Map first 3 PCA components to RGB
GB_PCA_df$RGB <- GB_PCA_df %>%
  dplyr::select(PC1, PC2, PC3) %>%
  sweep(2, apply(., 2, function(x){ 2 * max(abs(x)) }), "/") %>%
  sweep(2, 0.5, "+") %>%
  rgb(alpha = 1)


filename <- file.path(OUTPUTDIR, "PCA_RGB_worldmap.tiff")
filename_png <- file.path(OUTPUTDIR, "PCA_RGB_worldmap.png")

rgbmap <- basemap + geom_jitter(aes(x=Longitude, y=Latitude), color = GB_PCA_df$RGB)
cat(paste("writing", filename, "\n"))
ggsave(filename, rgbmap, height = 5, width = 8)
ggsave(filename_png, rgbmap, height = 5, width = 8)
