source("requirements.R")

unusualness_df <-   read.delim(file = "output/unusualness/tables/unusualness.tsv", sep = "\t") %>%
  dplyr::select(Language_ID, prob_lca) 

languages_df <- read_csv(GRAMBANK_LANGUAGES, col_types=LANGUAGES_COLSPEC) %>%		
  dplyr::select(Language_ID = Language_level_ID, Longitude, Latitude, Family_name, Name) %>% 
  distinct(Language_ID, .keep_all = T) %>% 
  mutate(Family_group = ifelse(is.na(Family_name), "Singleton", Family_name)) %>% 
  group_by(Family_name) %>% 
  mutate(n = n()) %>% 
  ungroup %>%
  mutate(Family_group = if_else(n == 1, "Singleton", Family_group)) 

limit <- quantile(unusualness_df$prob_lca, 0.98, na.rm = T)  # top 2 %

h <- unusualness_df %>% 
  ggplot(aes(x = prob_lca, fill=..x..)) +
  geom_histogram(bins = 50) +
  annotate("segment",col="black", alpha = 0.6, x = limit, xend = limit, y = 0, yend = 50, size = 0.5, linetype = "dashed") +
  scale_fill_viridis('prob_lca', option="D", direction = -1) +
  guides(fill="none") +
  xlab("prob_lca") + ylab("Number of Languages") +
  theme_classic() +
  theme(    axis.title = element_text(size = 8),
            panel.background = element_rect(fill = "transparent", colour = NA),  
            plot.background = element_rect(fill = "transparent", colour = NA))

#worldmaps
#rendering a worldmap that is pacific centered
world <- map_data('world', wrap=c(-25,335), ylim=c(-56,80), margin=T)

lakes <- map_data("lakes", wrap=c(-25,335), col="white", border="gray", ylim=c(-55,65), margin=T)

#shifting the longlat of the dataframe to match the pacific centered map
scores_joined_for_plotting <- unusualness_df %>% 
  left_join(languages_df, by = "Language_ID") %>%
  mutate(Longitude = if_else(Longitude <= -25, Longitude + 360, Longitude))

map <- ggplot(scores_joined_for_plotting) +
  geom_polygon(data=world, aes(x=long, y=lat, group=group),
               colour="gray87",
               fill="gray87", size = 0.5) +
  geom_polygon(data=lakes, aes(x=long, y=lat, group=group),
               colour="gray87",
               fill="white", size = 0.3) +
  geom_point(aes(x=Longitude, y=Latitude, color=prob_lca), alpha = 0.6) +
  geom_text_repel(
    data=scores_joined_for_plotting[scores_joined_for_plotting$prob_lca >= limit, ],
    aes(x=Longitude, y=Latitude, label=Name),color="black", max.overlaps = 20
  ) +
  scale_color_viridis('prob_lca', option="D", direction = -1) +
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
  expand_limits(x = scores_joined_for_plotting$Longitude, y = scores_joined_for_plotting$Latitude)

mh <- map + inset_element(h, right = 0.5, bottom = 0.0, left = 0.2, top = 0.4)

png("output/unusualness/plots/surprisal_worldmap.png", height = 5, width = 8, units = "in", res = 300)
plot(mh)
x <- dev.off()