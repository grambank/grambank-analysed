source("requirements.R")
source("unusualness/analysis/unusualness_fun.R")

#script written by Simon Greenhill and Luke Maurits

cat("Calcualting the unusualness score per language and make a map plot.\n")

OUTPUTDIR <- file.path("unusualness")

# Add autotyp areas to language tibble
autotyp_areas <- read_tsv(file.path("non_GB_datasets", "glottolog_AUTOTYP_areas.tsv"), col_types = cols()) %>%
    dplyr::select(ID = Language_ID, AUTOTYP_area) 

languages <- read_csv(GRAMBANK_LANGUAGES, col_types=LANGUAGES_COLSPEC) %>%		
    dplyr::select(ID = Language_level_ID, Longitude, Latitude, Family_name, Name) %>% 
    distinct(ID, .keep_all = T) %>% 
    mutate(Family_group = ifelse(is.na(Family_name), "Singleton", Family_name)) %>% 
    group_by(Family_name) %>% 
    mutate(n = n()) %>% 
    ungroup %>%
    mutate(Family_group = if_else(n == 1, "Singleton", Family_group)) %>%
    left_join(autotyp_areas, by = "ID")

GB <- read_tsv(file.path("GB_wide", "GB_wide_imputed_binarized.tsv"), col_types=WIDE_COLSPEC) %>%
    column_to_rownames("Language_ID")

family_probs <- get_family_probs(GB, languages)
area_probs <- get_area_probs(GB, languages)

scores <- sapply(
    rownames(GB),
    function(lang) { get_lang_score(lang, GB, languages, family_probs, area_probs, c(0.83, 0.17)) }, #combining the score per language based on area and family to one, using a ratio based on the SP-modelling's estimate of the influence of phylogeny vs area.
    simplify=TRUE
)
scores <- as.data.frame(scores)
colnames(scores) <- c('score')
scores$ID <- rownames(scores)

scores_joined <- left_join(scores, languages, by = "ID")

scores_joined %>% 
    write_tsv("unusualness/tables/scores.tsv")

cor_lat_score <- cor.test(abs(scores_joined$Latitude), scores_joined$score, method = "pearson")

cat("The Pearson-correlation between absolute latitude and unsualness score is", round(cor_lat_score$estimate, 2), " with a p-value at", round(cor_lat_score$p.value, 7),"(rounded to 7-decimal places).\n")

cat("Unusualness scores calculated.\n")

#plotting
# figure out limit
limit <- quantile(scores$score, 0.02, na.rm = T)  # top 2 %

h <- ggplot(scores, aes(score, fill=..x..)) +
    geom_histogram(bins = 30) +
    annotate("segment",col="black", alpha = 0.6, x = limit, xend = limit, y = 0, yend = 100, size = 0.5, linetype = "dashed") +
    scale_fill_viridis_c('Score', option="A", direction=-1) +
    xlab("Unusualness Score") + ylab("Number of Languages") +
    guides(fill="none") +
    scale_x_reverse() + #reverse to make more intuitive
    theme_classic() +
    theme(    axis.title = element_text(size = 8),
              panel.background = element_rect(fill = "transparent", colour = NA),  
              plot.background = element_rect(fill = "transparent", colour = NA))

#worldmaps
#rendering a worldmap that is pacific centered
world <- map_data('world', wrap=c(-25,335), ylim=c(-56,80), margin=T)

lakes <- map_data("lakes", wrap=c(-25,335), col="white", border="gray", ylim=c(-55,65), margin=T)

#shifting the longlat of the dataframe to match the pacific centered map
scores_joined_for_plotting <- scores_joined %>%
    mutate(Longitude = if_else(Longitude <= -25, Longitude + 360, Longitude))

map <- ggplot(scores_joined_for_plotting) +
    geom_polygon(data=world, aes(x=long, y=lat, group=group),
                 colour="gray87",
                 fill="gray87", size = 0.5) +
    geom_polygon(data=lakes, aes(x=long, y=lat, group=group),
                 colour="gray87",
                 fill="white", size = 0.3) +
    geom_point(aes(x=Longitude, y=Latitude, color=score), alpha = 0.6) +
    geom_text_repel(
        data=scores_joined_for_plotting[scores_joined_for_plotting$score <= limit, ],
        aes(x=Longitude, y=Latitude, label=Name),color="black", max.overlaps = 20
    ) +
    scale_color_viridis_c('Score', option="A") +
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
    expand_limits(x = scores$Longitude, y = scores$Latitude)

mh <- map + inset_element(h, right = 0.5, bottom = 0.0, left = 0.2, top = 0.4)

plot(mh)

ggsave(file.path(OUTPUTDIR, "plots", 'worldmap.tiff'), mh, height = 5, width = 8)
ggsave(file.path(OUTPUTDIR, "plots", 'worldmap.png'), mh, height = 5, width = 8)

fams <- scores_joined_for_plotting %>% 
    filter(!is.na(Family_name))

p <- ggplot(fams, aes(x=score, y=Family_name, label=Name, color=score)) +
    geom_point() +
    geom_text_repel(
        data = fams[fams$score <= limit, ], max.overlaps = 30 ) +
    scale_color_viridis_c(option="A") +
    ylab(NULL) + xlab("Unusualness Score") +
    guides(color="none") +
    scale_x_reverse() +
    theme_classic()

plot(p)

ggsave(file.path(OUTPUTDIR, "plots", 'perfamily.tiff'), p, width=10, height=16)
