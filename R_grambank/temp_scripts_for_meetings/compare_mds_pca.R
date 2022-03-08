source("requirements.R")

#reading in lg meta data
Language_meta_data <-  read_csv(GRAMBANK_LANGUAGES, col_types=LANGUAGES_COLSPEC) %>%		
  dplyr::select(Language_ID = Language_level_ID, Family_name, Name, Macroarea) %>% 
  distinct(Language_ID, .keep_all = T) %>% 
  mutate(Family_name = ifelse(is.na(Family_name), "Isolate", Family_name))

GB_mds <- read_tsv(file.path("MDS", "MDS_table.tsv"))


GB_PCA_df <- suppressMessages(read_tsv(file.path("PCA", 'PCA_language_values.tsv'))) %>%
  dplyr::select(Language_ID, PC1, PC2, PC3) %>% 
  left_join(Language_meta_data, by = "Language_ID" ) %>%
  dplyr::select(Language_ID, everything())

joined <- full_join(GB_mds, GB_PCA_df)

joined %>% 
  ggplot(aes(x = V1, y = PC1)) +
  geom_point() +
  geom_point(color = "turquoise3") +
  ggpubr::stat_cor(method = "pearson", p.digits = 2, geom = "label", color = "blue",
                   label.y.npc="top", label.x.npc = "left", alpha = 0.8) +
  geom_smooth(method='lm', formula = 'y ~ x') +
  theme_classic() +
  labs(title="",		
       x ="MDS - V1", y = "PCA - PC1")

ggsave("MDS/mds_PCA_1_plot.png")

joined %>% 
  ggplot(aes(x = V2, y = PC2)) +
  geom_point() +
  geom_point(color = "sienna1") +
  ggpubr::stat_cor(method = "pearson", p.digits = 2, geom = "label", color = "blue",
                   label.y.npc="top", label.x.npc = "left", alpha = 0.8) +
  geom_smooth(method='lm', formula = 'y ~ x') +
  theme_classic() +
  labs(title="",		
       x ="MDS - V2", y = "PCA - PC2")

ggsave("MDS/mds_PCA_2_plot.png")


