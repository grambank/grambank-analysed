source("requirements.R")

OUTPUTDIR <- "output/gramgaps_data/"
if(!dir.exists(OUTPUTDIR)){
  dir.create(OUTPUTDIR)
}

#this is a script for making plots that are based on the mean scores per language

load("output/Gramgaps_data/gaps_with_scores_added.rdata")

glottolog_df <- read_tsv("output/non_GB_datasets/glottolog-cldf_wide_df.tsv", show_col_types = F) %>% 
  dplyr::select(Language_ID, aes)

df <- values_with_scores_added %>% 
  group_by(Language_ID) %>% 
  summarise(mean_answerability = mean(answerability, na.rm = T),
            mean_straightforwardness = mean(straightforwardness, na.rm = T),
            mean_explicitness = mean(explicitness, na.rm = T), 
            mean_described = mean(described, na.rm = T)) 

glottolog_df$aes <- factor(glottolog_df$aes, levels = c("not_endangered", "shifting", "threatened", "moribund", "nearly_extinct", "extinct", NA))

df %>% 
  left_join(glottolog_df, by = "Language_ID") %>%
  ggplot(aes(x = aes, y =mean_described, col= aes )) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle =70, hjust = 1), 
        legend.position = "None") 

ggsave(filename = paste0(OUTPUTDIR, "/descriptive_level_per_aes.png"))

df %>% 
  left_join(glottolog_df, by = "Language_ID") %>%
  ggplot(aes(x = aes, y =mean_answerability, col= aes )) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle =70, hjust = 1), 
        legend.position = "None") 

ggsave(filename = paste0(OUTPUTDIR, "/answerability_level_per_aes.png"))


df %>% 
  left_join(glottolog_df, by = "Language_ID") %>%
  ggplot(aes(x = aes, y =mean_straightforwardness, col= aes )) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle =70, hjust = 1), 
        legend.position = "None") 

ggsave(filename = paste0(OUTPUTDIR, "/straightforwardness_level_per_aes.png"))

df %>% 
  left_join(glottolog_df, by = "Language_ID") %>%
  ggplot(aes(x = aes, y = mean_explicitness, col= aes )) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle =70, hjust = 1), 
        legend.position = "None") 
  
ggsave(filename = paste0(OUTPUTDIR, "/explicitness_level_per_aes.png"))




