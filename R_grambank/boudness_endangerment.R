source("requirements.R")

boudness <- read_tsv("output/fusion_score/fusion_score.tsv")

bromham_supp_data <- read_tsv("output/non_GB_datasets/bromham_supp_data.tsv") %>% 
  filter(period == "80")

glottolog_df <- read_tsv("output/non_GB_datasets/glottolog-cldf_wide_df.tsv", show_col_types = F) %>% 
  dplyr::select(ISO639P3code, Language_ID, aes)

glottolog_df$aes <- factor(glottolog_df$aes, levels = c("not_endangered", "shifting", "threatened", "moribund", "nearly_extinct", "extinct", NA))

glottolog_df$aes_num <- as.numeric(glottolog_df$aes)


joined_df <- left_join(boudness, bromham_supp_data) %>% 
  left_join(glottolog_df)

joined_df  %>% 
ggplot(aes(x = aes_num, y = mean_morph)) +
  geom_point() +
  ggpubr::stat_cor() +
  geom_smooth(method='lm', formula = 'y ~ x') 

ggsave(filename = "output/fusion_score/glottolog_aes_level_fusion.png")

joined_df  %>% 
  ggplot(aes(x = level, y = mean_morph)) +
  geom_point() +
  theme_classic() +
  ggpubr::stat_cor(method = "pearson", p.digits = 2, geom = "label", color = "blue",
                   label.y.npc="top", label.x.npc = "left", alpha = 0.8) +
  geom_smooth(method='lm', formula = 'y ~ x') 
  

ggsave(filename = "output/fusion_score/bromham_80_level_fusion.png")
