source("requirements.R")
#reading in GB

GB <- read.delim(file.path("output", "GB_wide", "GB_cropped_for_missing.tsv"), sep ="\t") %>% 
  dplyr::select(Language_ID)

#glottolog df 
glottolog_df <- read_tsv("output/non_GB_datasets/glottolog-cldf_wide_df.tsv", show_col_types = F) %>% 
  dplyr::select(ISO639P3code, Language_ID, aes)

#Reading in and wrangling data from Bromham et al 2022.
#data taken from supplementary material excel file from Nature's website:
#https://www.nature.com/articles/s41559-021-01604-y#MOESM4

bromham_supp_data <- readxl::read_xlsx("dists/Bromham_et_al_2022_supplementary_data.xlsx", sheet = 4, skip = 2) %>% 
  rename(ISO639P3code = 1 ) %>% 
  dplyr::select(- `1=1-6a, 2=6b, 3=7, 4=8a, 5=8b, 6=9, 7=10`) %>% 
  reshape2::melt(id.vars ="ISO639P3code") %>% 
  left_join(glottolog_df, by = "ISO639P3code" ) %>% 
  left_join(GB, by = "Language_ID") %>% 
  dplyr::select(Language_ID, variable, value) %>% 
  mutate(variable = as.character(variable)) %>% 
  tidyr::separate(col = variable, into = c("level", "period"), sep = "\\.", fill = "right") %>% 
  mutate(period = ifelse(is.na(period), "present", period)) %>% 
  mutate(level = str_extract(level, pattern = "[:digit:]") %>%  as.numeric()) %>%
  group_by(Language_ID, period) %>% 
  summarise(level = weighted.mean(level, value))


glottolog_df$aes <- factor(glottolog_df$aes, levels = c("not_endangered", "shifting", "threatened", "moribund", "nearly_extinct", "extinct", NA))

glottolog_df %>% 
  inner_join(GB) %>% 
  group_by(aes) %>% 
  summarise(n = n()) %>% 
  ggplot() +
  geom_bar(aes(x = aes, y = n, fill = aes), stat = "identity") +
  theme_classic() +
  theme(legend.position = "none")

ggsave("output/dists/endangerment_plot_glottolog.png")


#bromham et al predictions
 #present
present_df <- bromham_supp_data %>% 
  inner_join(GB) %>% 
  filter(period == "present") 

present_df$Language_ID <- fct_reorder(present_df$Language_ID, present_df$level)

present_df %>% 
  ggplot() +
  geom_bar(aes(x = Language_ID, y = level), fill = "orange", stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  ggtitle(label = "Bromham et al endangerment level (present)")

ggsave("output/dists/endangerment_plot_bromham_present.png")

#40

forty_df <- bromham_supp_data %>% 
  inner_join(GB) %>% 
  filter(period == "40") 

forty_df$Language_ID <- fct_reorder(forty_df$Language_ID, forty_df$level)

forty_df %>% 
  ggplot() +
  geom_bar(aes(x = Language_ID, y = level), fill = "darkgreen", stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  ggtitle(label = "Bromham et al endangerment level (40 years in future)")

ggsave("output/dists/endangerment_plot_bromham_40.png")

# 80
eighty_df <- bromham_supp_data %>% 
  inner_join(GB) %>% 
  filter(period == "80") 

eighty_df$Language_ID <- fct_reorder(eighty_df$Language_ID, eighty_df$level)

eighty_df %>% 
  ggplot() +
  geom_bar(aes(x = Language_ID, y = level), fill = "purple", stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  ggtitle(label = "Bromham et al endangerment level (80 years in future)")

ggsave("output/dists/endangerment_plot_bromham_80.png")
