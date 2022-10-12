source("requirements.R")

Language_meta_data <-  read_tsv("output/non_GB_datasets/glottolog-cldf_wide_df.tsv", show_col_types = F) %>% 
  dplyr::select(Language_ID, level, Family_ID, Name, Macroarea) %>% 
  mutate(americas = ifelse(str_detect(Macroarea, "merica"), "americas", "not americas"))

#reading in GB
GB <- read.delim("output/GB_wide/GB_cropped_for_missing.tsv", sep ="\t") %>% 
  reshape2::melt(id.vars = "Language_ID") %>% 
  left_join(Language_meta_data) %>% 
  filter(!is.na(value)) %>% 
  filter(value == 1) %>% 
  group_by(americas, variable, value) %>% 
  summarise(n = n()) %>% 
  group_by(variable) %>% 
  mutate(sum = sum(n)) %>% 
  ungroup() %>% 
  mutate(prop = n / sum) %>% 
  group_by(variable) %>% 
  mutate(max = which.max(n))
