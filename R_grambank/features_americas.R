source("requirements.R")

Language_meta_data <-  read_tsv("output/non_GB_datasets/glottolog-cldf_wide_df.tsv", show_col_types = F) %>% 
  dplyr::select(Language_ID, level, Family_ID, Name, Macroarea) %>% 
  mutate(americas = ifelse(str_detect(Macroarea, "merica"), "americas", "not americas"))

#reading in GB
GB_americas_prop <- read.delim("output/GB_wide/GB_cropped_for_missing.tsv", sep ="\t") %>% 
  reshape2::melt(id.vars = "Language_ID") %>% 
  left_join(Language_meta_data) %>% 
  filter(!is.na(value)) %>% 
  filter(value == 1) %>% 
  group_by(americas, variable, value) %>% 
  summarise(n = n()) %>% 
  group_by(variable) %>% 
  mutate(sum = sum(n)) %>% 
  ungroup() %>% 
  mutate(prop = n / sum)


#Shorter names for the Grambank parameters that can be used in graphs
Grambank_id_abbrev <- read_tsv(
  file.path("output", "GB_wide", "parameters_binary.tsv"),
  show_col_types=FALSE) %>% 
  rename(variable = ID)


GB_americas_prop %>% 
  group_by(variable) %>%
  slice_max(n) %>% 
  filter(americas == "americas") %>% 
  left_join(Grambank_id_abbrev) %>% 
  dplyr::select(Feature_ID = variable, Grambank_ID_desc, Name)
