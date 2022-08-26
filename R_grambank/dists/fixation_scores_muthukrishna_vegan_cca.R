source("requirements.R")
h_load("vegan")

Language_meta_data <-  read_tsv("output/non_GB_datasets/glottolog-cldf_wide_df.tsv", show_col_types = F) %>% 
  dplyr::select(Language_ID, level, Family_ID, Name, Macroarea) %>% 
  mutate(Family_ID = ifelse(is.na(Family_ID), "Isolate", Family_ID)) 
#  filter(!is.na(Macroarea)) %>% 
#  filter(level == "language") %>% 
#    group_by(Family_ID, Macroarea) %>%
#  summarise(n = n(), .groups = "drop") %>% 
#  arrange(desc(n)) %>% 
#  distinct(Family_ID, .keep_all = T)

#sometimes languages in an autotyp-area are not classified all as the same macroarea. we take the macroarea per autotyp-area that is the most common. for example, N Coast Asia has 58 languages glottolog puts in Eurasia and 2 in North Americas. We assign N Coast Asia to Eurasia macroarea

autotyp_area <- read_tsv("output/non_GB_datasets/glottolog_AUTOTYP_areas.tsv", col_types = cols()) %>%
  dplyr::select(Language_ID, AUTOTYP_area)

macroarea_per_autotyp_area <- Language_meta_data %>% 
  left_join(autotyp_area, by = "Language_ID") %>% 
  group_by(AUTOTYP_area, Macroarea) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  full_join(autotyp_area, by = "AUTOTYP_area") %>% 
  group_by(Language_ID, AUTOTYP_area) %>% 
  slice_max(n) %>% 
  dplyr::select(Language_ID, AUTOTYP_area, Macroarea)

Language_meta_data <- Language_meta_data  %>% 
  dplyr::select(-Macroarea) %>% 
  left_join(macroarea_per_autotyp_area, by = "Language_ID") %>% 
  mutate(americas = ifelse(str_detect(Macroarea, "merica"), "americas", "not americas"))

df_long <- read_tsv("output/dists/cfx_AUTOTYP_area_cut_off_0_list.tsv", show_col_types = F) %>% 
  separate("Vars", into = c("Var1", "Var2"), sep = " - ")


veg
