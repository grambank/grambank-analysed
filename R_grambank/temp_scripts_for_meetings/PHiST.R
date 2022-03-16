source("requirements.R")

p_load(haplotypes, 
       optparse, 
       phangorn)

#reading in GB
GB <- read.delim(file.path("GB_wide", "GB_cropped_for_missing.tsv"), sep ="\t") 

GB_matrix <- GB %>%
  column_to_rownames("Language_ID") %>%
  as.matrix()

#reading in lg meta data
#areas
if (!dir.exists("non_GB_datasets/glottolog_AUTOTYP_areas.tsv")) { source("unusualness/processing/assigning_AUTOTYP_areas.R") }		
autotyp_area <- read_tsv("non_GB_datasets/glottolog_AUTOTYP_areas.tsv") %>%
  dplyr::select(Language_ID, AUTOTYP_area)

#glottolog-cldf
Language_meta_data <-  read_tsv("non_GB_datasets/glottolog-cldf_wide_df.tsv", col_types = cols()) %>% 
  mutate(Language_level_ID = ifelse(is.na(Language_level_ID), Language_ID, Language_level_ID)) %>% 
  dplyr::select(-Language_ID) %>% 
  dplyr::select(Language_ID  = Language_level_ID, Family_ID, Name, Macroarea) %>% 
  distinct(Language_ID, .keep_all = T) %>% 
  mutate(Family_ID = ifelse(is.na(Family_ID), "Isolate", Family_ID))

#join to ensure exact same order
Language_meta_data <- GB %>% 
  left_join(Language_meta_data, by = "Language_ID") %>% 
  left_join(autotyp_area, by = "Language_ID")

#calculating gower distances %>% 
GB_dist <- GB_matrix %>% 
  cluster::daisy(metric = "gower", warnBin = F) %>% 
  as.matrix()

rownames(GB_dist) <- rownames(GB)
colnames(GB_dist) <- rownames(GB)

GB_dist_2 = GB_dist^2

#PHiST family

phist_family = haplotypes::pairPhiST(x = GB_dist,
                              Language_meta_data$Family_ID,
                              nperm = 99,
                              showprogbar = TRUE
)

phist_family_matrix = phist_family$PhiST

phist_family_list <-phist_family_matrix %>% 
  reshape2::melt() %>% 
  filter(!is.na(value)) %>% 
  unite(Var1, Var2, col = "Vars", sep = " - ") 

phist_family_list$Vars <- fct_reorder(phist_family_list$Vars, phist_family_list$value)

phist_family_list %>% 
  ggplot() +
  geom_point(aes(x = Vars, y = value))

phist_family %>% 
saveRDS("temp_scripts_for_meetings/phist_family.rdata")

#macroarea

phist_macroarea = haplotypes::pairPhiST(x = GB_dist,
                              Language_meta_data$Macroarea,
                              nperm = 99,
                              showprogbar = TRUE
)


phist_macroarea_matrix = phist_macroarea$PhiST

phist_macroarea_list <-phist_macroarea_matrix %>% 
  reshape2::melt() %>% 
  filter(!is.na(value)) %>% 
  unite(Var1, Var2, col = "Vars", sep = " - ") 

phist_macroarea_list$Vars <- fct_reorder(phist_macroarea_list$Vars, phist_macroarea_list$value)

phist_macroarea_list %>% 
  ggplot() +
  geom_point(aes(x = Vars, y = value)) +
  theme_classic() +
  theme(axis.text = element_text(angle = 70, hjust=1))

ggsave("temp_scripts_for_meetings/phist_macroareas.png")

phist_macroarea %>% 
  saveRDS("temp_scripts_for_meetings/phist_macroarea.rdata")


