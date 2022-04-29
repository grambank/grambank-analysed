source("requirements.R")

p_load(haplotypes, 
       optparse, 
       phangorn)

#reading in GB
GB <- read.delim(file.path("output", "GB_wide", "GB_cropped_for_missing.tsv"), sep ="\t") 

GB_matrix <- GB %>%
  column_to_rownames("Language_ID") %>%
  as.matrix()

#reading in lg meta data
#areas
if (!file.exists("output/non_GB_datasets/glottolog_AUTOTYP_areas.tsv")) { source("unusualness/processing/assigning_AUTOTYP_areas.R") }		
autotyp_area <- read_tsv("output/non_GB_datasets/glottolog_AUTOTYP_areas.tsv", col_types = cols()) %>%
  dplyr::select(Language_ID, AUTOTYP_area)

#glottolog-cldf
glottolog_df <-  read_tsv("output/non_GB_datasets/glottolog-cldf_wide_df.tsv", col_types = cols()) %>% 
  mutate(Language_level_ID = ifelse(is.na(Language_level_ID), Language_ID, Language_level_ID)) %>% 
  dplyr::select(-Language_ID) %>% 
  dplyr::select(Language_ID  = Language_level_ID, Family_ID, Name, Macroarea) %>% 
  distinct(Language_ID, .keep_all = T) %>% 
  mutate(Family_ID = ifelse(is.na(Family_ID), "Isolate", Family_ID))

#join to ensure exact same order
Language_meta_data <- GB %>% 
  left_join(glottolog_df, by = "Language_ID") %>% 
  left_join(autotyp_area, by = "Language_ID")

#calculating gower distances %>% 
GB_dist <- GB_matrix %>% 
  cluster::daisy(metric = "gower", warnBin = F) %>% 
  as.matrix()

rownames(GB_dist) <- rownames(GB)
colnames(GB_dist) <- rownames(GB)


#macroarea

phist_macroarea = haplotypes::pairPhiST(x = GB_dist,
                                        Language_meta_data$Macroarea,
                                        nperm = 99,
                                        showprogbar = TRUE
)


phist_macroarea_matrix <- phist_macroarea$PhiST

phist_macroarea_list <-phist_macroarea_matrix %>% 
  reshape2::melt() %>% 
  filter(!is.na(value)) %>% 
  unite(Var1, Var2, col = "Vars", sep = " - ") 

mean(phist_macroarea_list$value)

phist_macroarea_list$Vars <- fct_reorder(phist_macroarea_list$Vars, phist_macroarea_list$value)

phist_macroarea_list %>% 
  write_tsv("output/dist_fixation_scores/PHiST_macroarea.tsv")

phist_macroarea_list %>% 
  ggplot() +
  geom_bar(aes(x = Vars, y = value, fill = 1 - value), stat = "identity") +
  theme_classic() +
  theme(axis.text.x  = element_text(angle = 70, hjust=1), 
        legend.position = "none", 
        text = element_text(size = 20), 
        axis.title = element_blank())

ggsave("output/dist_fixation_scores/phist_macroareas.png",  width = 7.87 , height =  8.26)

phist_macroarea %>% 
  saveRDS("output/dist_fixation_scores/phist_macroarea.rdata")


#AUTOTYP_area

phist_AUTOTYP_area = haplotypes::pairPhiST(x = GB_dist,
                                        Language_meta_data$AUTOTYP_area,
                                        nperm = 99,
                                        showprogbar = TRUE
)


phist_AUTOTYP_area_matrix <- phist_AUTOTYP_area$PhiST

phist_AUTOTYP_area_list <-phist_AUTOTYP_area_matrix %>% 
  reshape2::melt() %>% 
  filter(!is.na(value)) %>% 
  unite(Var1, Var2, col = "Vars", sep = " - ") 

mean(phist_AUTOTYP_area_list$value)

phist_AUTOTYP_area_list$Vars <- fct_reorder(phist_AUTOTYP_area_list$Vars, phist_AUTOTYP_area_list$value)

phist_AUTOTYP_area_list %>% 
  write_tsv("output/dist_fixation_scores/PHiST_autotyp_area.tsv")

phist_AUTOTYP_area_list %>% 
  ggplot() +
  geom_bar(aes(x = Vars, y = value, fill = 1 - value), stat = "identity") +
  theme_classic() +
  theme(axis.text.x  = element_blank(), 
        legend.position = "none", 
        text = element_text(size = 20), 
        axis.title = element_blank())

ggsave("output/dist_fixation_scores/phist_AUTOTYP_areas.png",  width = 7.87 , height =  8.26)

phist_AUTOTYP_area %>% 
  saveRDS("output/dist_fixation_scores/phist_AUTOTYP_area.rdata")


#Family_ID

phist_Family_ID = haplotypes::pairPhiST(x = GB_dist,
                                           Language_meta_data$Family_ID,
                                           nperm = 99,
                                           showprogbar = TRUE
)


phist_Family_ID_matrix <- phist_Family_ID$PhiST

phist_Family_ID_list <-phist_Family_ID_matrix %>% 
  reshape2::melt() %>% 
  filter(!is.na(value)) %>% 
  unite(Var1, Var2, col = "Vars", sep = " - ") 

mean(phist_Family_ID_list$value)

phist_Family_ID_list$Vars <- fct_reorder(phist_Family_ID_list$Vars, phist_Family_ID_list$value)

phist_Family_ID_list %>% 
  write_tsv("output/dist_fixation_scores/PHiST_Family_ID.tsv")

phist_Family_ID_list %>% 
  ggplot() +
  geom_bar(aes(x = Vars, y = value, fill = 1 - value), stat = "identity") +
  theme_classic() +
  theme(axis.text.x  = element_blank(), 
        legend.position = "none", 
        text = element_text(size = 20), 
        axis.title = element_blank())

ggsave("output/dist_fixation_scores/phist_Family_IDs.png",  width = 7.87 , height =  8.26)

phist_Family_ID %>% 
  saveRDS("temp_scripts_for_meetings/phist_Family_ID.rdata")
