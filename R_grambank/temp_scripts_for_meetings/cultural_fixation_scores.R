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
Language_meta_data <-  read_tsv("output/non_GB_datasets/glottolog-cldf_wide_df.tsv", col_types = cols()) %>% 
  mutate(Language_level_ID = ifelse(is.na(Language_level_ID), Language_ID, Language_level_ID)) %>% 
  dplyr::select(-Language_ID) %>% 
  dplyr::select(Language_ID  = Language_level_ID, Family_ID, Name, Macroarea) %>% 
  distinct(Language_ID, .keep_all = T) %>% 
  mutate(Family_ID = ifelse(is.na(Family_ID), "Isolate", Family_ID))

#join to ensure exact same order
Language_meta_data <- GB %>% 
  left_join(Language_meta_data, by = "Language_ID") %>% 
  left_join(autotyp_area, by = "Language_ID")

#using Muthukrishna's et al's approach

source("temp_scripts_for_meetings/Muthukrishna_2020_CultureFst.r")

group_df <- Language_meta_data %>% 
  dplyr::select(Language_ID, Macroarea)

GB_cropped <- GB %>% 
  left_join(group_df, by = "Language_ID") %>% 
  dplyr::select(-Language_ID) %>% 
  dplyr::select(Macroarea, everything())

features <- GB_cropped[,-1] %>% colnames()
types <- rep(0, length(features))
names(types) <- features

cfx_object <- CultureFst(d = GB_cropped, loci = features, type = types, bootstrap = T, no.samples = 100, label = "output/CFx_test") 

cfx_macroarea_matrix <- cfx_object$mean.fst %>% as.matrix()
cfx_macroarea_matrix[upper.tri(x = cfx_macroarea_matrix, diag = T)] <- NA

cfx_macroarea_list <-cfx_macroarea_matrix %>% 
  reshape2::melt() %>% 
  filter(!is.na(value)) %>% 
  unite(Var1, Var2, col = "Vars", sep = " - ") %>% 
  rename(Value_cfx = value)


joined <- phist_macroarea_list %>% 
  rename(Value_PHiST = value) %>% 
  full_join(cfx_macroarea_list)

joined %>% 
  ggplot() +
  geom_point(aes(x = Value_PHiST, y = Value_cfx)) +
  theme_minimal()

ggsave("comparison_phist_cfx.png")

cor.test(joined$Value_PHiST, joined$Value_cfx)
