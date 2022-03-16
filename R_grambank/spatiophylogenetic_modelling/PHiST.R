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
Language_meta_data <-  read_csv(GRAMBANK_LANGUAGES, col_types=LANGUAGES_COLSPEC) %>%		
  dplyr::select(Language_ID = Language_level_ID, Family_name, Name, Macroarea) %>% 
  distinct(Language_ID, .keep_all = T) %>% 
  mutate(Family_name = ifelse(is.na(Family_name), "Isolate", Family_name))

#join to ensure exact same order
Language_meta_data <- GB %>% 
  left_join(Language_meta_data, by = "Language_ID")

#calculating gower distances %>% 
GB_dist <- GB_matrix %>% 
  cluster::daisy(metric = "gower", warnBin = F) %>% 
  as.matrix()

rownames(GB_dist) <- rownames(GB)
colnames(GB_dist) <- rownames(GB)

GB_dist_2 = GB_dist^2

#PHiST

phist = haplotypes::pairPhiST(x = GB_dist,
                              Language_meta_data$Family_name,
                              nperm = 99,
                              showprogbar = TRUE
)

phist_matrix = phist$PhiST

phist_list <-phist_matrix %>% 
  reshape2::melt() %>% 
  filter(!is.na(value)) %>% 
  unite(Var1, Var2, col = "Vars", sep = " - ") 

phist_list$Vars <- fct_reorder(phist_list$Vars, phist_list$value)

phist_list %>% 
  ggplot() +
  geom_point(aes(x = Vars, y = value))

phist %>% 
saveRDS("temp_scripts_for_meetings/PHiST.rdata")

