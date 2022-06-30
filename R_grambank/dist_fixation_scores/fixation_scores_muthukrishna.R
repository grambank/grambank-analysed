#source("requirements.R")

source("fun_def_h_load.R")

h_load(pkg = c("tidyverse", "reshape2", "ggplot2"))

OUTPUTDIR <- "output/dist_fixation_scores/"
if (!dir.exists(OUTPUTDIR)) {dir.create(OUTPUTDIR)}

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

source("dist_fixation_scores/fun_def_Muthukrishna_2020_CultureFst.r")

#by macroarea
group_df <- Language_meta_data %>% 
  dplyr::select(Language_ID, group = Macroarea)

GB_cropped <- GB %>% 
  left_join(group_df, by = "Language_ID") %>% 
  dplyr::select(-Language_ID) %>% 
  dplyr::select(group, everything())

cat(paste0("Below are the number of languages per Macroarea, smallest 6 groups and then a boxplot of all.\n"))

plot_df <- GB_cropped %>% 
  group_by(group) %>% 
  summarise(n = n())  %>% 
  arrange(n)

plot_df[1:6,]

h_load("txtplot")

txtplot::txtboxplot(plot_df$n, width = 50)

features <- GB_cropped[,-1] %>% colnames()
types <- rep(0, length(features))
names(types) <- features

cfx_object <- CultureFst(d = GB_cropped, loci = features, type = types, bootstrap = T, no.samples = 100, label = NULL) 

cfx_macroarea_matrix <- cfx_object$mean.fst %>% as.matrix()
cfx_macroarea_matrix[upper.tri(x = cfx_macroarea_matrix, diag = T)] <- NA

cfx_macroarea_list <-cfx_macroarea_matrix %>% 
  reshape2::melt() %>% 
  filter(!is.na(value)) %>% 
  unite(Var1, Var2, col = "Vars", sep = " - ") %>% 
  rename(Value_cfx = value)

cfx_macroarea_list$Vars <- fct_reorder(cfx_macroarea_list$Vars, cfx_macroarea_list$Value_cfx)

cfx_macroarea_list %>% 
  ggplot(aes(x = Vars, y = Value_cfx)) +
  geom_bar(aes(fill = 1 - Value_cfx), stat = "identity") +
  theme_classic() +
  theme(text = element_text(angle = 70, hjust = 1, size = 20), 
        legend.position = "None", 
        axis.title.x = element_blank())

mean(cfx_macroarea_list$Value_cfx)
  
ggsave(filename = file.path(OUTPUTDIR, "cfx_barplot_macroarea.png"), height =  7.89, width =  8.61)

cfx_macroarea_list %>% 
  write_tsv(file = file.path(OUTPUTDIR, "cfx_AUTOTYP_macroarea_list.tsv"))

#by autotyp area
group_df <- Language_meta_data %>% 
  dplyr::select(Language_ID, group = AUTOTYP_area)

GB_cropped <- GB %>% 
  left_join(group_df, by = "Language_ID") %>% 
  dplyr::select(-Language_ID) %>% 
  dplyr::select(group, everything())

cat(paste0("Below are the number of languages per AUTOTYP_area, smallest 6 groups and then a boxplot of all.\n"))

plot_df <- GB_cropped %>% 
  group_by(group) %>% 
  summarise(n = n())  %>% 
  arrange(n)

plot_df[1:6,]

txtplot::txtboxplot(plot_df$n, width = 50)

features <- GB_cropped[,-1] %>% colnames()
types <- rep(0, length(features))
names(types) <- features

cfx_object <- CultureFst(d = GB_cropped, loci = features, type = types, bootstrap = T, no.samples = 100, label = NULL) 

cfx_AUTOTYP_area_matrix <- cfx_object$mean.fst %>% as.matrix()
cfx_AUTOTYP_area_matrix[upper.tri(x = cfx_AUTOTYP_area_matrix, diag = T)] <- NA

cfx_AUTOTYP_area_list <-cfx_AUTOTYP_area_matrix %>% 
  reshape2::melt() %>% 
  filter(!is.na(value)) %>% 
  unite(Var1, Var2, col = "Vars", sep = " - ") %>% 
  rename(Value_cfx = value)

cfx_AUTOTYP_area_list$Vars <- fct_reorder(cfx_AUTOTYP_area_list$Vars, cfx_AUTOTYP_area_list$Value_cfx)

cfx_AUTOTYP_area_list %>% 
  ggplot(aes(x = Vars, y = Value_cfx)) +
  geom_bar(aes(fill = 1 - Value_cfx), stat = "identity") +
  theme_classic() +
  theme(text = element_text(angle = 70, hjust = 1, size = 20), 
        legend.position = "None", 
        axis.title.x = element_blank())

mean(cfx_AUTOTYP_area_list$Value_cfx)

ggsave(filename = file.path(OUTPUTDIR, "cfx_barplot_AUTOTYP_area.png"), height =  7.89, width =  8.61)

cfx_AUTOTYP_area_list %>% 
  write_tsv(file = file.path(OUTPUTDIR, "cfx_AUTOTYP_area_list.tsv"))

#families

cut_off_vec <- c(0, 1, 3, 5, 10, 20, 50)

for(i in cut_off_vec) {
  cat("Running the Muthukrishna Cfst for families with a cut-off at ", i, ".\n")
#i <- cut_off_vec[2]
    cut_off <- i

group_df <- Language_meta_data %>% 
  dplyr::select(Language_ID, group = Family_ID) %>% 
  group_by(group) %>% 
  mutate(n = n()) %>% 
  filter(n > cut_off) %>% 
  dplyr::select(-n)

GB_cropped <- GB %>% 
  inner_join(group_df, by = "Language_ID") %>% 
  dplyr::select(-Language_ID) %>% 
  dplyr::select(group, everything())
 
features <- GB_cropped[,-1] %>% colnames()
types <- rep(0, length(features))
names(types) <- features

cfx_object <- CultureFst(d = GB_cropped, loci = features, type = types, bootstrap = T, no.samples = 100, label = NULL) 

cfx_Family_ID_matrix <- cfx_object$mean.fst %>% as.matrix()
cfx_Family_ID_matrix[upper.tri(x = cfx_Family_ID_matrix, diag = T)] <- NA

cfx_Family_ID_list <-cfx_Family_ID_matrix %>% 
  reshape2::melt() %>% 
  filter(!is.na(value)) %>% 
  unite(Var1, Var2, col = "Vars", sep = " - ") %>% 
  rename(Value_cfx = value)

cfx_Family_ID_list$Vars <- fct_reorder(cfx_Family_ID_list$Vars, cfx_Family_ID_list$Value_cfx)

plot_title <- paste0("muthukrishan_cfx_family_id_cut_off_", cut_off)

cfx_Family_ID_list %>% 
  ggplot(aes(x = Vars, y = Value_cfx)) +
  geom_bar(aes(fill = 1 - Value_cfx), stat = "identity") +
  theme_classic() +
  theme(legend.position = "None", 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank()) +
  ggtitle(plot_title) 

mean(cfx_Family_ID_list$Value_cfx)

ggsave(filename = paste0(OUTPUTDIR, "cfx_barplot_Family_ID_cut_off_", cut_off, ".png"), height =  7.89, width =  8.61)

cfx_Family_ID_list %>% 
  write_tsv(file = paste0(OUTPUTDIR, "cfx_Family_ID_list_cut_off_", cut_off, ".tsv"))
}
