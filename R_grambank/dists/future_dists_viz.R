source("requirements.R")

OUTPUTDIR <- "output/dists/"
if (!dir.exists(OUTPUTDIR)) {dir.create(OUTPUTDIR)}

#reading in GB
GB <- read.delim(file.path("output", "GB_wide", "GB_cropped_for_missing.tsv"), sep ="\t") 

GB_matrix <- GB %>%
  column_to_rownames("Language_ID") %>%
  as.matrix()

#calculating gower distances %>% 
GB_dist <- GB_matrix %>% 
  cluster::daisy(metric = "gower", warnBin = F) %>% 
  as.matrix()

#insert back the names
rownames(GB_dist) <- rownames(GB_matrix)
colnames(GB_dist) <- rownames(GB_matrix)

GB_dist[upper.tri(GB_dist, diag = T)] <- NA

#glottolog df 
glottolog_df <- read_tsv("output/non_GB_datasets/glottolog-cldf_wide_df.tsv", show_col_types = F) %>% 
  dplyr::select(ISO639P3code, Language_ID)

#Reading in and wrangling data from Bromham et al 2022.
#data taken from supplementary material excel file from Nature's website:
#https://www.nature.com/articles/s41559-021-01604-y#MOESM4

bromham_supp_data <- readxl::read_xlsx("dists/Bromham_et_al_2022_supplementary_data.xlsx", sheet = 4, skip = 2) %>% 
  rename(ISO639P3code = 1 ) %>% 
  dplyr::select(- `1=1-6a, 2=6b, 3=7, 4=8a, 5=8b, 6=9, 7=10`) %>% 
  reshape2::melt(id.vars ="ISO639P3code") %>% 
  left_join(glottolog_df, by = "ISO639P3code" ) %>% 
  inner_join(GB, by = "Language_ID") %>% 
  dplyr::select(Language_ID, variable, value) %>% 
  mutate(variable = as.character(variable)) %>% 
  tidyr::separate(col = variable, into = c("level", "period"), sep = "\\.", fill = "right") %>% 
  mutate(period = ifelse(is.na(period), "present", period)) %>% 
  mutate(level = str_extract(level, pattern = "[:digit:]") %>%  as.numeric()) %>%
  group_by(Language_ID, period) %>% 
  summarise(level = weighted.mean(level, value))

#subsetting for languages that have a value lower than 6 for each time period
languages_present <- bromham_supp_data %>% 
  filter(period == "present") %>% 
  filter(level < 6) %>% 
  dplyr::select(Language_ID)

GB_dist_present <- GB_dist[languages_present$Language_ID, languages_present$Language_ID] %>% 
  reshape2::melt() %>% 
  filter(!is.na(value)) %>% 
  dplyr::select(value)

languages_40_years_future <- bromham_supp_data %>% 
  filter(period == "40") %>% 
  filter(level < 6) %>% 
  dplyr::select(Language_ID)

GB_dist_40_years_future <- GB_dist[languages_40_years_future$Language_ID, languages_40_years_future$Language_ID] %>% 
  reshape2::melt() %>% 
  filter(!is.na(value)) %>% 
  dplyr::select(value)

languages_80_years_future <- bromham_supp_data %>% 
  filter(period == "80") %>% 
  filter(level < 6) %>% 
  dplyr::select(Language_ID)

GB_dist_80_years_future <- GB_dist[languages_80_years_future$Language_ID, languages_80_years_future$Language_ID] %>%
  reshape2::melt() %>% 
  filter(!is.na(value)) %>% 
  dplyr::select(value)

mean(GB_dist_80_years_future$value)

#visualisations

joined


GB_dist_present %>% 
  ggplot() +
  geom_histogram(mapping = aes(x = value)) +
  theme_classic()








