library(tidyverse)
source("global_variables.R")


#The CLDF-format is long, but for PCA, imputation and distances it is better with wide formatted data. 
# This script takes the CLDF dataset and makes it wide. 
# It makes two versions, "with question" and "strict". The "*_with_question" table contains ?-values, "*_strict" does not. See README.md for details.
# If there is more than one dialect for the same language, the dialect with the most features filled out is kept and the other one discarded. This is to optimise for matches to the Jäger-tree used later in the analysis.

#script written by Hedvig Skirgård

#language table for aggregating to language level
languages_df <- read_csv(GRAMBANK_LANGUAGES, col_types = LANGUAGES_COLSPEC) %>% 
  dplyr::select(Language_ID =ID, Language_level_ID, level) 

#output directory
OUTPUTDIR <- file.path('.', "output", 'GB_wide')
if (!dir.exists(OUTPUTDIR)){dir.create(OUTPUTDIR)}

#reading in cldf-formatted long data for grambank
GB <- read_csv(GRAMBANK_VALUES, col_types=VALUES_COLSPEC)

GB_wide_with_question <- GB %>% 
  dplyr::select(Language_ID, Parameter_ID, Value) %>% 
  spread(key = Parameter_ID, value = Value, drop = FALSE) %>% 
  column_to_rownames("Language_ID")

GB_wide_strict <- GB %>% 
  dplyr::select(Language_ID, Parameter_ID, Value) %>% 
  mutate(Value = ifelse(!Value %in% c(0, 1, 2, 3, 4), NA, Value)) %>% 
  spread(key = Parameter_ID, value = Value, drop = FALSE) %>% 
  column_to_rownames("Language_ID")

n_dialects <-   GB_wide_strict %>% 
  rownames_to_column("Language_ID") %>% 
  inner_join(languages_df, by = "Language_ID") %>% 
  filter(level == "dialect") %>% 
  nrow()

n_family <-   GB_wide_strict %>% 
  rownames_to_column("Language_ID") %>% 
  inner_join(languages_df, by = "Language_ID") %>% 
  filter(level == "family") %>% 
  nrow()

cat("There were", nrow(GB_wide_strict), "languages and dialects read in. There were specifically", n_dialects, "dialects and", n_family, "proto-languages.\n") 

GB_wide_strict$na_prop <- apply(GB_wide_strict, 1, function(x) mean(is.na(x)))
GB_wide_with_question$na_prop <- apply(GB_wide_with_question, 1, function(x) mean(is.na(x)))

#for when people want a table with some meta data and data coverage status.
#languages_df <- read_csv(GRAMBANK_LANGUAGES, col_types = LANGUAGES_COLSPEC) %>% 
#  dplyr::select(Language_ID =ID, Language_level_ID, level, Longitude, Latitude, Macroarea, Family_name, Family_level_ID, level) 

#GB_wide_strict %>% 
#rownames_to_column("Language_ID") %>% 
#  dplyr::select(Language_ID, na_prop) %>% 
#  left_join(languages_df) %>% write_tsv("output/GB_wide/GB_languages.tsv", na = "")

GB_wide_strict %>% 
  rownames_to_column("Language_ID") %>%
  arrange(na_prop) %>% 
  left_join(languages_df, by = "Language_ID") %>% 
  dplyr::select(Language_ID, Language_level_ID,na_prop, everything()) %>% 
  distinct(Language_level_ID, .keep_all = T) %>% 
  mutate(Language_ID = Language_level_ID) %>% 
  dplyr::select(-Language_level_ID, -level) %>% 
  write_tsv(file.path(OUTPUTDIR, "GB_wide_strict.tsv"))

cat("Wrote ", OUTPUTDIR, "/GB_wide_strict.tsv \n", sep = "") 

GB_wide_with_question_lg_levelled <- GB_wide_with_question %>% 
  rownames_to_column("Language_ID") %>%
  arrange(na_prop) %>% 
  left_join(languages_df, by = "Language_ID") %>% 
  dplyr::select(Language_ID, Language_level_ID,na_prop, everything()) %>% 
  distinct(Language_level_ID, .keep_all = T) %>% 
  mutate(Language_ID = Language_level_ID) %>% 
  dplyr::select(-Language_level_ID, -level)

GB_wide_with_question_lg_levelled %>% 
  write_tsv(file.path(OUTPUTDIR, "GB_wide_with_question.tsv"))

cat("Wrote ", OUTPUTDIR, "/GB_wide_with_question.tsv \n", sep = "")

cat("All but one dialect has been dropped per language for the analysis. The dialect with the least amount of missing data was kept. There are now", nrow(GB_wide_with_question_lg_levelled),"languages in the (non-imputed) tables in the GB_wide folder.\n") 