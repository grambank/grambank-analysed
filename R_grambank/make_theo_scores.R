source("requirements.R")		

GB_morph_counts <- read_tsv(file = "output/fusion_score/fusion_score.tsv", show_col_types = F) %>% 
  dplyr::select(Language_ID, "Fusion" = mean_morph)


GB_wide <- read_tsv(file.path("output", "GB_wide", "GB_wide_binarized.tsv"), col_types=WIDE_COLSPEC) %>%
  filter(na_prop <= 0.25 )

feature_scores <- data.table::fread(file.path("output", "GB_wide", "parameters_binary.tsv") ,
                                    encoding = 'UTF-8', 
                                    quote = "\"", header = TRUE, 
                                    sep = "\t") %>% 
  filter(Binary_Multistate != "multi") %>% 
  dplyr::select(Parameter_ID = ID, Name, Flexivity,`locus of marking`, `word order`, `Gender/noun class`, informativity) %>% 
  mutate(`Gender/noun class` = as.numeric(`Gender/noun class`)) %>% 
  mutate(Flexivity = as.numeric(Flexivity)) %>% 
  mutate(`locus of marking` = as.numeric(`locus of marking`)) %>% 
  mutate( `word order` = as.numeric( `word order`))

GB_long_for_calc <- GB_wide  %>% 
  dplyr::select(-na_prop) %>% 
  reshape2::melt(id.vars = "Language_ID") %>% 
  dplyr::rename(Parameter_ID = variable) %>%
  left_join(feature_scores, by = "Parameter_ID")

##Flexivity scores
lg_df_for_flex_count <- GB_long_for_calc %>% 
  filter(!is.na(Flexivity)) %>% 
  filter(!is.na(value)) %>% 
  mutate(value_weighted = if_else(Flexivity == 0, abs(value-1), value)) %>% #  reversing the values of the features that have a score of 0
  group_by(Language_ID) %>% 
  dplyr::summarise(`Flexivity` = mean(value_weighted), .groups = "drop_last")

##`locus of marking`s
lg_df_for_HM_DM_count <- GB_long_for_calc %>% 
  filter(!is.na(`locus of marking`)) %>% 
  filter(!is.na(value)) %>% 
  mutate(value_weighted = if_else(`locus of marking` == 0, abs(value-1), value)) %>% # reversing the values of the features that have a score of 0
  group_by(Language_ID) %>% 
  dplyr::summarise(`locus of\nmarking` = mean(value_weighted), .groups = "drop_last")

##`Gender/noun class`_scores
lg_df_for_gender_nc_count <- GB_long_for_calc %>% 
  filter(!is.na(`Gender/noun class`)) %>% 
  filter(!is.na(value)) %>% 
  mutate(value_weighted = if_else(`Gender/noun class` == 0, abs(value-1), value)) %>% #  reversing the values of the features that have a score of 0
  group_by(Language_ID) %>% 
  dplyr::summarise(`Gender/\nnoun class` = mean(value_weighted), .groups = "drop_last")

##OV_VO scores
lg_df_for_OV_VO_count <- GB_long_for_calc %>% 
  filter(!is.na( `word order`)) %>% 
  filter(!is.na(value)) %>% 
  mutate(value_weighted = if_else( `word order` == 0, abs(value-1), value)) %>% # reversing the values of the features that refer to free-standing markers 
  group_by(Language_ID) %>% 
  dplyr::summarise( `word order` = mean(value_weighted), .groups = "drop_last")

##informativity score

lg_df_informativity_score <-  GB_long_for_calc %>% 
  filter(!is.na(informativity)) %>% 
  mutate(value = if_else(Parameter_ID == "GB140", abs(value-1), value)) %>% # reversing GB140 because 0 is the informative state
  group_by(Language_ID, informativity) %>% #grouping per language and per informativity category 
  dplyr::summarise(sum_informativity = sum(value, na.rm = TRUE), #for each informativity cateogry for each langauge, how many are answered 1 ("yes")
                   sum_na = sum(is.na(value)), .groups = "drop_last" ) %>%  #how many of the values per informativity category are missing
  mutate(sum_informativity = ifelse(sum_na >= 1 & sum_informativity == 0, NA, sum_informativity)) %>% #if there is at least one NA and the sum of values for the entire category is 0, the informaitivyt score should be NA because there could be a 1 hiding under the NA value.
  mutate(informativity_score = ifelse(sum_informativity >= 1, 1, sum_informativity)) %>% 
  ungroup() %>% 
  group_by(Language_ID) %>% 
  dplyr::summarise(`Informativity`= mean(informativity_score, na.rm = TRUE), .groups = "drop_last")

#all together now

lg_df_all_scores <- lg_df_for_OV_VO_count %>% 
  full_join(lg_df_for_flex_count, by = "Language_ID") %>% 
  full_join(lg_df_for_gender_nc_count, by = "Language_ID") %>% 
  full_join(lg_df_for_HM_DM_count, by = "Language_ID") %>% 
  full_join(GB_morph_counts, by = "Language_ID") %>% 
  full_join(lg_df_informativity_score, by = "Language_ID") 

lg_df_all_scores %>% 
  write_delim("output/PCA/theo_scores.tsv", na = "", delim = "\t")