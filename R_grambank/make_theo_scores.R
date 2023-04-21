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
  dplyr::select(Parameter_ID = ID, Name, Flexivity, Locus_of_Marking, Word_Order, Gender_or_Noun_Class, Informativity) %>% 
  mutate(Gender_or_Noun_Class = as.numeric(Gender_or_Noun_Class)) %>% 
  mutate(Flexivity = as.numeric(Flexivity)) %>% 
  mutate(Locus_of_Marking = as.numeric(Locus_of_Marking)) %>% 
  mutate(Word_Order = as.numeric(Word_Order))

GB_long_for_calc <- GB_wide  %>% 
  dplyr::select(-na_prop) %>% 
  reshape2::melt(id.vars = "Language_ID") %>% 
  dplyr::rename(Parameter_ID = variable) %>%
  left_join(feature_scores, by = "Parameter_ID")

##Flexivity scores
lg_df_for_flex_count <- GB_long_for_calc %>% 
  filter(!is.na(Flexivity)) %>% 
  filter(!is.na(value)) %>% 
  #  reversing the values of the features that have a score of 0
  mutate(value_weighted = if_else(Flexivity == 0, abs(value-1), value)) %>%
  group_by(Language_ID) %>% 
  dplyr::summarise(Flexivity = mean(value_weighted), .groups = "drop_last")

##`locus of marking`s
lg_df_for_HM_DM_count <- GB_long_for_calc %>% 
  filter(!is.na(Locus_of_Marking)) %>% 
  filter(!is.na(value)) %>% 
  # reversing the values of the features that have a score of 0
  mutate(value_weighted = if_else(Locus_of_Marking == 0, abs(value-1), value)) %>%
  group_by(Language_ID) %>% 
  dplyr::summarise(Locus_of_Marking = mean(value_weighted), .groups = "drop_last")

##Gender_or_Noun_Class scores
lg_df_for_gender_nc_count <- GB_long_for_calc %>% 
  filter(!is.na(Gender_or_Noun_Class)) %>% 
  filter(!is.na(value)) %>% 
  #  reversing the values of the features that have a score of 0
  mutate(value_weighted = if_else(Gender_or_Noun_Class == 0, abs(value-1), value)) %>%
  group_by(Language_ID) %>% 
  dplyr::summarise(Gender_or_Noun_Class = mean(value_weighted), .groups = "drop_last")

##OV_VO scores
lg_df_for_OV_VO_count <- GB_long_for_calc %>% 
  filter(!is.na(Word_Order)) %>% 
  filter(!is.na(value)) %>% 
  # reversing the values of the features that refer to free-standing markers 
  mutate(value_weighted = if_else(Word_Order == 0, abs(value-1), value)) %>%
  group_by(Language_ID) %>% 
  dplyr::summarise(Word_Order = mean(value_weighted), .groups = "drop_last")

##informativity score

lg_df_informativity_score <-  GB_long_for_calc %>% 
  filter(!is.na(Informativity)) %>% 
  # reversing GB140 because 0 is the informative state
  mutate(value = if_else(Parameter_ID == "GB140", abs(value-1), value)) %>%
  #grouping per language and per informativity category 
  group_by(Language_ID, Informativity) %>%
  #for each informativity cateogry for each langauge, how many are answered 1 ("yes")
  #how many of the values per informativity category are missing
  dplyr::summarise(sum_informativity = sum(value, na.rm = TRUE),
                   sum_na = sum(is.na(value)), .groups = "drop_last" ) %>%
  #if there is at least one NA and the sum of values for the entire category is 0, the 
  # informativity score should be NA because there could be a 1 hiding under the NA value.
  mutate(sum_informativity = ifelse(sum_na >= 1 & sum_informativity == 0, NA, sum_informativity)) %>%
  mutate(informativity_score = ifelse(sum_informativity >= 1, 1, sum_informativity)) %>%
  ungroup() %>%
  group_by(Language_ID) %>%
  dplyr::summarise(Informativity = mean(informativity_score, na.rm = TRUE), .groups = "drop_last")

#all together now

lg_df_all_scores <- lg_df_for_OV_VO_count %>%
  full_join(lg_df_for_flex_count, by = "Language_ID") %>%
  full_join(lg_df_for_gender_nc_count, by = "Language_ID") %>%
  full_join(lg_df_for_HM_DM_count, by = "Language_ID") %>%
  full_join(GB_morph_counts, by = "Language_ID") %>%
  full_join(lg_df_informativity_score, by = "Language_ID")

lg_df_all_scores %>% 
  write_delim("output/PCA/theo_scores.tsv", na = "", delim = "\t")