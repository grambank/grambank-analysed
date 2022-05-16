source("requirements.R")

#this script takes the data where the comments have been classified and outputs a rdata object of that, as well as one that is summarised per feature

OUTPUTDIR <- "output/gramgaps_data/"
if(!dir.exists(OUTPUTDIR)){
  dir.create(OUTPUTDIR)
}

load("output/Gramgaps_data/values_tagged_for_comment_class.RData")
parameters <- read_csv("../grambank/cldf/parameters.csv", show_col_types = F) %>% 
  rename(Parameter_ID = ID) 
groupings <- read_csv("feature_grouping_for_analysis.csv", show_col_types = F) %>% rename(Parameter_ID = Feature_ID)

languages_df <- read_csv("../grambank/cldf/languages.csv", show_col_types = F) %>% 
  filter(level != "Family") %>% #probably best to remove the three proto-langauges in the set
  dplyr::select(Language_ID = ID) 

# join gaps with parameters and groupings
gaps %>%
  inner_join(languages_df, by = "Language_ID") %>% #probably best to remove the three proto-langauges in the set
  left_join(parameters, by = "Parameter_ID") -> gaps

# ANSWERABILITY
gaps %>% 
  mutate(answerability = ifelse(Value >= 0, 1, 0)) -> gaps


gaps %>% mutate(straightforwardness = case_when(is.na(Comment_class) & Value >= 0 ~ 1,
                                                Value >= 0 & Comment_class == "note on references or variety" ~ 1,
                                                Value >= 0 & Comment_class == "no category" ~ 1,
                                                Comment_class == "specific" ~ 0)) -> gaps

# TODO Jakob can you check some NA's here to be on the safe side

# explicitness
# If Comment_class = "not mentioned" → 0
# If value = ? and Comment_class = "specific" → 0
# If value >=0 and Comment_class = specific → 1
# If value >=0 and Comment_class = note on references → 1
# If value >=0 and Comment_class = NA → 1
# The rest of the combinations (8 remaining) don't affect explicitness 
gaps %>% mutate(explicitness = case_when(Comment_class == "not mentioned" ~ 0,
                                         Value == "?" & Comment_class == "specific" ~ 0,
                                         Value >= 0 & Comment_class == "specific" ~ 1,
                                         Value >= 0 & Comment_class == "note on references" ~ 1,
                                         is.na(Comment_class) & Value >= 0~ 1)) -> gaps

# DegreeColumn 4: described
# 
# - If Comment_class = "not mentioned" → 0
# - If value >=0 and Comment_class = specific → 1
# - If value >=0 and Comment_class = note on references → 1
# - If value >=0 and Comment_class = NA → 1
# - The rest of the combinations (10 remaining) don't affect "described"
         
gaps %>% mutate(described = case_when(Comment_class == "not mentioned" ~ 0,
                                         Value >= 0 & Comment_class == "specific" ~ 1,
                                         Value >= 0 & Comment_class == "note on references" ~ 1,
                                         is.na(Comment_class) & Value >= 0~ 1)) -> gaps
values_with_scores_added <- gaps

save(values_with_scores_added, file = paste0(OUTPUTDIR, "/gaps_with_scores_added.rdata"))

values_with_scores_summarised_per_feature <- values_with_scores_added %>% 
  group_by(Parameter_ID) %>% 
  summarise(mean_answerability = mean(answerability, na.rm = T),
            mean_straightforwardness = mean(straightforwardness, na.rm = T),
            mean_explicitness = mean(explicitness, na.rm = T), 
            mean_described = mean(described, na.rm = T)) 

values_with_scores_summarised_per_feature %>% 
  write_tsv(paste0(OUTPUTDIR, "/values_with_scores_summarised_per_feature.tsv"))