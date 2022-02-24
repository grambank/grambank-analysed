source("requirements.R")

#Glottolog df for language level merging and codes
glottolog_df<- read_tsv("non_GB_datasets/glottolog-cldf_wide_df.tsv", col_types = cols()) %>% 
  mutate(Language_level_ID = ifelse(is.na(Language_level_ID), Language_ID, Language_level_ID)) %>% 
    mutate(Family_ID = ifelse(is.na(Family_ID), Language_level_ID, Family_ID)) %>% 
  dplyr::select(Glottocode, "ISO_639"= ISO639P3code, Language_level_ID, level, Family_ID,Longitude, Latitude) 

e24 <- read_tsv("ethnologue_data/Table_of_Languages.tab", show_col_types = F) %>% 
  dplyr::select("ISO_639", L1_Users, All_Users, EGIDS, Is_Written, Institutional) %>% 
  left_join(glottolog_df)
  
#reading in unusualness scores.
# a low number (-86 compared to -40) makes you RARER
unusualness_score <- read_tsv("unusualness/tables/scores.tsv") %>% 
  dplyr::select(Glottocode = ID, unusualness_score=score)

data <- e24 %>% 
  dplyr::select(Glottocode, L1 = L1_Users , All_Users, EGIDS, Is_Written, Institutional, Longitude, Latitude, Family_ID) %>% 
  left_join(unusualness_score, by = "Glottocode") %>% 
  mutate(L2 = All_Users - L1) %>% #calculating the number of L2 users by subtracting the number of L1 from All users
  mutate(L2_prop = L2/ All_Users) %>% #calculating the proportion of L2 users out of the entire population
  mutate(L1_log10 = log10(L1+1)) %>% 
  mutate(L2_log10 = log10(L2+1)) %>% 
  dplyr::select(Glottocode, unusualness_score, L1, All_Users, L2, L2_prop, L1_log10 ,L2_log10, EGIDS, Is_Written, Institutional, Longitude, Latitude, Family_ID) %>% 
  mutate(Official = if_else(str_detect(EGIDS, "0") | #making a new variable for every language that has EGIDS of 0, 1 or 2
                            str_detect(EGIDS, "1") |  
                            str_detect(EGIDS, "2"), "Official", "Not Official" )) %>% 
  mutate(Official = if_else(str_detect(EGIDS, "10"), "Not Official", Official)) %>%  #making it so that cases where EGIDS is 10 are set to "Not Official"
  filter(!is.na(unusualness_score)) %>% 
  mutate(EGIDS = str_replace_all(EGIDS, "a", "" )) %>% 
  mutate(EGIDS = str_replace_all(EGIDS, "b", "" )) %>% 
  mutate(EGIDS = str_replace_all(EGIDS, "x", "" )) %>% 
  mutate(EGIDS = as.numeric(EGIDS))
  
   
#Having a look at the data overall in a scatter plot matrix
pairs.panels(data[, 2:12], 
             method = "pearson", # correlation method
             hist.col = "#a3afd1",# "#a9d1a3","",""),
             density = TRUE,  # show density plots
             ellipses = F, # show correlation ellipses
             cex.labels= 1,
             #           smoother= T,
             cor=T,
             lm=T,
             ci = T, cex.cor = 0.9,stars = T)

# trees
tree_filename = 'spatiophylogenetic_modelling/processed_data/jaeger_pruned.tree'
phylogenetic_tree = read.tree(tree_filename) 

#pruning data to only obs also in the tree
inner_joined_df <-  phylogenetic_tree$tip.label %>% 
  as.data.frame() %>% 
  dplyr::rename(Glottocode = ".") %>% 
  inner_join(data) 

#pruning tree to only observations in the data
phylogenetic_tree <- phylogenetic_tree %>% keep.tip(inner_joined_df$Glottocode)

#brms
#making a covariance matrix of the tree
vcv_tree <- vcv.phylo(phylogenetic_tree)


###BRMS

formula_for_brms <- unusualness_score ~ L1_log10 + L2_log10 + Is_Written + Official +
                                           (1 | gr(Glottocode, cov = vcv_tree)) +
                                           (L1_log10 + L2_log10 + Is_Written + Official | Family_ID)

full_model <- brms::brm(formula = formula_for_brms,
                    data = filter(inner_joined_df, !is.na(L1), !is.na(L2)),
                    data2 = list(vcv_tree= vcv_tree),
                    iter = 7500,
                    iter = 10000,
                    cores = 4,
                    control = list(adapt_delta =0.99, max_treedepth=15)
) %>% add_criterion("waic")
full_model %>% broom.mixed::tidy() %>% write_csv("unusualness/analysis/full_model.csv")

simplified_model <- brms::brm(unusualness_score ~ 1 + (1 | gr(Glottocode, cov = vcv_tree)),
                    data = filter(inner_joined_df, !is.na(L1), !is.na(L2)),
                    data2 = list(vcv_tree= vcv_tree),
                    iter = 7500,
                    iter = 25000,
                    control = list(adapt_delta =0.99, max_treedepth=15)
) %>% add_criterion("waic")
simplified_model %>% broom.mixed::tidy() %>% write_csv("unusualness/analysis/simplified_model.csv")

loo_compare(full_model, simplified_model, criterion="waic") %>% as.tibble() %>% write_csv("unusualness/analysis/model_comparison.csv")
