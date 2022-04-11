source("requirements.R")		

#script written by Hedvig Skirgård

cat("Creating plots comparing the loadings of components to theoretical scores.\n")

GB_morph_counts <- read_tsv(file.path("output", "Bound_morph", "Bound_morph_score.tsv"), col_types = cols()) %>% 
  dplyr::select(Language_ID, "Boundness" = mean_morph)

#comparison to PC1
PCA_df <- read_tsv(file.path("output", "PCA", 'PCA_language_values.tsv'), col_types = cols()) %>% 
  dplyr::select(Language_ID, PC1, PC2, PC3)

df_morph_count_PCA <- GB_morph_counts %>% 
  dplyr::select(Language_ID, "Boundness") %>% 
  full_join(PCA_df, by = "Language_ID")

df_morph_count_PCA__plot <- df_morph_count_PCA %>% 
  ggplot(aes(Boundness, PC1)) +
  geom_point(color = "turquoise3") +
  ggpubr::stat_cor(method = "pearson", p.digits = 2, geom = "label", color = "blue",
                   label.y.npc="top", label.x.npc = "left", alpha = 0.8) +
  geom_smooth(method='lm', formula = 'y ~ x') +
  theme_classic() +
  labs(title="",		
       x ="Boundness score", y = "PC1")

tiff(file.path("output", "PCA", "PC1_bound_morph_cor_plot.tiff"), width = 5, height = 4,  units = "in", res = 300)
plot(df_morph_count_PCA__plot)
x <- dev.off()

png(file.path("output", "PCA", "PC1_bound_morph_cor_plot.png"), width = 5, height = 4,  units = "in", res = 300 )
plot(df_morph_count_PCA__plot)
x <- dev.off()

##Calculating the other scores per language

GB_wide <- read_tsv(file.path("output", "GB_wide", "GB_wide_binarized.tsv"), col_types=WIDE_COLSPEC) %>%
  filter(na_prop <= 0.25 )

feature_scores <- data.table::fread(file.path("output", "GB_wide", "parameters_binary.tsv") ,
                                    encoding = 'UTF-8', 
                                    quote = "\"", header = TRUE, 
                                    sep = "\t") %>% 
    filter(Binary_Multistate != "multi") %>% 
  dplyr::select(Parameter_ID = ID, Flexivity,`locus of marking`, `word order`, `Gender/noun class`, informativity) %>% 
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
  mutate(value_weighted = if_else(Flexivity == 0, abs(value-1), value)) %>% # reversing the values of the features that refer to free-standing markers 
  group_by(Language_ID) %>% 
  dplyr::summarise(`Flexivity` = mean(value_weighted), .groups = "drop_last")

##`locus of marking`s
lg_df_for_HM_DM_count <- GB_long_for_calc %>% 
  filter(!is.na(`locus of marking`)) %>% 
  filter(!is.na(value)) %>% 
  mutate(value_weighted = if_else(`locus of marking` == 0, abs(value-1), value)) %>% # reversing the values of the features that refer to free-standing markers 
  group_by(Language_ID) %>% 
  dplyr::summarise(`locus of\nmarking` = mean(value_weighted), .groups = "drop_last")

##`Gender/noun class`_scores
lg_df_for_gender_nc_count <- GB_long_for_calc %>% 
  filter(!is.na(`Gender/noun class`)) %>% 
  filter(!is.na(value)) %>% 
  mutate(value_weighted = if_else(`Gender/noun class` == 0, abs(value-1), value)) %>% # reversing the values of the features that refer to free-standing markers 
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
  full_join(lg_df_informativity_score, by = "Language_ID") %>% 
  full_join(PCA_df, by = "Language_ID")  
    
##SPLOM for overview
tiff("output/PCA/splom_all_scores.tiff", height =30, width = 30, units = "cm", res = 400)

pairs.panels(lg_df_all_scores[,2:10], 
             method = "pearson", # correlation method
             hist.col = "#a3afd1",# "#a9d1a3","",""),
             density = TRUE,  # show density plots
             ellipses = F, # show correlation ellipses
             cex.labels= 1,
             #           smoother= T,
             cor=T,
             lm=T,
             ci = T, cex.cor = 0.9,stars = T
)
x <- dev.off()


png("output/PCA/splom_all_scores.png", height =30, width = 30, units = "cm", res = 400)

pairs.panels(lg_df_all_scores[,2:10], 
             method = "pearson", # correlation method
             hist.col = "#a3afd1",# "#a9d1a3","",""),
             density = TRUE,  # show density plots
             ellipses = F, # show correlation ellipses
             cex.labels= 1,
             #           smoother= T,
             cor=T,
             lm=T,
             ci = T, cex.cor = 0.9,stars = T
)
x <- dev.off()

cat("Scatterplot of PC and theoretical scores made.\n")


#specific scatterplots for MS

PCA2_vs_gender_plot <- lg_df_all_scores  %>% 
  ggplot(aes(`Gender/\nnoun class` , PC2)) +
  geom_point(color = "turquoise3") +
  ggpubr::stat_cor(method = "pearson", p.digits = 2, geom = "label", color = "blue",
                   label.y.npc="top", label.x.npc = "left", alpha = 0.8) +
  geom_smooth(method='lm', formula = 'y ~ x') +
  theme_classic() +
  labs(title="",		
       x ="Gender/\nnoun class", y = "PC2") + 		
  xlim(c(0,max(lg_df_all_scores$`Gender/
noun class`)))

tiff(file.path("output", "PCA", "PC2_gender_cor_plot.tiff"), width = 5, height = 4,  units = "in", res = 300)
plot(PCA2_vs_gender_plot)
x <- dev.off()

png(file.path("output", "PCA", "PC2_gender_cor_plot.png"), width = 5, height = 4,  units = "in", res = 300)
plot(PCA2_vs_gender_plot)
x <- dev.off()