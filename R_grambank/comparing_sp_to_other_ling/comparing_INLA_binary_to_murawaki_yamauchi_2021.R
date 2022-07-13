source("requirements.R")

#reading in gb results
dual_summary_effects <- read_tsv("output/spatiophylogenetic_modelling/featurewise/dual_summary_effects.tsv", show_col_types = F)

colnames(dual_summary_effects)

#reading in results from Murawaki and Yamauchi 202
S2 <- read_tsv("comparing_sp_to_other_ling/murawaki_yamauchi_2021_s2.tsv", show_col_types = F) %>% dplyr::select(WALS_Feature_ID, `vi (autologistic-basic)`)
S3 <- read_tsv("comparing_sp_to_other_ling/murawaki_yamauchi_2021_s3.tsv", show_col_types = F) %>% dplyr::select(WALS_Feature_ID, `hi (autologistic-basic)`)
S4 <- read_tsv("comparing_sp_to_other_ling/murawaki_yamauchi_2021_s4.tsv", show_col_types = F) %>% dplyr::select(WALS_Feature_ID, `vi (autologistic-Disjoint)`)
S5 <- read_tsv("comparing_sp_to_other_ling/murawaki_yamauchi_2021_S5.tsv", show_col_types = F) %>% dplyr::select(WALS_Feature_ID, `hi (autologistic-Disjoint)`)

wals_GB_mapping <- read_tsv("comparing_sp_to_other_ling/WALS_to_GB_match.tsv", show_col_types = F) 

murwaki_yamauchi_combined <- S2 %>% 
  full_join(S3,  by = "WALS_Feature_ID") %>% 
  full_join(S4,  by = "WALS_Feature_ID") %>% 
  full_join(S5,  by = "WALS_Feature_ID") %>% 
  full_join(wals_GB_mapping,  by = "WALS_Feature_ID") 

comparison_df <- inner_join(murwaki_yamauchi_combined, dual_summary_effects , by = "Feature_ID") %>% 
  filter(!is.na(Feature_ID)) %>% 
  dplyr::select("Feature_ID", "GB INLA\nPhylogenetic effect (mean)" = "Phylogenetic effect (mean)",  "vi\n(autologistic-basic)" = "vi (autologistic-basic)", "vi\n(autologistic-Disjoint)"= "vi (autologistic-Disjoint)", "GB INLA\nSpatial effect (mean)" = "Spatial effect (mean)", "hi\n(autologistic-basic)"= "hi (autologistic-basic)" , "hi\n(autologistic-Disjoint)" = "hi (autologistic-Disjoint)")

png("output/non_GB_datasets/murwwaki_yamauchi_comparison_INLA_splom.png", height =30, width = 30, units = "cm", res = 400)

comparison_df[,2:7] %>% 
  psych::pairs.panels(             method = "pearson", # correlation method
                                   hist.col = "#a3afd1",# "#a9d1a3","",""),
                                   density = TRUE,  # show density plots
                                   ellipses = F, # show correlation ellipses
                                   cex.labels= 1,
                                   #           smoother= T,
                                   cor=T,
                                   lm=T,
                                   ci = T, cex.cor = 0.9,stars = T)

x <- dev.off()
