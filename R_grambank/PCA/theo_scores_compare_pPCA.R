source("requirements.R")		

#script written by Hedvig Skirgård

cat("Creating plots comparing the loadings of components to theoretical scores.\n")

PCA_df <- read_tsv(file.path("output", "PCA", 'PCA_language_values.tsv'), col_types = cols()) %>% 
  dplyr::select(Language_ID, PC1, PC2, PC3)

pPCA_df <- read_tsv("output/PCA/ppca_language_values.tsv", show_col_types = F) %>% 
  dplyr::select(Language_ID, pPC1 = PC1, pPC2 = PC2, pPC3 = PC3)

theo_scores <- read_tsv("output/PCA/theo_scores.tsv", show_col_types = F)

lg_df_all_scores <- full_join( PCA_df,theo_scores, by = "Language_ID") %>% 
  full_join(pPCA_df, by = "Language_ID")


##SPLOM for overview
tiff("output/PCA/splom_all_scores_PCA_pPCA.tiff", height =30, width = 30, units = "cm", res = 400)

pairs.panels(lg_df_all_scores[,2:13], 
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

