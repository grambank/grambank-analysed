source("requirements.R")

#comparing pc values per lang
pca_langs <- read_tsv("output/PCA/PCA_language_values.tsv", show_col_types = F) %>% 
  dplyr::select(Language_ID, PC1, PC2, PC3)

ppca_langs <- read_tsv("output/PCA/ppca_language_values.tsv", show_col_types = F) %>% 
  dplyr::select(Language_ID, pPC1 = PC1, pPC2 = PC2, pPC3 = PC3)

pca_joined_df <- full_join(pca_langs, ppca_langs, by = "Language_ID") %>% 
  dplyr::select(Language_ID, PC1, pPC1, PC2, pPC2, PC3, pPC3)


png("output/PCA/splom_languages_scores_pca_ppca.png", height =30, width = 30, units = "cm", res = 400)

psych::pairs.panels(pca_joined_df[,-1],
                    method = "pearson", # correlation method
                    hist.col = "#a3afd1",# "#a9d1a3","",""),
                    density = TRUE,  # show density plots
                    ellipses = F, # show correlation ellipses
                    cex.labels= 3,
                    #           smoother= T,
                    cor=T,
                    lm=T,
                    main = "pairwise comparison languages scores PCA vs pPCA",
                    ci = T, cex.cor = 0.9,stars = T
)
x <- dev.off()

#comparing feature loadings
PCA_loadings <- read_tsv("output/PCA/pca_rotations.tsv", show_col_types = F) %>% 
  dplyr::select(PC, Parameter_ID, Contribution) %>% 
  reshape2::dcast(Parameter_ID ~ PC, value.var = "Contribution")

pPCA_loadings <- read_tsv("output/PCA/ppca_rotations.tsv", show_col_types = F) %>% 
  dplyr::select(Parameter_ID, pPC1 = PC1, pPC2 = PC2, pPC3 = PC3)

PCA_loadings_joined <- full_join(PCA_loadings, pPCA_loadings, by = "Parameter_ID") %>% 
  dplyr::select(Parameter_ID, PC1, pPC1, PC2, pPC2, PC3, pPC3)


png("output/PCA/splom_loadings_pca_ppca.png", height =30, width = 30, units = "cm", res = 400)
psych::pairs.panels(PCA_loadings_joined[,-1],
                    method = "pearson", # correlation method
                    hist.col = "#a3afd1",# "#a9d1a3","",""),
                    density = TRUE,  # show density plots
                    ellipses = F, # show correlation ellipses
                    cex.labels= 3,
                    #           smoother= T,
                    cor=T,
                    lm=T,
                    main = "pairwise comparison feature loadings PCA vs pPCA",
                    ci = T, cex.cor = 0.9,stars = T
)

x <- dev.off()



