source("requirements.R")		

#script written by Hedvig Skirgård

cat("Creating plots comparing the loadings of components to theoretical scores.\n")

PCA_df <- read_tsv(file.path("output", "PCA", 'PCA_language_values.tsv'), col_types = cols()) %>% 
  dplyr::select(Language_ID, PC1, PC2, PC3)

theo_scores <- read_tsv("output/PCA/theo_scores.tsv", show_col_types = F)

lg_df_all_scores <- full_join( PCA_df,theo_scores, by = "Language_ID")

#comparison to PC1

df_morph_count_PCA__plot <- lg_df_all_scores  %>% 
  ggplot(aes(Fusion, PC1)) +
  geom_point(color = "turquoise3") +
  ggpubr::stat_cor(method = "pearson", p.digits = 2, geom = "label", color = "blue",
                   label.y.npc="top", label.x.npc = "left", alpha = 0.8) +
  geom_smooth(method='lm', formula = 'y ~ x') +
  theme_classic() +
  labs(title="",		
       x ="Fusion score", y = "PC1")

tiff(file.path("output", "PCA", "PC1_fusion__morph_cor_plot.tiff"), width = 5, height = 4,  units = "in", res = 300)
plot(df_morph_count_PCA__plot)
x <- dev.off()

png(file.path("output", "PCA", "PC1_fusion__morph_cor_plot.png"), width = 5, height = 4,  units = "in", res = 300 )
plot(df_morph_count_PCA__plot)
x <- dev.off()


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