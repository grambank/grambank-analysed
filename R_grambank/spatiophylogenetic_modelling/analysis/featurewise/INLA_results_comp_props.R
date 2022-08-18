source("spatiophylogenetic_modelling/analysis/featurewise/INLA_results_most_prep.R")

#comparing phylo and spatial effects with prop of 1s

comp_df <- GB_df %>%
  reshape2::melt(id.vars = "Language_ID") %>%
  group_by(variable) %>%
  summarise(one_prop = mean(value, na.rm = T)) %>%
  rename(Feature_ID = variable) %>%
  full_join(dual_summary, by = "Feature_ID") %>%
  dplyr::select(Feature_ID, one_prop, mean_phylogenetic, mean_spatial)

png("output/spatiophylogenetic_modelling/INLA_spec_results_plots/splom_dual.png")

psych::pairs.panels(comp_df[,-1],
                    method = "pearson", # correlation method
                    hist.col = "#a3afd1",# "#a9d1a3","",""),
                    density = TRUE,  # show density plots
                    ellipses = F, # show correlation ellipses
                    cex.labels= 1,
                    #           smoother= T,
                    cor=T,
                    lm=T,
                    ci = T, cex.cor = 0.9,stars = T)

x <- dev.off()
