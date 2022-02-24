source("requirements.R")

#script written by Sam Passmore

# Check INLA is installed
if (!is_installed("INLA")) { source(file.path("spatiophylogenetic_modelling", "install_inla.R")) } else {
  cat("Great, INLA was already installed.\n") }
suppressPackageStartupMessages(library(INLA, quietly = T, warn.conflicts = F, verbose = F))


load('spatiophylogenetic_modelling/results/jaeger_models.RData')

cat("Plotting results of sp-models.\n")

get_marginalvalues = function(values, along = seq(0, 1, by = 0.01)){
  tmarg = inla.tmarginal(function(x) 1/sqrt(x), # convert to SD 
                 values, 
                 method = "linear") 
  
  qmarg = inla.qmarginal(along, tmarg)
  qmarg
}

singlemodel_values = lapply(marginals.hyperpar[1:6], '[[', 2) 
singlemodel_values = lapply(singlemodel_values, function(x) get_marginalvalues(x, along = c(0.025, 0.5, 0.975)))


value_length = length(singlemodel_values[[1]])
plot_df = data.frame(values = unlist(singlemodel_values), 
           PC = rep(c(1:3, 1:3), each = value_length), 
           effect = rep(c("Phylogenetic", "Spatial"), each = value_length * 3),
           model = rep(c("Phylogenetic", "Spatial"), each = value_length * 3))


doublemodel_values = lapply(marginals.hyperpar[7:9], function(x){
  phy_effect = get_marginalvalues(x[[2]], along = c(0.025, 0.5, 0.975))
  sp_effect = get_marginalvalues(x[[3]], along = c(0.025, 0.5, 0.975))
  data.frame(values = c(phy_effect, sp_effect),
             PC = rep(1:3, each = length(phy_effect) * 2),
             effect = rep(c("Phylogenetic", "Spatial"), each = length(phy_effect)),
             model = "Spatiophylogenetic")
})

doublemodel_values = do.call(rbind, doublemodel_values)

plot_df = rbind(plot_df, doublemodel_values)

plot_df_summ = groupwiseMean(values ~ PC + effect + model,
                             data   = plot_df,
                             conf   = 0.95,
                             digits = 3)

plot_df_summ$PC = paste0("PC", plot_df_summ$PC)

col_vector <- c("purple4", "turquoise3")

plot_1 = ggplot(plot_df_summ, aes(y = model, x = Mean, fill = effect, col = effect)) + 
  geom_errorbar(width = 0.2, 
                aes(xmin = Mean - Trad.lower, 
                    xmax = Mean + Trad.upper),
                position = position_dodge(width = 0.4)) + 
  geom_point(shape=21, size=2, position = position_dodge(width = 0.4)) + 
  scale_colour_manual(values = col_vector) + 
  scale_fill_manual(values = col_vector) + 
    theme_classic() + 
  facet_grid(~PC) + 
  ylab("") + xlab("Standard deviations") +
  theme(legend.title = element_blank(), legend.position="bottom") 

plot_1

ggsave(filename = 'spatiophylogenetic_modelling/figures/spatiophylo_plot.tiff',
       plot_1, height = 4, width = 7)
ggsave(filename = 'spatiophylogenetic_modelling/figures/spatiophylo_plot.png',
       plot_1, height = 4, width = 7)

plot_df$PC = paste0("PC", plot_df$PC)

plot_2 = ggplot(plot_df, aes(y = model, x = values, fill = effect, col = effect)) + 
  geom_density_ridges(alpha = 0.5, bandwidth = 0.015) + 
  facet_grid(~PC) + 
  ylab("") + xlab("Estimates") +
  scale_colour_manual(values = col_vector) + 
  scale_fill_manual(values = col_vector) + 
  theme_classic() + 
    theme(legend.title = element_blank(), legend.position="bottom") 

plot_2

ggsave(filename = 'spatiophylogenetic_modelling/figures/spatiophylo_ridgeplot.tiff',
       plot_2, height = 4, width = 7)