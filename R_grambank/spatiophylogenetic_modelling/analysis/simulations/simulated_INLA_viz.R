source("fun_def_h_load.R")
h_load(pkg = c("qs", "unglue", "tidyverse"))

#written by Russell Dinnage

col_vector <- c("purple4", "turquoise3")

OUTPUTDIR <- "output/spatiophylogenetic_modelling/simulation_plots/"
if(!dir.exists(OUTPUTDIR)){
  dir.create(OUTPUTDIR, showWarnings = FALSE)
}

qs_files <- list.files("output/spatiophylogenetic_modelling/simulated_output", pattern = ".qs", full.names = TRUE)
params <- list.files("output/spatiophylogenetic_modelling/simulated_output", pattern = ".qs") %>%
  unglue_data("Prop{prop}_Lambda{lamb}Iter{iter}.qs",
              convert = TRUE)

summary_mean <- qs_files %>%
  purrr::map(qread) %>% 
  map_dfr(~tibble(only_phy = mean(.x[[1]]$icc_posterior),
                  only_spat = mean(.x[[2]]$icc_posterior),
                  both_spat = mean(.x[[3]]$icc_posterior[ , 1]),
                  both_phy = mean(.x[[3]]$icc_posterior[ , 2])))

all_mean <- bind_cols(params, summary_mean)


plot_dat_mean <- all_mean %>%
  select(prop, lamb, only_phy, only_spat, both_phy, both_spat) %>%
  pivot_longer(cols = c(-prop, -lamb), names_to = "label", values_to = "value") %>%
  separate(label, c("model", "var"), "_") %>%
  mutate(lamb = as.factor(lamb),
         model = factor(model, levels = c("only", "both")),
         Variable = ifelse(var == "phy", "Phylogenetic", "Spatial"),
         )

plot_dat_mean = plot_dat_mean %>% 
  mutate(model2 = 
            ifelse(model == "only", var,
                   ifelse(model == "both", model, NA)))

plot_dat_mean$model2 = recode(plot_dat_mean$model2, 
         "2" = "Dual Process",
         phy = "Phylogeny Only",
         spat = "Spatial Only")

p1 = ggplot(plot_dat_mean, aes(x = lamb, y = value)) +
  geom_boxplot(aes(fill = Variable), lwd = 0.2) +
  geom_point(aes(fill = Variable), alpha = 0.25, pch = 21,
             colour = "black",
             position = position_dodge(width = 0.75)) +
  facet_grid(prop ~ model2, scales = "free") +
  ylab("Estimated Posterior Mean of 'Lambda'") +
  xlab("Simulated Value of Phylogenetic 'Lambda'") +
  scale_y_continuous(sec.axis = dup_axis(name = "Proportion of trait presence", labels = NULL)) + 
  theme_minimal(base_size = 8, base_line_size = 0.2) + 
  theme(legend.position = "bottom", legend.title = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.2)) 
 
plot(p1)

ggsave(paste0(OUTPUTDIR, "sim_results_means_with_points_SP.pdf"), 
       width = 120, units = "mm")

##### save pdf ########
pdf(file = paste0(OUTPUTDIR, "sim_results_means_with_points.pdf"))

ggplot(plot_dat_mean, aes(lamb, value)) +
  geom_boxplot(aes(fill = Variable)) +
  geom_point(aes(fill = Variable), alpha = 0.25, pch = 21,
             colour = "black",
             position = position_dodge(width = 1)) +
  facet_grid(prop ~ model, scales = "free",
             labeller = label_both) +
  ylab("Estimated Posterior Mean of 'Lambda'") +
  xlab("Simulated Value of Phylogenetic 'Lambda'") +
  theme_minimal()

x <- dev.off()

pdf(file = paste0(OUTPUTDIR,"sim_results_means_boxplots_only.pdf"))

ggplot(plot_dat_mean, aes(lamb, value)) +
  geom_boxplot(aes(fill = Variable)) +
  scale_fill_manual(values = col_vector) +
  facet_grid(prop ~ model, scales = "free",
             labeller = label_both) +
  ylab("Estimated Posterior Mean of 'Lambda'") +
  xlab("Simulated Value of Phylogenetic 'Lambda'") +
  theme_minimal()

x <- dev.off()