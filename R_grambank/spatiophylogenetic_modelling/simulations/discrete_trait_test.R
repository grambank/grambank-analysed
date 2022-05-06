library(ape)
library(geiger)
library(ggplot2)


iterations = 100
n_tips = 1265
tree = rtree(n = n_tips)

lambda_set = seq(0, 1, by = 0.2)

# Rates matrix
q = matrix(c(-0.5, 0.5, 0.5, -0.5), 2)


plot_df = data.frame(count = numeric(),
                     function_name = character(), 
                     lambda_value = numeric(),
                     iteration = numeric())
for(i in seq_along(lambda_set)){
  lambda_value = lambda_set[i]
  cat("Starting Lambda =", lambda_value, "...\n")
  for(i in 1:iterations){
    scaled_tree = geiger::rescale(tree, lambda_value, model = "lambda")
    
    # Discrete variable in ape
    ape_discrete = ape::rTraitDisc(phy = scaled_tree, 
                                   k = 2, 
                                   states = 0:1, 
                                   model = q)
    
    # Discrete variable in geiger
    geiger_discrete = geiger::sim.char(scaled_tree, 
                                       q, 
                                       model="discrete")
    
    p_df = data.frame(count = c(sum(ape_discrete == 1), 
                                sum(geiger_discrete == 1)),
                      function_name = c("ape::rTraitDisc", "geiger::sim.char"), 
                      lambda_value = lambda_value,
                      iteration = i)
    
    plot_df = rbind(plot_df, p_df)
  }
}

# 0.51 and 0.49 are equivalent in a discrete trait so make them the same
plot_df$props = plot_df$count / n_tips
plot_df$props = ifelse(plot_df$props <= 0.5, plot_df$props,
                       1 - plot_df$props )

ggplot(plot_df, aes(x = factor(lambda_value), 
                    y = props, 
                    fill = function_name)) + 
  geom_boxplot() +
  xlab("Lambda Value") + 
  ylab("Proportion of 1's") + 
  ggtitle("Simulated proportions of discrete traits in ape and geiger simulation tools",
          paste0("N tips = ", n_tips, " for ", iterations, " iterations.")
          ) +
  theme_light() + 
  theme(legend.title = element_blank(),
        legend.position = "right")

### Add distribution of Grambank ratios

GB_imputed_filename <- file.path("output", "GB_wide", "GB_wide_imputed_binarized.tsv")
if (!file.exists(GB_imputed_filename)) { 
  source("make_wide.R")
  source("make_wide_binarized.R")
  source("impute_missing_values.R")}		
GB_imputed <- read.delim(GB_imputed_filename, sep = "\t")

GB_counts = apply(GB_imputed[,2:ncol(GB_imputed)], 2, function(x) sum(x == 1))

gb_df = data.frame(count = GB_counts,
                   function_name = "Grambank",
                   lambda_value = "Grambank",
                   iteration = 1:length(GB_counts)
)

plot_df_gb = rbind(plot_df, gb_df)

plot_df_gb$props = plot_df_gb$count / nrow(GB_imputed)
plot_df_gb$props = ifelse(plot_df_gb$props <= 0.5, plot_df_gb$props,
                          plot_df_gb$props - 0.5)

ggplot(plot_df_gb, aes(x = factor(lambda_value), 
                    y = props, 
                    fill = function_name)) + 
  geom_boxplot() +
  xlab("Lambda Value") + 
  ylab("Count of simulated 1's") + 
  ggtitle("Simulated counts for ape and geiger simulation tools",
          paste0("N tips = ", n_tips, " for ", iterations, " iterations.")
  ) +
  theme_light() + 
  theme(legend.title = element_blank(),
        legend.position = "right")
