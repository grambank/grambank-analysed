source("requirements.R")
p_load(beepr)

fns <- list.files("output/spatiophylogenetic_modelling/results_debug_2022-04-08/phylo_only/", pattern = "*.rdata", full.names = T)


df<- data.frame(matrix(ncol = 3, nrow = 0))
colnames(df) <- c("Feature_ID", "Log precision for phy_id_generic", "Log precision for phy_id_iid_model") 
df$Feature_ID <- as.character(df$Feature_ID)
df$`Log precision for phy_id_generic` <- as.numeric(df$`Log precision for phy_id_generic`)
df$`Log precision for phy_id_iid_model` <- as.numeric(df$`Log precision for phy_id_iid_model`)




for(fn in fns) {

#  fn <- fns[1]
  cat("I'm on ", fn, ".\n")
    output <- readRDS(fn)
    
    feature <- fn %>% str_extract("GB[0-9]*[a|b]?")
    
df_feature <-    output$mode$theta %>% 
      t() %>% 
      as.data.frame() %>% 
      mutate(Feature_ID = feature) 
    
df <- df %>% 
  full_join(df_feature, by = c("Feature_ID", "Log precision for phy_id_generic", "Log precision for phy_id_iid_model"))
    
}
beep(1)

df %>% 
  write_tsv("output/spatiophylogenetic_modelling/fails/theta_values_phylo_only_tsv")

df %>% 
  ggplot() +
  geom_point(aes(x = `Log precision for phy_id_generic`, y = `Log precision for phy_id_iid_model`))

ggsave("output/spatiophylogenetic_modelling/fails/theta_values_phylo_only_plot.png")
