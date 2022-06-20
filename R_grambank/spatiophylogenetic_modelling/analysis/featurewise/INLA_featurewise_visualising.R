source("requirements.R")

fns <- list.files("output/spatiophylogenetic_modelling/featurewise/", pattern = "*qs", full.names = T)


  df <- data.frame(matrix(ncol = 13, nrow = 0))
colnames(df) <- c("phylogeny_only_effect", "spatial_only_effect", "AUTOTYP_area_effect", "dual_model_spatial_effect", "dual_model_phylo_effect",
                  "trial_model_phylo_effect", "trial_model_spatial_effect", "trial_model_autotyp_area_effect", "phylogeny_only_waic", "spatial_only_waic",
                  "AUTOTYP_area_waic", "dual_model_waic", "trial_model_waic")

df <- df %>% 
  mutate_all(as.numeric)

df$Feature_ID <- as.character(df$phylogeny_only) 

for(fn in fns) {
  qs <- qs::qread(fn) 
#  qs <- qs::qread( fns[10]) 
  
  Feature_ID <- basename(fn) %>% str_replace_all(".qs", "")
  
  #options(warn=-1) #everytime one of the output from strip_INLA is null there will be a warning. So we're turning off warnings temporarily to not clog the console
  phylogeny_only_effect <- qs[[1]]$icc_posterior %>% mean()
  spatial_only_effect <- qs[[2]]$icc_posterior %>% mean()
  AUTOTYP_area_effect <- qs[[3]]$icc_posterior %>% mean()
  dual_model_spatial_effect  <- qs[[4]]$icc_posterior[ , 1] %>% mean()
  dual_model_phylo_effect  <- qs[[4]]$icc_posterior[ , 2] %>% mean()
  trial_model_phylo_effect  <- qs[[5]]$icc_posterior[ , 1] %>% suppressWarnings(mean())
  trial_model_spatial_effect  <- qs[[5]]$icc_posterior[ , 2] %>% suppressWarnings(mean())
  trial_model_autotyp_area_effect  <- qs[[5]]$icc_posterior[ , 4] %>% suppressWarnings(mean())
  options(warn=0)

  phylogeny_only_waic <- qs[[1]]$waic$waic
  spatial_only_waic <- qs[[2]]$waic$waic
  AUTOTYP_area_waic <- qs[[3]]$waic$waic
  dual_model_waic  <- qs[[4]]$waic$waic
  trial_model_waic  <- qs[[5]]$waic$waic
  
  df_spec <- cbind(    phylogeny_only_effect = phylogeny_only_effect,
    spatial_only_effect = spatial_only_effect,
    AUTOTYP_area_effect = AUTOTYP_area_effect,
    dual_model_spatial_effect  = dual_model_spatial_effect,
    dual_model_phylo_effect  = dual_model_phylo_effect ,
    trial_model_phylo_effect  = trial_model_phylo_effect,
    trial_model_spatial_effect  = trial_model_spatial_effect  ,
    trial_model_autotyp_area_effect  =  trial_model_autotyp_area_effect, 
    phylogeny_only_waic = phylogeny_only_waic ,
    spatial_only_waic =  spatial_only_waic ,
    AUTOTYP_area_waic =    AUTOTYP_area_waic ,
    dual_model_waic  = dual_model_waic,
    trial_model_waic  =  ifelse(is.null(trial_model_waic), NA, trial_model_waic)
  ) %>% as.data.frame()
 
  df <- df_spec %>% 
    mutate_all(as.numeric) %>% 
    mutate(Feature_ID = Feature_ID) %>% 
    full_join(df, by = c("phylogeny_only_effect", "spatial_only_effect", "AUTOTYP_area_effect", "dual_model_spatial_effect", "dual_model_phylo_effect",
                         "phylogeny_only_waic", "spatial_only_waic", "AUTOTYP_area_waic", "dual_model_waic", "trial_model_waic", "Feature_ID"))
}

df %>% 
  dplyr::select(Feature_ID, 
                phylogeny_only_waic,
                spatial_only_waic ,
                AUTOTYP_area_waic ,
                dual_model_waic  ,
                trial_model_waic ) %>% 
  reshape2::melt(id.vars= "Feature_ID") %>%
  ggplot() +
  geom_boxplot(aes(x = variable, y = value)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle= 70, hjust = 1), legend.position = "None") +
  geom_jitter(aes(x = variable, y = value, color = variable))

ggsave("output/spatiophylogenetic_modelling/featurewise/boxplot_waic.png")
  

mean(df$phylogeny_only_waic)
mean(df$spatial_only_waic)
mean(df$AUTOTYP_area_waic)
mean(df$dual_model_waic)
mean(df$trial_model_waic)

  
