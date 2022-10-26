model_surprisal <- qs::qread("output/unusualness/model_surprisal.qs")

# Check summary and estimate Bayesian  R2
summary(model_surprisal)
bayes_R2(model_surprisal)

#########################################
## (6) Predictive model
#########################################

# Add predictions and residuals
surprisal_predictions<-predict(model_surprisal)
gb$Pred_Surprisal<-surprisal_predictions[,1]
gb$Pred_Error<-surprisal_predictions[,2]
gb$Res_Surprisal<-with(gb,Surprisal-Pred_Surprisal)
gb$Z_Surprisal<-with(gb,Res_Surprisal)




# Plot this
ggplot(gb,aes(x=Res_Surprisal))+geom_histogram(bins = 30)

language_meta <- read_tsv("output/non_GB_datasets/glottolog_AUTOTYP_areas.tsv", show_col_types = F) %>% 
  dplyr::select(Language_ID, AUTOTYP_area)


gb <- gb %>% 
  left_join(language_meta)

gb$AUTOTYP_area <- fct_reorder(gb$AUTOTYP_area, gb$Surprisal)

gb %>% 
  ggplot() +
  geom_density_ridges(aes(x = Surprisal, y = AUTOTYP_area), fill = "orange", bandwidth = 0.1,alpha = 0.6,
                      quantile_lines = F, quantile_fun = median, jittered_points = TRUE, point_size = 2, point_shape = 21  ,  position = position_points_jitter(height = 0))  +
  geom_density_ridges(aes(x = Pred_Surprisal, y = AUTOTYP_area), fill = "purple4", bandwidth = 0.1, alpha = 0.8,
                      quantile_lines = F, quantile_fun = median, jittered_points = TRUE, point_size = 2, point_shape = 21  ,  position = position_points_jitter(height = 0))  +
  theme_classic() +
  theme(text = element_text(size = 18), 
        axis.title.y = element_blank()) +
  xlab("Unusualness")


ggsave("output/unusualness/plots/unusualness_modeL-predictions.png")