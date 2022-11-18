source("requirements.R")

model_surprisal <- qs::qread("output/unusualness/model_surprisal.qs")

model_summary <- summary(model_surprisal)
bayes_R2(model_surprisal)

Coefficient <- c("Intercept", "SD", "SD (phylogeny)", "SD (spatial)")

Estimate <- c(model_summary$fixed$Estimate, 
              model_summary$spec_pars$Estimate, 
              model_summary$random$Language_ID$Estimate, 
              model_summary$random$Language_ID2$Estimate)


`Estimated error` <- c(model_summary$fixed$Est.Error,
                       model_summary$spec_pars$Est.Error, 
                       model_summary$random$Language_ID$Est.Error, 
                       model_summary$random$Language_ID2$Est.Error)

brms_table_unusualness <- cbind(Coefficient, Estimate, `Estimated error`)

brms_table_unusualness %>% 
  as.data.frame() %>% 
  mutate(Estimate = as.numeric(Estimate) %>% round(2)) %>% 
  mutate(`Estimated error` = as.numeric(`Estimated error`) %>% round(2)) %>% 
  write_tsv("output/unusualness/tables/unsualness_brms_predict_table.tsv", na = "")
  

gb <- read_tsv("output/unusualness/tables/model_df.tsv", na = "")

# Plot this
ggplot(gb,aes(x=Res_Surprisal))+geom_histogram(bins = 30)

language_meta <- read_tsv("output/non_GB_datasets/glottolog_AUTOTYP_areas.tsv", show_col_types = F) %>% 
  dplyr::select(Language_ID, AUTOTYP_area)

gb <- gb %>% 
  left_join(language_meta)

gb$AUTOTYP_area <- fct_reorder(gb$AUTOTYP_area, gb$Surprisal)

gb %>% 
  ggplot() +
  geom_violin(aes(x = Surprisal, y = AUTOTYP_area),fill = "turquoise3", color = "turquoise3", alpha = 0.3, size = 0.6, draw_quantiles = c(0.25, 0.5, 0.75))  +
  geom_violin(aes(x = Pred_Surprisal, y = AUTOTYP_area), color = "black", alpha = 0, fill = "turquoise3")  +
  theme_classic() +
  theme(text = element_text(size = 18), 
        axis.title.y = element_blank()) +
  xlab("Unusualness")

ggsave("output/unusualness/plots/violin_unusualness_modeL_predictions.png", height = 10, width = 7)
