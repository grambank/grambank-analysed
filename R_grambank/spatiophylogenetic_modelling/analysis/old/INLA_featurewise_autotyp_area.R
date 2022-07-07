cat("\n\n#### AUTOTYP_area only model ####\n\n")
source("spatiophylogenetic_modelling/analysis/featurewise/INLA_featurewise_set_up.R")
sink(file = file.path(  OUTPUTDIR , "INLA_featurewise_autotyp_area_log.txt"), split = T)


#make empty df to bind to
df_autotyp_area_only <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(df_autotyp_area_only) <- c("2.5%","50%", "97.5%", "Feature_ID", "effect", "model", "waic", "marginals.hyperpar.AUTOTYP_area_id_iid_model") 
df_autotyp_area_only$`2.5%` <- as.numeric(df_autotyp_area_only$`2.5%`)
df_autotyp_area_only$`50%` <- as.numeric(df_autotyp_area_only$`50%`)
df_autotyp_area_only$`97.5%` <- as.numeric(df_autotyp_area_only$`97.5%`)
df_autotyp_area_only$Feature_ID <- as.character(df_autotyp_area_only$Feature_ID)
df_autotyp_area_only$model <- as.character(df_autotyp_area_only$model)
df_autotyp_area_only$effect <- as.character(df_autotyp_area_only$effect)
df_autotyp_area_only$waic <- as.numeric(df_autotyp_area_only$waic)
df_autotyp_area_only$marginals.hyperpar.AUTOTYP_area_id_iid_model <- as.list(df_autotyp_area_only$marginals.hyperpar.AUTOTYP_area_id_iid_model)

index <- 0

cat("Starting INLA AUTOTYP-area only runs featurewise runs at", as.character(Sys.time()), ".\n")

for(feature in features){
  
  #feature <- features[1]
  
  cat(paste0("# Running the autotyp_area-only model on feature ", feature, ". That means I'm ", round(index/length(features) * 100, 2), "% done.\n The time is ", as.character(Sys.time()), ".\n"))
  index <- index + 1 
  
  formula = eval(substitute(this_feature ~ 
                              f(AUTOTYP_area_id_iid_model, 
                                hyper = pcprior,
                                model = "iid"),
                            list(this_feature=as.name(feature))))
  
  output <- try(expr ={INLA::inla(formula = formula,
                                  control.compute = list(waic=TRUE, dic = FALSE, mlik = FALSE, config = TRUE),
                                  control.predictor = list(compute=TRUE, link=1), #@Sam should we do this?
                                  control.family = list(control.link=list(model="logit")),   #@Sam should we do this?
                                  control.inla = list(tolerance = 1e-6, h = 0.001),
                                  data = grambank_df,family = "binomial")})
  
  if (class(output) != "try-error") {
    
    autotyp_area_effect = try(expr = {INLA::inla.tmarginal(function(x) 1/sqrt(x),
                                         output$marginals.hyperpar$`Precision for AUTOTYP_area_id_iid_model`, method = "linear") %>%
      INLA::inla.qmarginal(c(0.025, 0.5, 0.975), .)})
    
    if (class(  autotyp_area_effect) != "try-error") {
    
    df <- autotyp_area_effect %>% 
      as.data.frame() %>% 
      t() %>% 
      as.data.frame() %>% 
      rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
      mutate(Feature_ID = feature) %>% 
      mutate(effect = "autotyp_area_iid_model") %>% 
      mutate(model = "autotyp_area_only") %>% 
      mutate(waic = output$waic$waic)  %>% 
      mutate(marginals.hyperpar.AUTOTYP_area_id_iid_model = output$marginals.hyperpar[1])
    }else{
      
      cat(paste0("Couldn't extract area effect from feature ", feature, ", making empty df!\n"))
      
      df <- tibble(
        "2.5%" = c(NA),
        "50%" =c(NA),
        "97.5%" =c(NA)) %>% 
        mutate(Feature_ID = feature) %>% 
        mutate(effect = "autotyp_area_iid_model") %>% 
        mutate(model = "autotyp_area_only") %>% 
        mutate(waic = output$waic$waic)  %>% 
        mutate(marginals.hyperpar.AUTOTYP_area_id_iid_model = output$marginals.hyperpar[1])
    }

    df_autotyp_area_only <- df_autotyp_area_only  %>% 
      full_join(df, by = c(join_columns, "marginals.hyperpar.AUTOTYP_area_id_iid_model"))
    
    if(save_RDS_featurewise ==1){
    suppressWarnings(    saveRDS(output, file = paste0(OUTPUTDIR, "autotyp_area_only/autotyp_area_only_", feature, ".rdata")))
    #Don't be alarmed by the suppress warnings. saveRDS() is being kind and reminding us that the package stats may not be available when loading. However, this is not a necessary warning for us so we've wrapped saveRDS in suppressWarnings
  }}
  rm(output)  
}


df_autotyp_area_only %>% write_tsv(file = file.path(OUTPUTDIR, "autotyp_area_only/df_autotyp_area_only.tsv"))
df_autotyp_area_only %>% saveRDS(file = file.path(OUTPUTDIR, "autotyp_area_only/df_autotyp_area_only.Rdata"))

cat("All done with the autotyp_area only model, 100% done!\n The time is ", as.character(Sys.time()), ".\n")

sink(file = NULL)