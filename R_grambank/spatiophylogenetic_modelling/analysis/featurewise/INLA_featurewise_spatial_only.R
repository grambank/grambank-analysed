cat("\n\n#### Spatial only Model ####\n\n")
source("spatiophylogenetic_modelling/analysis/featurewise/INLA_featurewise_set_up.R")
sink(file = file.path(  OUTPUTDIR , "INLA_featurewise_spatial_only_log.txt"), split = T)

#make empty df to bind to
df_spatial_only <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(df_spatial_only) <- c("2.5%","50%", "97.5%", "Feature_ID", "effect", "model", "waic", "marginals.hyperpar.spatial_id_iid_model", "marginals.hyperpar.spatial_id_generic") 
df_spatial_only$`2.5%` <- as.numeric(df_spatial_only$`2.5%`)
df_spatial_only$`50%` <- as.numeric(df_spatial_only$`50%`)
df_spatial_only$`97.5%` <- as.numeric(df_spatial_only$`97.5%`)
df_spatial_only$Feature_ID <- as.character(df_spatial_only$Feature_ID)
df_spatial_only$model <- as.character(df_spatial_only$model)
df_spatial_only$effect <- as.character(df_spatial_only$effect)
df_spatial_only$waic <- as.numeric(df_spatial_only$waic)
df_spatial_only$marginals.hyperpar.spatial_id_iid_model <- as.list(df_spatial_only$marginals.hyperpar.spatial_id_iid_model)
df_spatial_only$marginals.hyperpar.spatial_id_generic <- as.list(df_spatial_only$marginals.hyperpar.spatial_id_generic)



index <- 0

cat("Starting INLA spatial-only featurewise runs at", as.character(Sys.time()), ".\n")

cat("#### spatial only model ####\n")
for(feature in features){
  
  #feature <- features[1]
  
  cat(paste0("# Running the spatial-only model on feature ", feature, ". That means I'm ", round(index/length(features) * 100, 2), "% done.\n The time is ", as.character(Sys.time()), ".\n"))
  index <- index + 1 
  
  
  formula <- eval(substitute(this_feature ~
                               f((spatial_id_generic), 
                                 model = "generic0",
                                 Cmatrix = spatial_prec_mat,
                                 constr = TRUE, 
                                 hyper = pcprior) + 
                               f(spatial_id_iid_model,
                                 model = "iid", 
                                 hyper = pcprior), list(this_feature=as.name(feature))))
  
  output <- try({INLA::inla(formula = formula,
                            control.compute = list(waic=TRUE, dic = FALSE, mlik = FALSE, config = TRUE),
                            control.inla = list(tolerance = 1e-6, h = 0.001),
                            control.predictor = list(compute=TRUE, link=1), #@Sam should we do this?
                            control.family = list(control.link=list(model="logit")),   #@Sam should we do this?
                            data = grambank_df,family = "binomial") })
  
  if (class(output) != "try-error") {
    
    #pulling out phy_id_generic effect
    #if the hessian has negative eigenvalues, then the hyperpar will contain inf values and the extract won't work, therefore there's an if statement testing for this.
    
    spatial_effect_generic = try(expr = {inla.tmarginal(function(x) 1/sqrt(x),
                                                      output$marginals.hyperpar$`Precision for spatial_id_generic`,
                                                      method = "linear") %>%
        inla.qmarginal(c(0.025, 0.5, 0.975), .)}
    )
    
    if (class(spatial_effect_generic) != "try-error") {
      df_spatial_only_generic  <- spatial_effect_generic %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
        mutate(Feature_ID = feature) %>% 
        mutate(effect = "spatial_only_generic") %>% 
        mutate(model = "spatial_only") %>% 
        mutate(waic = output$waic$waic)  %>% 
        mutate(marginals.hyperpar.spatial_id_generic = output$marginals.hyperpar[1])
    } else{
      
      
      cat(paste0("Couldn't extract spatial generic effect from feature ", feature, ", making empty df!\n"))
      
      df_spatial_only_generic <- tibble(
        "2.5%" = c(NA),
        "50%" =c(NA),
        "97.5%" =c(NA)) %>% 
        mutate(Feature_ID = feature) %>% 
        mutate(effect = "spatial_only_generic") %>% 
        mutate(model = "spatial_only") %>% 
        mutate(waic = output$waic$waic)  %>% 
        mutate(marginals.hyperpar.spatial_id_generic = output$marginals.hyperpar[1])
    }
    
    #pulling out spatial_id_iid_model effect
    #if the hessian has negative eigenvalues, then the hyperpar will contain inf values and the extract won't work, therefore there's an if statement testing for this.
    
    spatial_effect_iid_model = try(expr = {
      inla.tmarginal(function(x) 1/sqrt(x),
                     output$marginals.hyperpar$`Precision for spatial_id_iid_model`,
                     method = "linear") %>%
        inla.qmarginal(c(0.025, 0.5, 0.975), .) }
    )
    
    if (class(spatial_effect_iid_model) != "try-error") {
      df_spatial_only_iid_model <- spatial_effect_iid_model %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
        mutate(Feature_ID = feature) %>% 
        mutate(effect = "spatial_only_iid_model") %>% 
        mutate(model = "spatial_only") %>% 
        mutate(waic = output$waic$waic)  %>% 
        mutate(marginals.hyperpar.spatial_id_iid_model = output$marginals.hyperpar[2]) } else{
          
          cat(paste0("Couldn't extract spatial iid effect from feature ", feature, ", making empty df!\n") )
          
          df_spatial_only_iid_model<- tibble(
            "2.5%" = c(NA),
            "50%" =c(NA),
            "97.5%" =c(NA)) %>%  
            mutate(Feature_ID = feature) %>% 
            mutate(effect = "spatial_only_iid_model") %>% 
            mutate(model = "spatial_only") %>% 
            mutate(waic = output$waic$waic)  %>% 
            mutate(marginals.hyperpar.spatial_id_iid_model = output$marginals.hyperpar[2])
        }
    
    df_spatial_only <- df_spatial_only  %>% 
      full_join(df_spatial_only_iid_model, 
                by = c(join_columns, "marginals.hyperpar.spatial_id_iid_model")) %>% 
      full_join(df_spatial_only_generic, 
                by = c(join_columns, "marginals.hyperpar.spatial_id_generic"))
    
  
if(save_RDS_featurewise ==1){
      
qs::qsave(x = output, preset = "high", file = paste0(OUTPUTDIR, "spatial_only/spatial_only_", feature, ".rdata")) 
}
    }
  rm(output)  
}

df_spatial_only %>% write_tsv(file = file.path(OUTPUTDIR, "spatial_only/df_spatial_only.tsv"))
df_spatial_only %>% saveRDS(file = file.path(OUTPUTDIR, "spatial_only/df_spatial_only.Rdata"))


cat("All done with the spatial only model, 100% done!\n The time is ", as.character(Sys.time()), ".\n")

sink(file = NULL)