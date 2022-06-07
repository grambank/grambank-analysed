cat("\n\n#### Spatial & Phylo Model ####\n\n")

source("spatiophylogenetic_modelling/analysis/featurewise/INLA_featurewise_set_up.R")
sink(file = file.path(  OUTPUTDIR , "INLA_featurewise_dual_log.txt"), split = T)

#make empty df to bind to
df_dual_spatial_phylo <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(df_dual_spatial_phylo) <- c("2.5%","50%", "97.5%", "Feature_ID", "effect", "waic", "model", "marginals.hyperpar.phy_id_iid_model", "marginals.hyperpar.phy_id_generic", "marginals.hyperpar.spatial_id_iid_model", "marginals.hyperpar.spatial_id_generic") 
df_dual_spatial_phylo$`2.5%` <- as.numeric(df_dual_spatial_phylo$`2.5%`)
df_dual_spatial_phylo$`50%` <- as.numeric(df_dual_spatial_phylo$`50%`)
df_dual_spatial_phylo$Feature_ID <- as.character(df_dual_spatial_phylo$Feature_ID)
df_dual_spatial_phylo$effect <- as.character(df_dual_spatial_phylo$effect)
df_dual_spatial_phylo$model <- as.character(df_dual_spatial_phylo$model)
df_dual_spatial_phylo$waic <- as.numeric(df_dual_spatial_phylo$waic)
df_dual_spatial_phylo$marginals.hyperpar.phy_id_iid_model <- as.list(df_dual_spatial_phylo$marginals.hyperpar.phy_id_iid_model)
df_dual_spatial_phylo$marginals.hyperpar.phy_id_generic <- as.list(df_dual_spatial_phylo$marginals.hyperpar.phy_id_generic)
df_dual_spatial_phylo$marginals.hyperpar.spatial_id_iid_model <- as.list(df_dual_spatial_phylo$marginals.hyperpar.spatial_id_iid_model)
df_dual_spatial_phylo$marginals.hyperpar.spatial_id_generic <- as.list(df_dual_spatial_phylo$marginals.hyperpar.spatial_id_generic)

index <- 0

#something was going awry with calculating the marginal effects of a specific features, GB051, so the for loop above has been split in twain: one which saves the entire output of inla() as an rdata object in a directory and one that calculates the effects and renders the same kind of df as above. This way running the for loop with inla() can still happen and we can debug the particulars after. 

cat("Starting INLA dual process (spatial and phylo) runs featurewise runs at", as.character(Sys.time()), ".\n")

for(feature in features){
  
  #feature <- features[2]
  
  cat(paste0("# Running the spatial-phylo (double-process) model on feature ", feature, ". That means I'm ", round(index/length(features) * 100, 2), "% done.\nThe time is ", as.character(Sys.time()), ".\n"))
  index <- index + 1 
  
  formula = eval(substitute(this_feature ~
                              f((phy_id_generic), 
                                model = "generic0",
                                Cmatrix = phylo_prec_mat,
                                constr = TRUE, 
                                hyper = pcprior) + 
                              f(phy_id_iid_model,
                                model = "iid", 
                                hyper = pcprior,
                                constr = TRUE) +
                              f((spatial_id_generic), 
                                model = "generic0",
                                Cmatrix = spatial_prec_mat,
                                constr = TRUE, 
                                hyper = pcprior) + 
                              f(spatial_id_iid_model,
                                model = "iid", 
                                hyper = pcprior,
                                constr = TRUE),
                            list(this_feature=as.name(feature))))
  
  output <- try({INLA::inla(formula = formula,
                            control.compute = list(waic=TRUE, dic = FALSE, mlik = FALSE, config = TRUE),
                            control.inla = list(tolerance = 1e-6, h = 0.001),
                            control.predictor = list(compute=TRUE, link=1), 
                            control.family = list(control.link=list(model="logit")),  
                            data = grambank_df, family = "binomial") })
  
  if (class(output) != "try-error") {
    
    phylo_effect_generic = try(expr = {INLA::inla.tmarginal(function(x) 1/sqrt(x),
                                          output$marginals.hyperpar$`Precision for phy_id_generic`,
                                          method = "linear") %>%
      INLA::inla.qmarginal(c(0.025, 0.5, 0.975), .)})
    
    if (class(phylo_effect_generic) != "try-error") {
      df_phylo_generic  <- phylo_effect_generic %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
        mutate(Feature_ID = feature) %>% 
        mutate(effect = "phylo_in_dual_generic") %>% 
        mutate(model = "dual") %>% 
        mutate(waic = output$waic$waic)  %>% 
        mutate(marginals.hyperpar.phy_id_generic = output$marginals.hyperpar[1])
    } else{
      
      
      cat(paste0("Couldn't extract phy generic effect from feature ", feature, ", making empty df!\n"))
      
      df_phylo_generic <- tibble(
        "2.5%" = c(NA),
        "50%" =c(NA),
        "97.5%" =c(NA)) %>% 
        mutate(Feature_ID = feature) %>% 
        mutate(effect = "phylo_in_dual_generic") %>% 
        mutate(model = "dual") %>% 
        mutate(waic = output$waic$waic)  %>% 
        mutate(marginals.hyperpar.phy_id_generic = output$marginals.hyperpar[1])
    }
    

    phylo_effect_iid_model = try(expr = {INLA::inla.tmarginal(function(x) 1/sqrt(x),
                                            output$marginals.hyperpar$`Precision for phy_id_iid_model`,
                                            method = "linear") %>%
      INLA::inla.qmarginal(c(0.025, 0.5, 0.975), .)})
    
    
    if (class(phylo_effect_iid_model) != "try-error") {

 df_phylo_iid_model  <- phylo_effect_iid_model %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
        mutate(Feature_ID = feature) %>% 
        mutate(effect = "phylo_in_dual_iid_model") %>% 
        mutate(model = "dual") %>% 
        mutate(waic = output$waic$waic)  %>% 
        mutate(marginals.hyperpar.phy_id_iid_model = output$marginals.hyperpar[2])
    } else{
      
      
      cat(paste0("Couldn't extract phy iid_model effect from feature ", feature, ", making empty df!\n"))
      
      df_phylo_iid_model <- tibble(
        "2.5%" = c(NA),
        "50%" =c(NA),
        "97.5%" =c(NA)) %>% 
        mutate(Feature_ID = feature) %>% 
        mutate(effect = "phylo_in_dual_iid_model") %>% 
        mutate(model = "dual") %>% 
        mutate(waic = output$waic$waic)  %>% 
        mutate(marginals.hyperpar.phy_id_iid_model = output$marginals.hyperpar[2])
    }
    
    spatial_effect_generic = try(expr = {INLA::inla.tmarginal(function(x) 1/sqrt(x),
                                            output$marginals.hyperpar$`Precision for spatial_id_generic`,
                                            method = "linear") %>%
      INLA::inla.qmarginal(c(0.025, 0.5, 0.975), .)})
    
    if (class(spatial_effect_generic) != "try-error") {
      df_spatial_generic  <- spatial_effect_generic %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
        mutate(Feature_ID = feature) %>% 
        mutate(effect = "spatial_in_dual_generic") %>% 
        mutate(model = "dual") %>% 
        mutate(waic = output$waic$waic)  %>% 
        mutate(marginals.hyperpar.spatial_id_generic = output$marginals.hyperpar[3])
    } else{
      
      cat(paste0("Couldn't extract phy generic effect from feature ", feature, ", making empty df!\n"))
      
      df_spatial_generic <- tibble(
        "2.5%" = c(NA),
        "50%" =c(NA),
        "97.5%" =c(NA)) %>% 
        mutate(Feature_ID = feature) %>% 
        mutate(effect = "spatial_in_dual_generic") %>% 
        mutate(model = "dual") %>% 
        mutate(waic = output$waic$waic)  %>% 
        mutate(marginals.hyperpar.spatial_id_generic = output$marginals.hyperpar[3])
    }

    
    spatial_effect_iid_model = try(expr = {INLA::inla.tmarginal(function(x) 1/sqrt(x),
                                              output$marginals.hyperpar$`Precision for spatial_id_iid_model`,
                                              method = "linear") %>%
      INLA::inla.qmarginal(c(0.025, 0.5, 0.975), .)})
    
    if (class(spatial_effect_iid_model) != "try-error") {
      df_spatial_iid_model  <- spatial_effect_iid_model %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
        mutate(Feature_ID = feature) %>% 
        mutate(effect = "spatial_in_dual_iid_model") %>% 
        mutate(model = "dual") %>% 
        mutate(waic = output$waic$waic)  %>% 
        mutate(marginals.hyperpar.spatial_id_iid_model = output$marginals.hyperpar[4])
    } else{
      
      cat(paste0("Couldn't extract phy iid_model effect from feature ", feature, ", making empty df!\n"))
      
      df_spatial_iid_model <- tibble(
        "2.5%" = c(NA),
        "50%" =c(NA),
        "97.5%" =c(NA)) %>% 
        mutate(Feature_ID = feature) %>% 
        mutate(effect = "spatial_in_dual_iid_model") %>% 
        mutate(model = "dual") %>% 
        mutate(waic = output$waic$waic)  %>% 
        mutate(marginals.hyperpar.spatial_id_iid_model = output$marginals.hyperpar[4])
    }
    
    df_dual_spatial_phylo <- df_dual_spatial_phylo %>%
      full_join(df_phylo_generic, by = c(join_columns,
                                         "marginals.hyperpar.phy_id_generic")) %>% 
      full_join(df_phylo_iid_model, by = c(join_columns,
                                           "marginals.hyperpar.phy_id_iid_model")) %>% 
      full_join(df_spatial_generic, by = c(join_columns, 
                                           "marginals.hyperpar.spatial_id_generic")) %>% 
      full_join(df_spatial_iid_model, by = c(join_columns,
                                             "marginals.hyperpar.spatial_id_iid_model"))

    if(save_RDS_featurewise ==1){
    suppressWarnings(    saveRDS(output, file = paste0(OUTPUTDIR, "dual_process_rdata/spatial_phylo_", feature, ".rdata")))
    #Don't be alarmed by the suppress warnings. saveRDS() is being kind and reminding us that the package stats may not be available when loading. However, this is not a necessary warning for us so we've wrapped saveRDS in suppressWarnings
    }
  }
  
  rm(output)  
}

df_dual_spatial_phylo %>% write_tsv(file = file.path(OUTPUTDIR, "df_spatial_phylo.tsv"))
df_dual_spatial_phylo %>% saveRDS(file = file.path(OUTPUTDIR, "df_spatial_phylo.Rdata"))

cat("All done with the dual model, 100% done!\n The time is ", as.character(Sys.time()), ".\n")

sink(file = NULL)