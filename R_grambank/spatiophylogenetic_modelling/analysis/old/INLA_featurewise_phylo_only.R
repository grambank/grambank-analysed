cat("\n\n#### Phylogenetic only model ####\n\n")

source("spatiophylogenetic_modelling/analysis/featurewise/INLA_featurewise_set_up.R")
sink(file = file.path(  OUTPUTDIR , "phylo_only/INLA_featurewise_phylo_only_log.txt"), split = T)

#make empty df to bind to
df_phylo_only <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(df_phylo_only) <- c("2.5%","50%", "97.5%", "Feature_ID", "effect", "model", "waic", "marginals.hyperpar.phy_id_iid_model", "marginals.hyperpar.phy_id_generic") 
df_phylo_only$`2.5%` <- as.numeric(df_phylo_only$`2.5%`)
df_phylo_only$`50%` <- as.numeric(df_phylo_only$`50%`)
df_phylo_only$`97.5%` <- as.numeric(df_phylo_only$`97.5%`)
df_phylo_only$Feature_ID <- as.character(df_phylo_only$Feature_ID)
df_phylo_only$model <- as.character(df_phylo_only$model)
df_phylo_only$effect <- as.character(df_phylo_only$effect)
df_phylo_only$waic <- as.numeric(df_phylo_only$waic)
df_phylo_only$marginals.hyperpar.phy_id_iid_model <- as.list(df_phylo_only$marginals.hyperpar.phy_id_iid_model)
df_phylo_only$marginals.hyperpar.phy_id_generic <- as.list(df_phylo_only$marginals.hyperpar.phy_id_generic)

index <- 0

cat("Starting INLA phylo-only featurewise runs at", as.character(Sys.time()), ".\n")

for(feature in features){
  
  #feature <- features[1]
  
  cat(paste0("# Running the phylo-only model on feature ", 
             feature, 
             ". That means I'm ", 
             round(index/length(features) * 100, 
                   2), 
             "% done.\nThe time is ", as.character(Sys.time()), ".\n"))
  
  index <- index + 1 
  
  
  formula <- eval(substitute(this_feature ~
                               f((phy_id_generic), 
                                 model = "generic0",
                                 Cmatrix = phylo_prec_mat,
                                 constr = TRUE, 
                                 hyper = pcprior) + 
                               f(phy_id_iid_model,
                                 model = "iid", 
                                 hyper = pcprior), list(this_feature=as.name(feature))))
  
  output <- try(expr = {
    INLA::inla(formula = formula,
               control.compute = list(waic=TRUE, dic = FALSE, mlik = FALSE, config = TRUE),
               control.inla = list(tolerance = 1e-6, h = 0.001),
               control.predictor = list(compute=TRUE, link=1), 
               control.family = list(control.link=list(model="logit")),   
               data = grambank_df,family = "binomial")
  })
  
  if (class(output) != "try-error") {
    #pulling out phy_id_generic effect
    #if the hessian has negative eigenvalues, then the hyperpar will contain inf values and the extract won't work, therefore there's an if statement testing for this.
    
    phylo_effect_generic = try(expr = {INLA::inla.tmarginal(function(x) 1/sqrt(x),
                                                      output$marginals.hyperpar$`Precision for phy_id_generic`,
                                                      method = "linear") %>%
        INLA::inla.qmarginal(c(0.025, 0.5, 0.975), .)}
    )
    
    if (class(phylo_effect_generic) != "try-error") {
      df_phylo_only_generic  <- phylo_effect_generic %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
        mutate(Feature_ID = feature) %>% 
        mutate(effect = "phylo_only_generic") %>% 
        mutate(model = "phylo_only") %>% 
        mutate(waic = output$waic$waic)  %>% 
        mutate(marginals.hyperpar.phy_id_generic = output$marginals.hyperpar[1])
    } else{
      
      
      cat(paste0("Couldn't extract phy generic effect from feature ", feature, ", making empty df!\n"))
      
      df_phylo_only_generic <- tibble(
        "2.5%" = c(NA),
        "50%" =c(NA),
        "97.5%" =c(NA)) %>% 
        mutate(Feature_ID = feature) %>% 
        mutate(effect = "phylo_only_generic") %>% 
        mutate(model = "phylo_only") %>% 
        mutate(waic = output$waic$waic)  %>% 
        mutate(marginals.hyperpar.phy_id_generic = output$marginals.hyperpar[1])
    }
  
    #pulling out phy_id_iid_model effect
    #if the hessian has negative eigenvalues, then the hyperpar will contain inf values and the extract won't work, therefore there's an if statement testing for this.

    phylo_effect_iid_model = try(expr = {
      INLA::inla.tmarginal(function(x) 1/sqrt(x),
                     output$marginals.hyperpar$`Precision for phy_id_iid_model`,
                     method = "linear") %>%
        INLA::inla.qmarginal(c(0.025, 0.5, 0.975), .) }
    )
    
    if (class(phylo_effect_iid_model) != "try-error") {
      df_phylo_only_iid_model <- phylo_effect_iid_model %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        rename("2.5%" = V1, "50%" = V2, "97.5%" = V3) %>% 
        mutate(Feature_ID = feature) %>% 
        mutate(effect = "phylo_only_iid_model") %>% 
        mutate(model = "phylo_only") %>% 
        mutate(waic = output$waic$waic)  %>% 
        mutate(marginals.hyperpar.phy_id_iid_model = output$marginals.hyperpar[2]) } else{
          
          cat(paste0("Couldn't extract phy iid effect from feature ", feature, ", making empty df!\n") )
          
          df_phylo_only_iid_model<- tibble(
            "2.5%" = c(NA),
            "50%" =c(NA),
            "97.5%" =c(NA)) %>%  
            mutate(Feature_ID = feature) %>% 
            mutate(effect = "phylo_only_iid_model") %>% 
            mutate(model = "phylo_only") %>% 
            mutate(waic = output$waic$waic)  %>% 
            mutate(marginals.hyperpar.phy_id_iid_model = output$marginals.hyperpar[2])
        }
    
    df_phylo_only <- df_phylo_only  %>% 
      full_join(df_phylo_only_iid_model, 
                by = c(join_columns, "marginals.hyperpar.phy_id_iid_model")) %>% 
      full_join(df_phylo_only_generic, 
                by = c(join_columns, "marginals.hyperpar.phy_id_generic"))
    
    if(save_RDS_featurewise ==1){
    
        output$.args$.parent.frame
      qs:qsave(x = output, file = paste0(OUTPUTDIR, "phylo_only/phylo_only_", feature, ".rdata"), preset = "high") 
    }
    }
  rm(output)
}

df_phylo_only %>% write_tsv(file = file.path(OUTPUTDIR, "phylo_only/df_phylo_only.tsv"))
df_phylo_only %>% saveRDS(file = file.path(OUTPUTDIR, "phylo_only/df_phylo_only.Rdata"))

cat("All done with the phylo only model, 100% done!\n The time is ", as.character(Sys.time()), ".\n")
sink(file = NULL)