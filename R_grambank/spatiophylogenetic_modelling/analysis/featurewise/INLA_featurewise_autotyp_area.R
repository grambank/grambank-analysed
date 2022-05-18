

###

cat("#### AUTOTYP_area only model ####\n")
source("spatiophylogenetic_modelling/analysis/featurewise/INLA_featurewise_set_up.R")
sink(file = file.path(  OUTPUTDIR , "INLA_featurewise_autotyp_area_log.txt"), split = T)

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
    
    suppressWarnings(    saveRDS(output, file = paste0(OUTPUTDIR, "autotyp_area_only/autotyp_area_only_", feature, ".rdata")))
    #Don't be alarmed by the suppress warnings. saveRDS() is being kind and reminding us that the package stats may not be available when loading. However, this is not a necessary warning for us so we've wrapped saveRDS in suppressWarnings
  }
  rm(output)  
}
cat("All done with the autotyp_area only model, 100% done!")

sink(file = NULL)