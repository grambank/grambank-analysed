

###Spatial Only


cat("#### Spatial only Model ####\n")
source("spatiophylogenetic_modelling/analysis/featurewise/INLA_featurewise_set_up.R")
sink(file = file.path(  OUTPUTDIR , "INLA_featurewise_spatial_only_log.txt"), split = T)

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
  
if(save_RDS_featurewise ==1){
      
    suppressWarnings(  saveRDS(output, file = paste0(OUTPUTDIR, "spatial_only/spatial_only_", feature, ".rdata")) )
    #Don't be alarmed by the suppress warnings. saveRDS() is being kind and reminding us that the package stats may not be available when loading. However, this is not a necessary warning for us so we've wrapped saveRDS in suppressWarnings.
  }}
  rm(output)  
}

cat("All done with the spatial only model, 100% done!")

sink(file = NULL)