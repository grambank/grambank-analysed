
cat("#### Spatial & Phylo Model ####\n")
source("spatiophylogenetic_modelling/analysis/featurewise/INLA_featurewise_set_up.R")
sink(file = file.path(  OUTPUTDIR , "INLA_featurewise_dual_log.txt"), split = T)

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
                            control.predictor = list(compute=TRUE, link=1), #@Sam should we do this?
                            control.family = list(control.link=list(model="logit")),  #@Sam should we do this?
                            data = grambank_df, family = "binomial") })
  
  if (class(output) != "try-error") {
    suppressWarnings(    saveRDS(output, file = paste0(OUTPUTDIR, "dual_process_rdata/spatial_phylo_", feature, ".rdata")))
    #Don't be alarmed by the suppress warnings. saveRDS() is being kind and reminding us that the package stats may not be available when loading. However, this is not a necessary warning for us so we've wrapped saveRDS in suppressWarnings
    
  }
  
  rm(output)  
}


sink(file = NULL)