#this script is to be applied to the output of INLA_all_models.R 
#This script is written by Sam Passmore, with minor adjustments by Hedvig Skirgård.

#The output could be either based on the real data or the simulated. Set the variable sim_or_real to indicate which is relevant. You can specify this from CLI if you want. If you don't specify, the default is real.

source("requirements.R")

args = commandArgs(trailingOnly = TRUE)
if(length(args) != 0){
  sim_or_real <- args[1] #you'd set this by going Rscript script.R "real"
} else { #if you're running this script chunkwise in Rstudio or similar instead of via command line, you'll read in the parameters this way:
  sim_or_real <- "real"
}

if(sim_or_real == "sim"){ #running the simulated data
  fns <- list.files("output/spatiophylogenetic_modelling/simulated_output/", pattern = "*.qs", full.names = T)
}

if(sim_or_real == "real"){ #running the on the real data
  fns <- list.files("output/spatiophylogenetic_modelling/featurewise/", pattern = "*.qs", full.names = T) 
}

#This script takes the output of strip_inla() and calculates the appropriate posteriors. 

#The output of strip_inla() is calculated per effect in the INLA-models. We specify five models: phylo only, spatial only, autotyp-area only, dual (spatial + phylo) and trial (phylo, spatial and autoty-area). The strip_inla() function pulls out the waic, dic, cpo and pit scores as well as the hypersample. The hypersample is what we use here to calculate the posteriors. 


for(fn in fns){
  
  #fn <- fns[1]
    qs <- qs::qread(fn)
    phylogeny_only_waic <- qs[[1]]$waic$waic
    spatial_only_waic <- qs[[2]]$waic$waic
    AUTOTYP_area_waic <- qs[[3]]$waic$waic
    dual_model_waic  <- qs[[4]]$waic$waic
    trial_model_waic  <- qs[[5]]$waic$waic

    phylogeny_only_pit <- qs[[1]]$pit  
    spatial_only_pit <- qs[[2]]$pit
    AUTOTYP_area_pit <- qs[[3]]$pit
    dual_model_pit  <- qs[[4]]$pit
    trial_model_pit  <- qs[[5]]$pit

    phylogeny_only_mlik <- qs[[1]]$mlik  
    spatial_only_mlik <- qs[[2]]$mlik
    AUTOTYP_area_mlik <- qs[[3]]$mlik
    dual_model_mlik  <- qs[[4]]$mlik
    trial_model_mlik  <- qs[[5]]$mlik
    
    phylogeny_only_cpo <- qs[[1]]$cpi  
    spatial_only_cpo <- qs[[2]]$cpi
    AUTOTYP_area_cpo <- qs[[3]]$cpi
    dual_model_cpo  <- qs[[4]]$cpi
    trial_model_cpo  <- qs[[5]]$cpi
    
    
  
  
}