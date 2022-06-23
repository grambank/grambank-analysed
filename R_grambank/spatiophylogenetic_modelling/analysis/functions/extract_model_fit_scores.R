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


#empty df to bind to in for loop
df <- data.frame(matrix(ncol = 25, nrow = 0))
colnames(df) <- c( "phylogeny_only_waic"   ,          "spatial_only_waic"      ,         "AUTOTYP_area_waic",              
                              "dual_model_waic"  ,               "trial_model_waic"     ,           "phylogeny_only_pit"   ,           "spatial_only_pit"            ,   
                               "AUTOTYP_area_pit"         ,       "dual_model_pit"             ,     "trial_model_pit"       ,         "phylogeny_only_mlik_integration",
                              "spatial_only_mlik_integration"   ,"AUTOTYP_area_mlik_integration" ,  "dual_model_mlik_integration"    , "trial_model_mlik_integration"   ,
                              "phylogeny_only_mlik_gaussian" ,  "spatial_only_mlik_gaussian"  ,    "AUTOTYP_area_mlik_gaussian"    ,  "dual_model_mlik_gaussian"    ,   
                             "trial_model_mlik_gaussian"    ,   "phylogeny_only_cpo"    ,          "spatial_only_cpo"   ,             "AUTOTYP_area_cpo"    ,           
                           "dual_model_cpo"          ,        "trial_model_cpo"            )

df <- df %>% 
  mutate_all(as.numeric)

df$Feature_ID <- as.character() 

index <- 0

for(fn in fns){
  
  #fn <- fns[49]
    qs <- qs::qread(fn)
    index <- index + 1
    
  cat(paste0("I'm on ", basename(fn), ". ", index,".\n"))
    
    df_spec <- data.frame(
    Feature_ID = basename(fn) %>% str_replace_all(".qs", ""),
    phylogeny_only_waic = qs[[1]]$waic$waic,
    spatial_only_waic = qs[[2]]$waic$waic,
    AUTOTYP_area_waic = qs[[3]]$waic$waic,
    dual_model_waic  = qs[[4]]$waic$waic,
    trial_model_waic  = qs[[5]]$waic$waic,

    phylogeny_only_pit = mean(qs[[1]]$pit ,na.rm = T),
    spatial_only_pit = mean(qs[[2]]$pit ,na.rm = T),
    AUTOTYP_area_pit = mean(qs[[3]]$pit ,na.rm = T),
    dual_model_pit  = mean(qs[[4]]$pit,na.rm = T),
    trial_model_pit  = mean(qs[[5]]$pit , na.rm = T) ,
    
    phylogeny_only_mlik_integration = mean(qs[[1]]$mlik[1] ,na.rm = T),
    spatial_only_mlik_integration = mean(qs[[2]]$mlik[1] ,na.rm = T),
    AUTOTYP_area_mlik_integration = mean(qs[[3]]$mlik[1] ,na.rm = T),
    dual_model_mlik_integration  = mean(qs[[4]]$mlik[1] ,na.rm = T),
    trial_model_mlik_integration  = mean(qs[[5]]$mlik[1] ,na.rm = T),

    phylogeny_only_mlik_gaussian = mean(qs[[1]]$mlik[2] ,na.rm = T),
    spatial_only_mlik_gaussian = mean(qs[[1]]$mlik[2] ,na.rm = T),
    AUTOTYP_area_mlik_gaussian = mean(qs[[1]]$mlik[2] ,na.rm = T),
    dual_model_mlik_gaussian  = mean(qs[[1]]$mlik[2] ,na.rm = T),
    trial_model_mlik_gaussian  = mean(qs[[1]]$mlik[2] ,na.rm = T),

    phylogeny_only_cpo = mean(qs[[1]]$cpo, na.rm = T),
    spatial_only_cpo = mean(qs[[2]]$cpo ,na.rm = T),
    AUTOTYP_area_cpo = mean(qs[[3]]$cpo,na.rm = T),
    dual_model_cpo  = mean(qs[[4]]$cpo ,na.rm = T),
    trial_model_cpo  = mean(qs[[5]]$cpo ,na.rm = T))

  df <- df_spec %>%
    full_join(df, by = c("Feature_ID", "phylogeny_only_waic", "spatial_only_waic", "AUTOTYP_area_waic", "dual_model_waic", "trial_model_waic",
                         "phylogeny_only_pit", "spatial_only_pit", "AUTOTYP_area_pit", "dual_model_pit", "trial_model_pit", "phylogeny_only_mlik_integration",
                         "spatial_only_mlik_integration", "AUTOTYP_area_mlik_integration", "dual_model_mlik_integration", "trial_model_mlik_integration",
                         "phylogeny_only_mlik_gaussian", "spatial_only_mlik_gaussian", "AUTOTYP_area_mlik_gaussian", "dual_model_mlik_gaussian",
                         "trial_model_mlik_gaussian", "phylogeny_only_cpo", "spatial_only_cpo", "AUTOTYP_area_cpo", "dual_model_cpo", "trial_model_cpo"))

}

df %>% 
  write_tsv("output/spatiophylogenetic_modelling/featurewise/model_scores.tsv", na = "")
