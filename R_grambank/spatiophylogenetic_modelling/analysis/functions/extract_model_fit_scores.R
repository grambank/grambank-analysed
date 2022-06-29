#this script is to be applied to the output of INLA_all_models.R
#This script is written by Sam Passmore, with minor adjustments by Hedvig Skirgård, and Russell Dinnage.

#The output could be either based on the real data or the simulated. Set the variable sim_or_real to indicate which is relevant. You can specify this from CLI if you want. If you don't specify, the default is real.

source("requirements.R")

# load precision matrices and calculate normalising constants
prec_mats <- readRDS("output/spatiophylogenetic_modelling/processed_data/precision_matrices.RDS")
phylo_det <- determinant(prec_mats$phylogenetic_precision, logarithm = TRUE)
spatial_det <- determinant(prec_mats$spatial_precision, logarithm = TRUE)

phylo_const <- 0.5 * phylo_det$modulus * phylo_det$sign
spatial_const <- 0.5 * spatial_det$modulus * spatial_det$sign

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
df <- data.frame(matrix(ncol = 40, nrow = 0))
colnames(df) <- c( "phylogeny_only_waic"   ,          "spatial_only_waic"      ,         "AUTOTYP_area_waic",
                              "dual_model_waic"  ,               "trial_model_waic"     ,           "phylogeny_only_pit"   ,           "spatial_only_pit"            ,
                               "AUTOTYP_area_pit"         ,       "dual_model_pit"             ,     "trial_model_pit"       ,         "phylogeny_only_mlik_integration",
                              "spatial_only_mlik_integration"   ,"AUTOTYP_area_mlik_integration" ,  "dual_model_mlik_integration"    , "trial_model_mlik_integration"   ,
                              "phylogeny_only_mlik_gaussian" ,  "spatial_only_mlik_gaussian"  ,    "AUTOTYP_area_mlik_gaussian"    ,  "dual_model_mlik_gaussian"    ,
                             "trial_model_mlik_gaussian"    ,   "phylogeny_only_cpo"    ,          "spatial_only_cpo"   ,             "AUTOTYP_area_cpo"    ,
                           "dual_model_cpo"          ,        "trial_model_cpo"            ,   "phylogeny_only_bf_integration",
                         "spatial_only_bf_integration", "AUTOTYP_area_bf_integration", "dual_model_bf_integration", "trial_model_bf_integration",
                         "phylogeny_only_bf_gaussian", "spatial_only_bf_gaussian", "AUTOTYP_area_bf_gaussian", "dual_model_bf_gaussian",
                         "trial_model_bf_gaussian"   ,  "phylogeny_only_dic"    ,          "spatial_only_dic"   ,             "AUTOTYP_area_dic"    ,           
                   "dual_model_dic"          ,        "trial_model_dic"          )

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
    
    phylogeny_only_dic = qs[[1]]$dic$dic,
    spatial_only_dic = qs[[2]]$dic$dic,
    AUTOTYP_area_dic = qs[[3]]$dic$dic,
    dual_model_dic  = qs[[4]]$dic$dic,
    trial_model_dic  = qs[[5]]$dic$dic,

    ## I don't think this summary means much. My understanding is that it
    ## is the distribution of PIT values that matter (that is, they should have a
    ## uniform distribution for a well-fitting model). This isn't easy to summarise
    ## really
    phylogeny_only_pit = mean(qs[[1]]$pit ,na.rm = T),
    spatial_only_pit = mean(qs[[2]]$pit ,na.rm = T),
    AUTOTYP_area_pit = mean(qs[[3]]$pit ,na.rm = T),
    dual_model_pit  = mean(qs[[4]]$pit,na.rm = T),
    trial_model_pit  = mean(qs[[5]]$pit , na.rm = T) ,

    phylogeny_only_mlik_integration = mean(qs[[1]]$mlik[1] ,na.rm = T) + phylo_const,
    spatial_only_mlik_integration = mean(qs[[2]]$mlik[1] ,na.rm = T) + spatial_const,
    AUTOTYP_area_mlik_integration = mean(qs[[3]]$mlik[1] ,na.rm = T),
    dual_model_mlik_integration  = mean(qs[[4]]$mlik[1] ,na.rm = T) + phylo_const + spatial_const,
    trial_model_mlik_integration  = mean(qs[[5]]$mlik[1] ,na.rm = T) + phylo_const + spatial_const,

    phylogeny_only_mlik_gaussian = mean(qs[[1]]$mlik[2] ,na.rm = T) + phylo_const,
    spatial_only_mlik_gaussian = mean(qs[[2]]$mlik[2] ,na.rm = T) + spatial_const,
    AUTOTYP_area_mlik_gaussian = mean(qs[[3]]$mlik[2] ,na.rm = T),
    dual_model_mlik_gaussian  = mean(qs[[4]]$mlik[2] ,na.rm = T) + phylo_const + spatial_const,
    trial_model_mlik_gaussian  = mean(qs[[5]]$mlik[2] ,na.rm = T) + phylo_const + spatial_const,

    ## these should be summarised as the sum of the log values
    phylogeny_only_cpo = sum(log(qs[[1]]$cpo), na.rm = T),
    spatial_only_cpo = sum(log(qs[[2]]$cpo), na.rm = T),
    AUTOTYP_area_cpo = sum(log(qs[[3]]$cpo), na.rm = T),
    dual_model_cpo  = sum(log(qs[[4]]$cpo), na.rm = T),
    trial_model_cpo  = sum(log(qs[[5]]$cpo), na.rm = T))

    ## bayes factors
    best_mlik_integration <- max(df_spec$phylogeny_only_mlik_integration,
                                 df_spec$spatial_only_mlik_integration,
                                 df_spec$AUTOTYP_area_mlik_integration,
                                 df_spec$dual_model_mlik_integration,
                                 df_spec$trial_model_mlik_integration)

    best_mlik_gaussian <- max(df_spec$phylogeny_only_mlik_gaussian,
                              df_spec$spatial_only_mlik_gaussian,
                              df_spec$AUTOTYP_area_mlik_gaussian,
                              df_spec$dual_model_mlik_gaussian,
                              df_spec$trial_model_mlik_gaussian)

  df <- df_spec %>%
    full_join(df,by = c("Feature_ID", "phylogeny_only_waic", "spatial_only_waic", "AUTOTYP_area_waic", "dual_model_waic", "trial_model_waic", "phylogeny_only_dic", "spatial_only_dic",
  "AUTOTYP_area_dic", "dual_model_dic", "trial_model_dic", "phylogeny_only_pit", "spatial_only_pit", "AUTOTYP_area_pit", "dual_model_pit", "trial_model_pit",
  "phylogeny_only_mlik_integration", "spatial_only_mlik_integration", "AUTOTYP_area_mlik_integration", "dual_model_mlik_integration", "trial_model_mlik_integration",
  "phylogeny_only_mlik_gaussian", "spatial_only_mlik_gaussian", "AUTOTYP_area_mlik_gaussian", "dual_model_mlik_gaussian", "trial_model_mlik_gaussian", "phylogeny_only_cpo",
  "spatial_only_cpo", "AUTOTYP_area_cpo", "dual_model_cpo", "trial_model_cpo"))



}

df %>% 
  write_tsv("output/spatiophylogenetic_modelling/featurewise/model_scores.tsv", na = "")
