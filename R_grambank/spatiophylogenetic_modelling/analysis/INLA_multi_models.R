# This script tests the following:
# Dual model (spatial + phylo)
# Trial (spatial + phylo + AUTOTYP-area)

source("requirements.R")

source("set_random_seed.R")

beep <- 0
if(beep == 1){
  h_load("beepr")
  beep(11) #starting sound
}

#wrangle CLI input
args = commandArgs(trailingOnly = TRUE)

if(length(args) != 0){
  sim_or_real <- args[1]
  start = args[2]
  end <- args[3]
  range <- start:end
  prec_matrices <- args[4]
  pcprior_choice <- args[5]
} else { #if you're running this script chunkwise in Rstudio or similar instead of via command line, you'll read in the parameters this way:
  sim_or_real <- "real"
  start <- 1
  end <- 113
  range <- start:end
  prec_matrices <- NULL
  pcprior_choice <- "default"
}

if(sim_or_real == "sim"){
  df_fn <- "output/spatiophylogenetic_modelling/simulated_data/simulated_data_df.tsv"
  if(!file.exists(df_fn)){
    cat(paste0("Simulating data.\n"))
    source("spatiophylogenetic_modelling/analysis/simulations/simulate_data.R")
  }
  df <- readr::read_tsv(file = df_fn, show_col_types = F)
  OUTPUTDIR_top <- "output/spatiophylogenetic_modelling/simulated_output/"
}

if(sim_or_real == "real"){
  df_fn <- "output/GB_wide/GB_cropped_for_missing.tsv"
  if(!file.exists(df_fn)){
    cat(paste0("Making GB data wide, binarising, cropping etc..\n"))
    source("make_wide.R")
    source("make_wide_binarized.R")
    source("impute_missing_values.R")
  }  
  df <- readr::read_tsv(file =   df_fn,show_col_types = F)
  OUTPUTDIR <- "output/spatiophylogenetic_modelling/featurewise/"
}

#dir setup
#make output dirs
if (!dir.exists("output/spatiophylogenetic_modelling/")) {
  dir.create("output/spatiophylogenetic_modelling/")
}

if (!dir.exists(OUTPUTDIR)) {
  dir.create(OUTPUTDIR)
}

#load parameters
source("global_variables.R")
source("spatiophylogenetic_modelling/analysis/INLA_parameters.R")

#load functions
source('spatiophylogenetic_modelling/analysis/functions/varcov_spatial.R')
source('spatiophylogenetic_modelling/analysis/functions/strip_inla.R')
source("spatiophylogenetic_modelling/install_inla.R")

time <- as.character(Sys.time())

sink(file = paste0(OUTPUTDIR , "INLA_log_range_", start, "_", end,"_", time, ".txt" ), split = T)

#set up features to go to go through
#checking which ones have an qs file already and excluding them

#loading inputs

## See the script make_precisionmatrices.R for details on how these were created
## Both matrices have their variance scaled to 1. 

cat("\n###\nLoading covariance matrices...\n")

if(is.null(prec_matrices)){
precision_matrices_fn <- "output/spatiophylogenetic_modelling/processed_data/precision_matrices_kappa_2_sigma_1.15.RDS"
if(!(file.exists(precision_matrices_fn))){
  source("spatiophylogenetic_modelling/analysis/make_precisionmatrices.R")} 
} else{
  precision_matrices_fn <- paste0("output/spatiophylogenetic_modelling/processed_data/", prec_matrices)
  if(!(file.exists(precision_matrices_fn))){
    source("spatiophylogenetic_modelling/analysis/make_precisionmatrices.R")} 
  }

precision_matrices = readRDS(precision_matrices_fn)
phylo_prec_mat = precision_matrices$phylogenetic_precision
spatial_prec_mat = precision_matrices$spatial_precision

cat("\n###\nDone with covariance matrices.\n")

#if a pcprior has not been specified from CLI, go with the default 
if(pcprior_choice == "default"){
  pcprior_vec <- c(prior_ten_percent)
  cat(paste0("I'm using just one pc prior. It is: \n", 
             pcprior_vec[[1]][2], "\n"
  ))
  } 

if(pcprior_choice == "loop_all_priors"){
  cat(paste0("I'm looping through all the suggested pc priors. These are: \n", 
            pcprior_vec[[1]][2], "\n",
            pcprior_vec[[2]][2], "\n",
            pcprior_vec[[3]][2], "\n",
            pcprior_vec[[4]][2], "\n"
      ))
  }

#reading in AUTOTYP-area
if (!file.exists("output/non_GB_datasets/glottolog_AUTOTYP_areas.tsv")) { 
  source("unusualness/processing/assigning_AUTOTYP_areas.R") }		
autotyp_area <- read.delim("output/non_GB_datasets/glottolog_AUTOTYP_areas.tsv", sep = "\t") %>%
  dplyr::select(Language_ID, AUTOTYP_area_id_iid_model = AUTOTYP_area)

#subsetting data to only entries which are also in the tree in case other slunk along
data <- lgs_in_analysis %>% 
  inner_join(df, by = "Language_ID") %>% 
  left_join(autotyp_area, by = "Language_ID")

x <- assert_that(all(data$Language_ID == lgs_in_analysis$Language_ID), msg = "Data doesn't match!")

## Since we are using a sparse phylogenetic matrix, we need to math taxa to the correct
## rows in the matrix
data$phylo_id = match(data$Language_ID, rownames(phylo_prec_mat))
data$spatial_id = match(data$Language_ID, rownames(spatial_prec_mat))
data$obs_id = 1:nrow(data)

#features to loop over
features <- df %>% 
  dplyr::select(-Language_ID) %>% 
  colnames() 

features <- features[range]


for(n in 1:length(pcprior_vec)){
#  n <- 1
  pcprior <- pcprior_vec[[n]]$param
   cat(paste0("###\nI'm running the analysis with pcprior = ", pcprior[2], ".\n###\n"))
   
index <- 0

for(feature in features){
  #  feature <- features[1]
  cat(paste0("I'm on ", feature, " and the time is ", Sys.time(), ".\n"))
  cat(paste0("Precision matrix = ", basename(precision_matrices_fn), ".\n"))
  cat(paste0("pcprior = ", pcprior[2], ".\n"))
  
  
  saved_file =   paste0(OUTPUTDIR, feature,"_",   substr(x = basename(precision_matrices_fn), 20, 37), "_pcprior",   pcprior[2], ".qs")
  
  if(file.exists(saved_file)){
    cat(paste0("The analysis file already exists!\n", 
               saved_file, "\n
               Moving on to next item.\n
               "))
    
  }else{
          
  # #### Dual Model ####
  # 
  # # The dual model contains the phylogenetic and spatial precision models using the same
  # # priors and scaled precision matrices. As well as a single effect for residuals (constrained to 1).
  
  formula <-  eval(substitute(expr = this_feature ~
                                f(spatial_id,
                                  model = "generic0",
                                  Cmatrix = spatial_prec_mat,
                                  hyper = pcprior) +
                                f(phylo_id,
                                  model = "generic0",
                                  Cmatrix = phylo_prec_mat,
                                  hyper = pcprior)  +
                                f(obs_id, model = "iid", hyper = obs_hyper),
                              list(this_feature=as.name(feature))))
  
  dual_model<- INLA::inla(formula = formula,
                          control.compute = list(waic = TRUE, cpo = TRUE, dic = TRUE), 
                          family = "binomial",
                          control.predictor=list(link=1),
                          data = data)
  
  cat(paste0("Finished running dual model on ", feature, " and the time is ", Sys.time(), ".\n"))
  
  
  
  # #### Trial Model ####
  # 
  # # The trial model contains the phylogenetic and spatial precision models using the same
  # # priors and scaled precision matrices. As well as a single effect for residuals (constrained to 1). In addition, there is also the categorical AUTOTYP-area included
  
  formula <-  eval(substitute(expr = this_feature ~
                                f(spatial_id,
                                  model = "generic0",
                                  Cmatrix = spatial_prec_mat,
                                  hyper = pcprior) +
                                f(phylo_id,
                                  model = "generic0",
                                  Cmatrix = phylo_prec_mat,
                                  hyper = pcprior)  +
                                f(obs_id, model = "iid", hyper = obs_hyper) +
                                f(AUTOTYP_area_id_iid_model, 
                                  hyper = pcprior,
                                  model = "iid"),
                              list(this_feature=as.name(feature))))
  
  
  trial_model<- INLA::inla(formula = formula,
                           control.compute = list(waic = TRUE, cpo = TRUE, dic = TRUE), 
                           family = "binomial",
                           control.predictor=list(link=1),
                           data = data)
  
  cat(paste0("Finished running trial on ", feature, " and the time is ", Sys.time(), ".\n"))
  
############
  #######strip interesting information from the INLA objects
  #########
  
  ## strip_inla function is located in another script. 
  #For each model we strip out the information we need to calculate heritability scores
  ## to save on storage space. 
  

  cat(paste0("Starting running strip inla on dual  of ", feature, " and the time is ", Sys.time(), ".\n"))
  
  dual_model_stripped = try(expr = {strip_inla(dual_model)})
  
  if (class(dual_model_stripped) == "try-error") {
    dual_model_stripped <- NULL
  }
  
  cat(paste0("Finished running strip inla on dual  of ", feature, " and the time is ", Sys.time(), ".\n"))

    cat(paste0("Starting running strip inla on trial model of  ", feature, " and the time is ", Sys.time(), ".\n"))
  
  trial_model_stripped <-try(expr = {strip_inla(trial_model)})
  
  if (class(trial_model_stripped) == "try-error") {
    trial_model_stripped <- NULL
  }
  
  cat(paste0("Finished running strip inla on trial model of  ", feature, " and the time is ", Sys.time(), ".\n"))
  
  
  
  #### Save output ####
  model_outputs = list(dual_model_stripped, 
                       trial_model_stripped)
  
    qs::qsave(model_outputs, 
            file = saved_file) }
  
  if(beep == 1){
    beep()
  }
  
  index <- index + 1
  
  cat(paste0("I've finished ", feature, " and the time is ", Sys.time(), ".\n That means I'm ", round((index/length(range)), 2)*100, "% done.\n"))
  rm(model_outputs, dual_model, trial_model)
} }


sink(file = NULL)