## Run Script

## Signal strength parameters
lambda_strong = 0.8
lambda_medium = 0.5
lambda_weak = 0.2

CLI <- "no" #set to "Yes" if you want it
brms <- "no" #set to "Yes" if you want it

if(CLI == "Yes") {
  command0 = paste("Rscript spatiophylogenetic_modelling/simulations/set_up.R")
  command1 = paste0("RScript spatiophylogenetic_modelling/simulations/single_lambda_test.R ", lambda_strong)
  command2 = paste0("RScript spatiophylogenetic_modelling/simulations/single_geography_test.R ", lambda_strong)
  command3 = paste0("RScript spatiophylogenetic_modelling/simulations/dual_test.R ", lambda_strong)
  
  system(command0)
  system(command1)
  system(command2)
  system(command3)
  
  ## Medium Signal
  command4 = paste0("RScript spatiophylogenetic_modelling/simulations/single_lambda_test.R ", lambda_medium)
  command5 = paste0("RScript spatiophylogenetic_modelling/simulations/single_geography_test.R ", lambda_medium)
  command6 = paste0("RScript spatiophylogenetic_modelling/simulations/dual_test.R ", lambda_medium)
  
  system(command4)
  system(command5)
  system(command6)
  
  ## Weak Signal
  command7 = paste0("RScript spatiophylogenetic_modelling/simulations/single_lambda_test.R ", lambda_weak)
  command8 = paste0("RScript spatiophylogenetic_modelling/simulations/single_geography_test.R ", lambda_weak)
  command9 = paste0("RScript spatiophylogenetic_modelling/simulations/dual_test.R ", lambda_weak)
  
  system(command7)
  system(command8)
  system(command9)
  
} else {

  start_all <- Sys.time()
  
  source("spatiophylogenetic_modelling/simulations/set_up.R")
  
  sink(file = "output/spatiophylogenetic_modelling/simulation/simulations.txt", split = T)
  lambda <- lambda_strong
  
  start_segment <- Sys.time()
  source("spatiophylogenetic_modelling/simulations/single_lambda_test.R")
  source("spatiophylogenetic_modelling/simulations/single_geography_test.R")
  source("spatiophylogenetic_modelling/simulations/dual_test.R")
  end_segment <- Sys.time()
  diff <- end_segment - start_segment
  cat("Computing the analysis for lamba =", lambda, "took this amount of time:\n")
  diff

  start_segment <- Sys.time()
  lambda <- lambda_medium
  source("spatiophylogenetic_modelling/simulations/single_lambda_test.R")
  source("spatiophylogenetic_modelling/simulations/single_geography_test.R")
  source("spatiophylogenetic_modelling/simulations/dual_test.R")
  end_segment <- Sys.time()
  diff <- end_segment - start_segment
  cat("Computing the analysis for lamba =", lambda, "took this amount of time:\n")
  diff
  
  start_segment <- Sys.time()
  lambda <- lambda_weak
  source("spatiophylogenetic_modelling/simulations/single_lambda_test.R")
  source("spatiophylogenetic_modelling/simulations/single_geography_test.R")
  source("spatiophylogenetic_modelling/simulations/dual_test.R")
end_segment <- Sys.time()
  diff <- end_segment - start_segment
  cat("Computing the analysis for lamba =", lambda, "took this amount of time:\n")
  diff
  sink(file = NULL)
  
}

