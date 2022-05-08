## Single INLA test
## Binomial INLA test
source("requirements.R")

#source("spatiophylogenetic_modelling/simulations/set_up.R")

## Parameters
if(CLI == "Yes") {
  args = commandArgs(trailingOnly=TRUE)
  lambda = as.numeric(args[1])
}

cat("Simulation for a dual process model with Lambda =", lambda, "...\n")

model_data = data.frame(longitude = longitude,
                        latitude = latitude,
                        phy_id_int = 1:nrow(grambank_metadata),
                        spat_id_int =  1:nrow(grambank_metadata),
                        error_int = 1:nrow(grambank_metadata),
                        glottocodes1 = grambank_metadata$Language_ID,
                        glottocodes2 = grambank_metadata$Language_ID,
                        glottocodes3 = grambank_metadata$Language_ID)



# rates matrix
q = matrix(c(-0.5, 0.5, 0.5, -0.5), 2)

output_list = list()
iter = 20
for(i in 1:iter){
  
  cat("I'm on iteration", 
      i, 
      "out of", 
      iter, 
      ". This is with lambda =", 
      lambda, "and dual process.\n"
  )
  
  y = geiger::sim.char(geiger::rescale(tree,
                                       lambda,
                                       model = "lambda"), 
                       q, 
                       model="discrete")[,1,]
  
  model_data$y = as.numeric(y) - 1
  
  print("fitDiscrete...")
  pagels_lambda = fitDiscrete(tree, 
                              factor(y), 
                              transform = "lambda")
  
  print("INLA...")
  inla_model = inla(y ~ 
         f(phy_id_int,
           model = "generic0",
           Cmatrix = phylo_prec_mat,
           constr = TRUE,
           hyper = pcprior_phy) + 
         f(spat_id_int,
           model = "generic0",
           Cmatrix = spatial_prec_mat,
           constr = TRUE,
           hyper = pcprior_phy) +
         f(error_int,
           model = "iid",
           constr = TRUE,
           hyper = pcprior_phy) , 
       data = model_data)
  
  if(brms != "no"){
    print("brms...")
    brms_model <- brm(
      y ~ 1 + 
        (1|gr(glottocodes1, cov = spatial)) + 
        (1|gr(glottocodes2, cov = phylogeny)) + 
        (1|glottocodes3), 
      data = model_data, 
      family = bernoulli(), 
      data2 = list(spatial = spatial_covar_mat,
                   phylogeny = phylo_covar_mat),
      prior = c(
        prior(normal(0, 50), "Intercept"),
        prior(student_t(3, 0, 20), "sd")
      ),
      chains = 1
    )
  }
  
  print("Phylo D...")
  phylo_d_results = phylo.d(data = model_data,
                              names.col = glottocodes1, 
                              phy = tree,
                              binvar = y)
  
  if(brms != "no"){
    output_list[[i]] = list(y = y,
                            pagels_lambda = pagels_lambda,
                            inla_model = lambda_model,
                            brms_model = brms_model,
                            phylo_d = phylo_d_results)
  }else{
    output_list[[i]] = list(y = y,
                            pagels_lambda = pagels_lambda,
                            inla_model = lambda_model,
                            phylo_d = phylo_d_results)
  }
}


suppressWarnings(
saveRDS(output_list, file = 
          paste0(OUTPUTDIR, "dual_",
            lambda,
            "_simulation.RDS")) )
