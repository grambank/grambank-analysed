## Single INLA test
## Binomial INLA test
source("requirements.R")

## Parameters
if(CLI == "Yes") {
  args = commandArgs(trailingOnly=TRUE)
  lambda = as.numeric(args[1])
}

#source("spatiophylogenetic_modelling/simulations/set_up.R")

model_data = data.frame(longitude = longitude,
                        latitude = latitude,
                        phy_id_int = 1:nrow(grambank_metadata),
                        phy_id2_int = 1:nrow(grambank_metadata),
                        spat_id_int =  1:nrow(grambank_metadata),
                        spat_id2_int = 1:nrow(grambank_metadata),
                        glottocodes = grambank_metadata$Language_ID,
                        glottocodes2 = grambank_metadata$Language_ID)



cat("Simulation for Lambda only model with Lambda =", lambda, "...\n")

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
      lambda, "\n."
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
  lambda_model = inla(formula = y ~
                           f(phy_id_int,
                             model = "generic0",
                             Cmatrix = phylo_prec_mat,
                             constr = TRUE,
                             hyper = pcprior_phy) +
                           f(phy_id2_int,
                             model = "iid",
                             hyper = pcprior_phy,
                             constr = TRUE),
                         family = "binomial",
                         control.compute = list(waic=TRUE),
                         control.inla =
                           list(tolerance = 1e-6, h = 0.001),
               control.mode(theta = c(2.02, 1.819)),
                         data = model_data)
  
  if(brms != "no"){
    print("brms...")
    brms_model <- brm(
      y ~ 1 + (1|gr(glottocodes, cov = A)) + (1|glottocodes2), 
      data = model_data, 
      family = bernoulli(), 
      data2 = list(A = phylo_covar_mat),
      prior = c(
        prior(normal(0, 50), "Intercept"),
        prior(student_t(3, 0, 20), "sd")
      ),
      chains = 1
      ) 
  }
  
  print("Phylo D...")
  phylo_d_results = phylo.d(data = model_data,
                              names.col = glottocodes, 
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

saveRDS(output_list, file = 
          paste0(
              OUTPUTDIR, "/lambda_",
            lambda,
            "_simulation.RDS"))
