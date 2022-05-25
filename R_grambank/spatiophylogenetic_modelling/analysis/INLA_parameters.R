#for INLA
kappa = 2 # smoothness parameter as recommended by Dinnage et al. (2020)
sigma = c(1, 1.15) # Sigma parameter. First value is not used. 

#for waic_test.R and spatial_parameters.R
kappa_vec = c(2, 4, 1, 2, 2)
sigma_vec =  c(1.15, 2, 40, 10, 20)
