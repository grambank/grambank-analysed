#for INLA
kappa = 2 # smoothness parameter as recommended by Dinnage et al. (2020)
sigma = c(1, 1.15) # Sigma parameter. First value is not used. 

#for waic_test.R and spatial_parameters.R
kappa_vec = c(2, 4, 1, 2, 2)
sigma_vec =  c(1.15, 2, 40, 10, 20)

#### Set up model priors ####

## Taken from Dinnage et al. (2020):
### As in the manuscript, we use a “PC” prior, which stand for “Penalizing Complexity”.
### This is a standard prior developed by the developers of INLA, which is “weakly informative”.
### It places slightly more of the prior probability density on values close to zero,
### but has a “long tail”, which allows the data to push the parameter away from zero if
### there is good evidence (i.e. the likelihood of the data is higher).
### The PC prior has two parameters, p1 and p2: p2 is the proportion of the prior probability
### density that falls above values greater than p1.
###
### Given we are modelling (somewhat) Gaussian data that we have standardised to a variance of 1,
### we would not expect any random factor to have a variance greater than one.
### So we will set our prior to only have about 10% of its prior probability density
### above 1. We will guestimate the total variance possibly explained by the
### phylogenetic effect based on it’s diagonal entries, which are all equal
### to 2.7373648. So to get to 1, the scaling factor would have to be about 0.36

pcprior = list(prec =list(prior="pc.prec", param = c(1, 0.1)))