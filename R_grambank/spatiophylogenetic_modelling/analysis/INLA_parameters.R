#for INLA
kappa = 2 # smoothness parameter as recommended by Dinnage et al. (2020)
sigma = c(1, 1.15) # Sigma parameter. First value is not used. 

#for testing different spatial settings
kappa_vec = c(2, 2, 2.5)
sigma_vec =  list(c(1, 1.15), c(1, 2), c(1, 3))

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

# we use a penalising complexity prior which are particular suited to the
# analyses of additive models. 
# We should test the sensitivity of priors on the full data model
pcprior = list(prec = list(
  prior="pc.prec",
  param = c(1, 0.1)) # This prior suggests that the probability that variance for the random effect is greater than 1 is 10%
)

# We need to fix the residual variance to one, since it is not an identifiable quantity
# within a binomial model. 
obs_hyper <- list(prec = list(initial = log(1), fixed = TRUE))

#making list of lgs in tree
tree_fn <- "output/spatiophylogenetic_modelling/processed_data/EDGE_pruned_tree.tree"
if(!(file.exists(tree_fn))){
  source("spatiophylogenetic_modelling/processing/pruning_EDGE_tree.R")}

tree_tips_df <-  read.tree(tree_fn) %>% 
  .$tip.label %>% 
  as.data.frame() %>% 
  rename("Language_ID"= ".")

gb_df_cropped <- readr::read_tsv(file = "output/GB_wide/GB_cropped_for_missing.tsv",show_col_types = F) 

lgs_in_analysis <- gb_df_cropped %>% 
    dplyr::select(Language_ID) %>% 
  inner_join(tree_tips_df, by = "Language_ID")

