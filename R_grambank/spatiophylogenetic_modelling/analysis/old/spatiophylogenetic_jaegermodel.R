source('requirements.R')

#script written by Sam Passmore

# load variational covariance matrix taken from geoR::varcov_spatial
source('spatiophylogenetic_modelling/analysis/varcov_spatial.R')

# Check INLA is installed
if (!is_installed("INLA")) { cat("INLA wasn't installed, installing now.\n") 
  source(file.path("spatiophylogenetic_modelling", "install_inla.R")) } else {
  cat("Great, INLA was already installed, loading now.\n") }
suppressPackageStartupMessages(
  library(INLA, quietly = T, warn.conflicts = F, verbose = F)
)

OUTPUTDIR_figures <- "output/spatiophylogenetic_modelling/figures"
if (!dir.exists(OUTPUTDIR_figures)){dir.create(OUTPUTDIR_figures)}

OUTPUTDIR_results <- "output/spatiophylogenetic_modelling/results"
if (!dir.exists(OUTPUTDIR_results)){dir.create(OUTPUTDIR_results)}

cat("#### Building Jaeger tree models ####\n")

#### Inputs ####
# language metadata
# pca
pca_filename = 'output/PCA/PCA_language_values.tsv'
pca_components = read_tsv(pca_filename, col_types = cols()) %>% 
  mutate(PC1 = scale(PC1), PC2 = scale(PC2), PC3 = scale(PC3)) 

languages <- read_csv(GRAMBANK_LANGUAGES, col_types=LANGUAGES_COLSPEC) %>%		
  dplyr::select(Language_ID = Language_level_ID, Family_name, Name, Longitude, Latitude, Macroarea) %>% 
  distinct(Language_ID, .keep_all = T) %>% 
  inner_join(pca_components, by = "Language_ID")

# Check that there are as many languages as we expect
x <- assert_that(nrow(pca_components) == n_imputed, msg = "The total number of languages has changed. This might mean the data doesn't align with the phylogeny")

# trees
tree_filename = 'output/spatiophylogenetic_modelling/processed_data/jaeger_pruned.tree'
phylogenetic_tree = read.tree(tree_filename)

# Subset PCA and languages to Jaeger set
pca_components = pca_components[pca_components$Language_ID %in% phylogenetic_tree$tip.label,]
languages = languages[languages$Language_ID %in% pca_components$Language_ID,]

# prune tree to dataset
taxa = pca_components$Language_ID
phylogenetic_tree = keep.tip(phylogenetic_tree, tip = taxa)

# Check the tree has the number of tips and data we expect. 
x <- assert_that(nrow(pca_components) == n_overlap_imputed_and_jaeger_tree, msg = "The number of languages in PCA have changed. It is no longer compatible with the phylogeny")
x <- assert_that(nrow(languages) == n_overlap_imputed_and_jaeger_tree, msg = "The number of languages in PCA have changed. It is no longer compatible with the phylogeny")

#check that there are as many observations as we expect
x <- assert_that(Ntip(phylogenetic_tree) == n_overlap_imputed_and_jaeger_tree, msg = "The number of tips in the phylogeny has changed. It will not be compatible with the dataset")

# All tips are in data
x <- assert_that(all(phylogenetic_tree$tip.label %in% pca_components$Language_ID), msg = "The data and phylogeny taxa do not match")

#### Parameters ####
pca_basename = tools::file_path_sans_ext(pca_filename) %>% basename()
tree_basename = tools::file_path_sans_ext(tree_filename) %>% basename()
out_name = paste(pca_basename, tree_basename, sep = "_")

#### Spatial Jittering ####
## There are some number of languages that have identical spatial coordinates, which we cannot allow for the spatial analysis.
## I have jittered the coordinates that are identical
## Jittering moves locations randomly by about 0.3 and 1 degree in Longitude & Latitude
duplicate_coords = languages[duplicated(languages[,c("Longitude", "Latitude")]) | duplicated(languages[,c("Longitude", "Latitude")], fromLast = TRUE),"Language_ID"]
duplicate_rowid = languages$Language_ID %in% duplicate_coords
languages$Latitude[duplicate_rowid] = jitter(languages$Latitude[duplicate_rowid], factor = 1)
languages$Longitude[duplicate_rowid] = jitter(languages$Longitude[duplicate_rowid], factor = 1)

#### Visualise repsonse ####
pdf('output/spatiophylogenetic_modelling/figures/PCA_normalitychecks.pdf', width =8.3, height = 11.7)
p1 = ggplot(data=pca_components) +
  geom_histogram(aes(x=PC1, color=I("black"),fill=I("orchid")), bins = 30) + 
  ggtitle("Histogram")

p2 = ggplot(data=pca_components, aes(sample = PC1)) + 
  stat_qq() + 
  stat_qq_line() + 
  ggtitle("QQ-Norm")

p1 + p2 + plot_annotation(title = paste('PC1', pca_basename, tree_basename, sep = ": "))

p3 = ggplot(data=pca_components) +
  geom_histogram(aes(x=PC2, color=I("black"),fill=I("orchid")), bins = 30) + 
  ggtitle("Histogram")

p4 = ggplot(data=pca_components, aes(sample = PC2)) + 
  stat_qq() + 
  stat_qq_line() + 
  ggtitle("QQ-Norm")

p3 + p4 + plot_annotation(title = paste('PC2', pca_basename, tree_basename, sep = ": "))

p5 = ggplot(data=pca_components) +
  geom_histogram(aes(x=PC3, color=I("black"),fill=I("orchid")), bins = 30) +
  ggtitle("Histogram")

p6 = ggplot(data=pca_components, aes(sample = PC3)) + stat_qq() + stat_qq_line() + ggtitle("QQ-Norm")

p5 + p6 + plot_annotation(title = paste('PC3', pca_basename, tree_basename, sep = ": "))

x <- dev.off()


#### Phylogenetic covariance matrix ####

cat("Calculating the phylogenetic variance covariance matrix.\n")

phylo_covar_mat <- ape::vcv(phylogenetic_tree)
phylo_covar_mat <- phylo_covar_mat / max(phylo_covar_mat)
# The diagonal of phylo_covar_mat should inform our prior
phylo_covar_mat <- phylo_covar_mat / exp(determinant(phylo_covar_mat)$modulus[1] /
                                          nrow(phylo_covar_mat))
phylo_prec_mat <-solve(phylo_covar_mat)

# Phylogenetic matrix is right dims
x <- assert_that(all(dim(phylo_prec_mat) == c(n_overlap_imputed_and_jaeger_tree, n_overlap_imputed_and_jaeger_tree)), msg = "The phylogeny has changed and will not match the data")

#### Spatial covariance matrix ####

cat("Calculating the spatial variance covariance matrix.\n")
## Ensure the order of languages matches the order within the phylogeny
languages = languages[order(match(languages$Language_ID, rownames(phylo_prec_mat))),]

spatial_covar_mat = varcov.spatial(languages[,c("Longitude", "Latitude")], cov.pars = sigma, kappa = kappa)$varcov
dimnames(spatial_covar_mat) = list(languages$Language_ID, languages$Language_ID)

cov2precision = function(spatial_covar_mat){
  spatial_covar_mat = spatial_covar_mat / exp(determinant(spatial_covar_mat)$modulus[1] /
                                                nrow(spatial_covar_mat))
  spatial_prec_mat = solve(spatial_covar_mat)
  spatial_prec_mat
}

spatial_prec_mat = cov2precision(spatial_covar_mat)

# Do Phylo and Spatial rownames match?
if(all(rownames(phylo_prec_mat) == rownames(spatial_covar_mat))){
  spatial_covar_mat = spatial_covar_mat / exp(determinant(spatial_covar_mat)$modulus[1] /
                                                nrow(spatial_covar_mat))
  spatial_prec_mat = solve(spatial_covar_mat)
  #spatial_prec_mat[1:5, 1:5]
} else {
  stop("Spatial and Phylo matrices do not align")
}

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

pcprior_phy = list(prec =list(prior="pc.prec", param = c(1, 0.1)))
pcprior_spa = list(prec =list(prior="pc.prec", param = c(1, 0.1)))

## Adding random effect ids
grambank_pca = pca_components %>%
  left_join(tibble(Language_ID = rownames(phylo_prec_mat),
                   phy_id = 1:nrow(phylo_prec_mat),
                   sp_id = 1:nrow(spatial_prec_mat)), by = "Language_ID")


cat("#### Phylogenetic only model ####\n")
cat("Building Phylogenetic models: PC1\n")
# PC1
phylogenetic_PC1 = inla(PC1 ~
                            f(phy_id, model = "generic0", Cmatrix = phylo_prec_mat,
                              constr = TRUE, hyper = pcprior_phy),
                        control.compute = list(waic=TRUE, dic = TRUE, mlik = FALSE, config = TRUE),
                        control.predictor = list(compute = TRUE),
                          data = grambank_pca)

# PC2
cat("Building Phylogenetic models: PC2\n")
phylogenetic_PC2 = inla(PC2 ~
                          f(phy_id, model = "generic0", Cmatrix = phylo_prec_mat,
                            constr = TRUE, hyper = pcprior_phy),
                        control.compute = list(waic=TRUE),
                        data = grambank_pca)

# PC3
cat("Building Phylogenetic models: PC3\n")
phylogenetic_PC3 = inla(PC3 ~
                          f(phy_id, model = "generic0", Cmatrix = phylo_prec_mat,
                            constr = TRUE, hyper = pcprior_phy),
                        control.compute = list(waic=TRUE),
                        data = grambank_pca)

## Phylogenetic effect
phylo_effect_varPC1 = inla.tmarginal(function(x) 1/sqrt(x),
                                  phylogenetic_PC1$marginals.hyperpar$`Precision for phy_id`,
                                  method = "linear") %>%
  inla.qmarginal(c(0.025, 0.5, 0.975), .)

# inla.tmarginal(function(x) 1/sqrt(exp(x)), phylogenetic_PC1$internal.marginals.hyperpar[[2]])

phylo_effect_varPC2 = inla.tmarginal(function(x) 1/sqrt(x),
                                     phylogenetic_PC2$marginals.hyperpar$`Precision for phy_id`,
                                     method = "linear") %>%
  inla.qmarginal(c(0.025, 0.5, 0.975), .)


phylo_effect_varPC3 = inla.tmarginal(function(x) 1/sqrt(x),
                                     phylogenetic_PC3$marginals.hyperpar$`Precision for phy_id`,
                                     method = "linear") %>%
  inla.qmarginal(c(0.025, 0.5, 0.975), .)

phylo_effect_var = rbind(phylo_effect_varPC1, phylo_effect_varPC2, phylo_effect_varPC3)
dimnames(phylo_effect_var) = list(c("PC1", "PC2", "PC3"), c("2.5%", "50%", "97.5%"))

write.csv(phylo_effect_var, file.path("output", "spatiophylogenetic_modelling", "results",
                                      paste0(out_name, "_phylogenyOnly.csv")))

cat("#### Spatial only Model ####\n")
# PC1
cat("Building Spatial models: PC1\n")
spatial_PC1 = inla(PC1 ~
                          f(sp_id, model = "generic0", Cmatrix = spatial_prec_mat,
                            constr = TRUE, hyper = pcprior_spa),
                        control.compute = list(waic=TRUE, dic = TRUE),
                        data = grambank_pca)

# PC2
cat("Building Spatial models: PC2\n")
spatial_PC2 = inla(PC2 ~
                          f(sp_id, model = "generic0", Cmatrix = spatial_prec_mat,
                            constr = TRUE, hyper = pcprior_spa),
                        control.compute = list(waic=TRUE),
                        data = grambank_pca)

# PC3
cat("Building Spatial models: PC3\n")
spatial_PC3 = inla(PC3 ~
                          f(sp_id, model = "generic0", Cmatrix = spatial_prec_mat,
                            constr = TRUE, hyper = pcprior_spa),
                        control.compute = list(waic=TRUE),
                        data = grambank_pca)


## Spatial effect
spatial_effect_varPC1 = inla.tmarginal(function(x) 1/sqrt(x),
                                     spatial_PC1$marginals.hyperpar$`Precision for sp_id`,
                                     method = "linear") %>%
  inla.qmarginal(c(0.025, 0.5, 0.975), .)


spatial_effect_varPC2 = inla.tmarginal(function(x) 1/sqrt(x),
                                     spatial_PC2$marginals.hyperpar$`Precision for sp_id`,
                                     method = "linear") %>%
  inla.qmarginal(c(0.025, 0.5, 0.975), .)


spatial_effect_varPC3 = inla.tmarginal(function(x) 1/sqrt(x),
                                     spatial_PC3$marginals.hyperpar$`Precision for sp_id`,
                                     method = "linear") %>%
  inla.qmarginal(c(0.025, 0.5, 0.975), .)

spatial_effect_var = rbind(spatial_effect_varPC1, spatial_effect_varPC2, spatial_effect_varPC3)
dimnames(spatial_effect_var) = list(c("PC1", "PC2", "PC3"), c("2.5%", "50%", "97.5%"))

write.csv(spatial_effect_var, file.path("output", "spatiophylogenetic_modelling", "results",
                                      paste0(out_name, "_spaceOnly.csv")))

cat("#### Spatial & Phylo Model ####\n")
cat("Building Spatiophylogenetic models: PC1\n")
spatiophylogenetic_PC1 = inla(PC1 ~ 
                          f(phy_id, model = "generic0", Cmatrix = phylo_prec_mat, constr = TRUE, hyper = pcprior_phy) + 
                          f(sp_id, model = "generic0", Cmatrix = spatial_prec_mat, constr = TRUE, hyper = pcprior_spa), 
                        control.compute = list(waic=TRUE),
                        data = grambank_pca)

cat("Building Spatiophylogenetic models: PC2\n")
spatiophylogenetic_PC2 = inla(PC2 ~  
                                f(phy_id, model = "generic0", Cmatrix = phylo_prec_mat, constr = TRUE, hyper = pcprior_phy) + 
                                f(sp_id, model = "generic0", Cmatrix = spatial_prec_mat, constr = TRUE, hyper = pcprior_spa), 
                              control.compute = list(waic=TRUE),
                              data = grambank_pca)

cat("Building Spatiophylogenetic models: PC3\n")
spatiophylogenetic_PC3 = inla(PC3 ~  
                                f(phy_id, model = "generic0", Cmatrix = phylo_prec_mat, constr = TRUE, hyper = pcprior_phy) + 
                                f(sp_id, model = "generic0", Cmatrix = spatial_prec_mat, constr = TRUE, hyper = pcprior_spa), 
                              control.compute = list(waic=TRUE),
                              data = grambank_pca)

#### Interpreting Random effects ####
#### From Dinnage et al. (2020)
#### If you are used to thinking in terms of variance and covariance (as I am), 
#### we now have to think in terms of precision and precision matrices.
#### Luckily, there is a fairly straightforward transformation to convert a scaling 
#### factor φ for a precision matrixinto a scaling factor on the equivalent covariance 
#### matrix"σ= 1/φ. We can transform the marginal posterior distribution of any 
#### parameter using theINLAfunctioninla.tmarginal(), then we can resummarise 
#### (usinginla.qmarginal()) to credible intervals to get something more 
#### interpretable if you are more familiar with variance thinking.

# PC1
spatiophylogenetic_speffect_varPC1 = inla.tmarginal(function(x) 1/sqrt(x), 
                                       spatiophylogenetic_PC1$marginals.hyperpar$`Precision for sp_id`, 
                                       method = "linear") %>%
  inla.qmarginal(c(0.025, 0.5, 0.975), .)

spatiophylogenetic_phyeffect_varPC1 = inla.tmarginal(function(x) 1/sqrt(x), 
                                                     spatiophylogenetic_PC1$marginals.hyperpar$`Precision for phy_id`, 
                                                     method = "linear") %>%
  inla.qmarginal(c(0.025, 0.5, 0.975), .)


# inla.tmarginal(function(x) 1/sqrt(exp(x)),result$internal.marginals.hyperpar[[2]])

# PC2
spatiophylogenetic_speffect_varPC2 = inla.tmarginal(function(x) 1/sqrt(x), 
                                                  spatiophylogenetic_PC2$marginals.hyperpar$`Precision for sp_id`, 
                                                  method = "linear") %>%
  inla.qmarginal(c(0.025, 0.5, 0.975), .)

spatiophylogenetic_phyeffect_varPC2 = inla.tmarginal(function(x) 1/sqrt(x), 
                                                   spatiophylogenetic_PC2$marginals.hyperpar$`Precision for phy_id`, 
                                                   method = "linear") %>%
  inla.qmarginal(c(0.025, 0.5, 0.975), .)

# PC3
spatiophylogenetic_speffect_varPC3 = inla.tmarginal(function(x) 1/sqrt(x), 
                                                  spatiophylogenetic_PC3$marginals.hyperpar$`Precision for sp_id`, 
                                                  method = "linear") %>%
  inla.qmarginal(c(0.025, 0.5, 0.975), .)

spatiophylogenetic_phyeffect_varPC3 = inla.tmarginal(function(x) 1/sqrt(x), 
                                                   spatiophylogenetic_PC3$marginals.hyperpar$`Precision for phy_id`, 
                                                   method = "linear") %>%
  inla.qmarginal(c(0.025, 0.5, 0.975), .)

spatiophylogenetic_models = rbind(
  c(spatiophylogenetic_phyeffect_varPC1, spatiophylogenetic_speffect_varPC1),
  c(spatiophylogenetic_phyeffect_varPC2, spatiophylogenetic_speffect_varPC2),
  c(spatiophylogenetic_phyeffect_varPC3, spatiophylogenetic_speffect_varPC3))

## output table

out_table = matrix(NA, ncol = 4, nrow = 3)
out_table[,1] = apply(phylo_effect_var, 1, function(x) paste(round(x, 2), collapse = ", "))
out_table[,2] = apply(spatial_effect_var, 1, function(x) paste(round(x, 2), collapse = ", "))
out_table[,3] = apply(spatiophylogenetic_models[,1:3], 1, function(x) paste(round(x, 2), collapse = ", "))
out_table[,4] = apply(spatiophylogenetic_models[,4:6], 1, function(x) paste(round(x, 2), collapse = ", "))

out_table = apply(out_table, 1:2, function(x) paste0("(", x, ")"))
dimnames(out_table) = list(c("PC1", "PC2", "PC3"),
                           c("Phylogenetic_only", "Spatial_Only", "Phylogenetic_Joint", "Spatial_Joint"))

write.csv(out_table, file.path("output", "spatiophylogenetic_modelling", "results", "jaeger_table.csv"))

#### Table 1 ####
table1 = matrix(NA, ncol = 4, nrow = 3)
dimnames(table1) = list(c("PC1", "PC2", "PC3"),
                        c("Standarddeviation_Phy",
                          "Percentage_Phy",
                          "Standarddevation_Sp",
                          "Percentage_SP"))

table1[,"Standarddeviation_Phy"] = 
  round(spatiophylogenetic_models[,2], 2)
table1[,"Standarddevation_Sp"] = 
  round(spatiophylogenetic_models[,5], 2)

### Get percentage explained results
get_percent = function(mean, sd){
  pct = pnorm(mean + sd / 2, lower.tail=TRUE) - 
    pnorm(mean - sd / 2, lower.tail=TRUE)
  round(pct, 2)
}

phy_percent = get_percent(
  colMeans(grambank_pca[,c("PC1", "PC2", "PC3")]),
  spatiophylogenetic_models[,2])

sp_percent = get_percent(
  colMeans(grambank_pca[,c("PC1", "PC2", "PC3")]),
           spatiophylogenetic_models[,5])

table1[,"Percentage_Phy"] = phy_percent
table1[,"Percentage_SP"]  = sp_percent

write.csv(table1, 
          "output/spatiophylogenetic_modelling/results/table_1.csv")

#### WAIC output ####
waic_output = matrix(NA, ncol = 3, nrow = 3)

waic_output[1,] = c(phylogenetic_PC1$waic$waic, 
                    spatial_PC1$waic$waic,
                    spatiophylogenetic_PC1$waic$waic)

waic_output[2,] = c(phylogenetic_PC2$waic$waic, 
                    spatial_PC2$waic$waic,
                    spatiophylogenetic_PC2$waic$waic)

waic_output[3,] = c(phylogenetic_PC3$waic$waic, 
                    spatial_PC3$waic$waic,
                    spatiophylogenetic_PC3$waic$waic)

dimnames(waic_output) = list(c("PC1", "PC2", "PC3"), c("Phylogeny", "Space", "Spatiophylogenetic"))
write.csv(waic_output, "output/spatiophylogenetic_modelling/results/jaeger_waic.csv")

#### Save models #### 
model_list = list(phylogenetic_PC1, phylogenetic_PC2, phylogenetic_PC3,
                  spatial_PC1, spatial_PC2, spatial_PC3, 
                  spatiophylogenetic_PC1, spatiophylogenetic_PC2, spatiophylogenetic_PC3)


marginals.hyperpar = lapply(model_list, function(x) x$marginals.hyperpar)

names(marginals.hyperpar) = c("Phy_PC1", "Phy_PC2", "Phy_PC3", 
                      "Space_PC1", "Space_PC2", "Space_PC3",
                      "SP_PC1", "SP_PC2", "SP_PC3")

save(marginals.hyperpar, 
     file = "output/spatiophylogenetic_modelling/results/jaeger_models.RData")
