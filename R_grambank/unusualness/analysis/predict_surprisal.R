
# Set working directory for output
#setup outpur dirs
OUTPUTDIR <- file.path("output/unusualness/")
if (!dir.exists(OUTPUTDIR)) { dir.create(OUTPUTDIR) }		

OUTPUTDIR_tables <- file.path("output/unusualness/tables")
if (!dir.exists(OUTPUTDIR_tables)) { dir.create(OUTPUTDIR_tables) }		

OUTPUTDIR_plots <- file.path("output/unusualness/plots")
if (!dir.exists(OUTPUTDIR_plots)) { dir.create(OUTPUTDIR_plots) }		

# Load pkgs
source("requirements.R")

surprisal_fn <- paste0(OUTPUTDIR_tables, "surprisal.tsv")
if(!file.exists(surprisal_fn)){
  source("unusualness/analysis/get_unusualness_bayesLCA.R")
}
gb <- read_tsv(file = surprisal_fn)

# Recode endangerment
gb<-gb %>%
  mutate(Endangerement=ifelse(aes %in% c("threatened","moribund","nearly_extinct"),"endangered",aes))

# Keep only the optimal Estimator - kernel 15
gb <- gb %>%
  filter(Estimator=="Kernel 15")

#########################################
## (5) Model unusualness in terms of genealogical, areal covariates, and endangerement status
#########################################

# Get phylogenetic covariance matrix
###################################
tree_fn <- "output/spatiophylogenetic_modelling/processed_data/EDGE_pruned_tree.tree"
if(!(file.exists(tree_fn))){
  source("spatiophylogenetic_modelling/processing/pruning_EDGE_tree.R")}
phylogenetic_tree = read.tree(tree_fn)

# Subsetting languages in GB to those with matches in the tree and subsetting the tree to those with matches in gb
phylogenetic_tree_tips_df <- phylogenetic_tree$tip.label %>% 
  as.data.frame() %>% 
  rename(Language_ID = ".")

gb <- gb %>% 
  inner_join(phylogenetic_tree_tips_df, by = "Language_ID")

phylo_covar_mat<-vcv.phylo(ape::keep.tip(phylogenetic_tree, gb$Language_ID))

# Output this
phylo_covar_mat  %>% 
  qs::qsave("output/unusualness/tables/phylo_covar_mat.qs")
###################################

# Get spatial covariance matrix
###################################

# Source functions relevant to the spatial variance-covariance matrices
source('spatiophylogenetic_modelling/analysis/varcov_spatial.R')

kappa = 2 # smoothness parameter as recommended by Dinnage et al. (2020)
sigma = c(1, 1.15) # Sigma parameter. First value is not used. 

# Get the longitude and latitude data from the simulated datasets
locations_df = read.delim('output/non_GB_datasets/glottolog-cldf_wide_df.tsv', sep = "\t") %>%
  inner_join(gb, by = "Language_ID") #subset to matches in tree and in cropped in GB.

spatial_covar_mat = varcov.spatial(locations_df[,c("Longitude", "Latitude")],
                                   cov.pars = sigma,
                                   kappa = kappa)$varcov

# Add language names to the matrix
rownames(spatial_covar_mat)<-gb$Language_ID
colnames(spatial_covar_mat)<-gb$Language_ID

spatial_covar_mat  %>% 
  qs::qsave("output/unusualness/tables/spatial_covar_mat.qs")
###################################

# Duplicates Language_ID just for the sake of the requirements in brms
gb$Language_ID2<-gb$Language_ID

# Function that obtains a predictive model of surprisal
regression_surprisal<-function(df,s_cov,p_cov) {
  m<-brm(formula = Surprisal~
           (1|gr(Language_ID, cov = p_cov))+
           (1|gr(Language_ID2, cov = s_cov)),
         data=df,
         data2=list(s_cov=s_cov,
                    p_cov=p_cov),
         chains = 4,
         iter = 6000,
         warmup = 2000,
         cores = 4,
         control = list(adapt_delta=0.99),
         backend="cmdstanr")
  return(m)}



model_surprisal<-regression_surprisal(gb,
                                      s_cov=spatial_covar_mat,
                                      p_cov=phylo_covar_mat)
# Estimate Bayesian  R2
bayes_R2(model_surprisal)

# Check summary
summary(model_surprisal_lca)
conditional_effects(model_surprisal_lca)