# Load pkgs
source("requirements.R")

# Set working directory for output
#setup outpur dirs
OUTPUTDIR <- file.path("output/unusualness/")
if (!dir.exists(OUTPUTDIR)) { dir.create(OUTPUTDIR) }		

OUTPUTDIR_tables <- file.path("output/unusualness/tables")
if (!dir.exists(OUTPUTDIR_tables)) { dir.create(OUTPUTDIR_tables) }		

OUTPUTDIR_plots <- file.path("output/unusualness/plots")
if (!dir.exists(OUTPUTDIR_plots)) { dir.create(OUTPUTDIR_plots) }		

surprisal_fn <- paste0(OUTPUTDIR_tables, "/surprisal.tsv")
if(!file.exists(surprisal_fn)){
  source("unusualness/analysis/get_unusualness_bayesLCA.R")
}
gb <- read.delim(file = surprisal_fn, sep = "\t") %>% 
  group_by(Language_ID, aes) %>% 
  summarise(Surprisal = mean(Surprisal, na.rm = T))

### NEXT PART REQUIRES MATRICES ETC

#########################################
## (5) Model unusualness in terms of genealogical, areal covariates, and endangerement status
#########################################

# Recode endangerment
gb<-gb %>%
  mutate(Endangerement=ifelse(aes %in% c("threatened","moribund","nearly_extinct"),"endangered",aes))

#In the spatiophylogenetic modelling of the features, we use the dataset cropped for missing data but without imputation. For the unsualness analsyis, we use the imputed data. They are different in the feature values, but it is the same subset of langauges in both. Therefore, we can use the same precision matrices for both the predict unsualness analysis and spatiophylogenetic modelling with INLA.
precision_matrices_fn <- "output/spatiophylogenetic_modelling/processed_data/precision_matrices.RDS"
if(!(file.exists(precision_matrices_fn))){
  source("spatiophylogenetic_modelling/analysis/simulations/make_precisionmatrices.R")}

precision_matrices = readRDS(precision_matrices_fn)
phylo_prec_mat = precision_matrices$phylogenetic_precision
spatial_prec_mat = precision_matrices$spatial_precision

#INLA phylo only
source("spatiophylogenetic_modelling/analysis/INLA_parameters.R")
source("spatiophylogenetic_modelling/install_inla.R")

gb <- tree_tips_df %>% 
  inner_join(gb, by =  "Language_ID")

dual_model = inla(formula =Surprisal ~
                    f(spatial_id,
                      model = "generic2",
                      Cmatrix = spatial_prec_mat,
                      hyper = pcprior) +
                    f(phylo_id,
                      model = "generic2",
                      Cmatrix = phylo_prec_mat,
                      hyper = pcprior),
                  data = gb,
                  control.family = list(hyper = list(prec = list(initial = log(1e+08), fixed = TRUE))))


formula <- Surprisal~
  (1 | gr(Glottocode, cov = spatial_covar_mat)) +
  (1 | gr(Glottocode, cov = phylo_covar_mat))  (1|Family)+
  Endangerement

# Function that obtains a predictive model of surprisal
model_surprisal<-function(df) {
  m<-brm(formula = formula,
         data=df,
         chains = 4,
         iter = 12000,
         warmup = 5000,
         cores = 4,
         control = list(adapt_delta=0.99),
         backend="cmdstanr")
  return(m)}

# From now on we centered our analyses on two estimates: LCA and kernel-20
gb_redux<-gb %>%
  filter(Estimator %in% c("prob_lca","prob_ker_20"))

# Get models
models_surprisal<-plyr::dlply(gb_redux,"Estimator",model_surprisal)

# Check summary
summary(model_surprisal_lca)
conditional_effects(model_surprisal_lca)