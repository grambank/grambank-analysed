
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
  source("unusualness/analysis/get_unsualness_bayesLCA.R")
}
gb <- read_tsv(file = surprisal_fn)

### NEXT PART REQUIRES MATRICES ETC

#########################################
## (5) Model unsualness in terms of genealogical, areal covariates, and endangerement status
#########################################

# Recode endangerment
gb<-gb %>%
  mutate(Endangerement=ifelse(aes %in% c("threatened","moribund","nearly_extinct"),"endangered",aes))


spatial_covar_mat_fn <- "output/spatiophylogenetic_modelling/spatial_covar_mat.tsv"
if(!file.exists(spatial_covar_mat_fn)){
  source("spatiophylogenetic_modelling/analysis/make_vcvs.R")
}

spatial_covar_mat <- read_tsv(spatial_covar_mat_fn, show_col_types = F) %>% 
  column_to_rownames("Language_ID") %>% 
  as.matrix()

phylo_covar_mat_fn <- "output/spatiophylogenetic_modelling/phylo_covar_mat.tsv"
if(!file.exists(phylo_covar_mat_fn)){
  source("spatiophylogenetic_modelling/analysis/make_vcvs.R")
}

phylo_covar_mat <- read_tsv(phylo_covar_mat_fn, show_col_types = F) %>% 
  column_to_rownames("Language_ID") %>% 
  as.matrix()


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