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

source("global_variables.R")
source("spatiophylogenetic_modelling/analysis/INLA_parameters.R")
source('spatiophylogenetic_modelling/analysis/functions/strip_inla.R')

#read in data
gb <- read.delim(file = surprisal_fn, sep = "\t") %>% 
  dplyr::select(Language_ID, aes, Surprisal, Estimator) %>% 
  filter(Estimator == "Kernel 30") %>% 
  inner_join(lgs_in_analysis, by = "Language_ID") %>% #subset to the lgs where we have phylo prec matrices
  mutate(Endangerement=ifelse(aes %in% c("threatened","moribund","nearly_extinct"),"endangered",aes)) # Recode endangerment

### NEXT PART REQUIRES MATRICES ETC

#########################################
## (5) Model unusualness in terms of genealogical, areal covariates, and endangerement status
#########################################

#In the spatiophylogenetic modelling of the features, we use the dataset cropped for missing data but without imputation. For the unsualness analsyis, we use the imputed data. They are different in the feature values, but it is the same subset of langauges in both. Therefore, we can use the same precision matrices for both the predict unsualness analysis and spatiophylogenetic modelling with INLA.
precision_matrices_fn <- "output/spatiophylogenetic_modelling/processed_data/precision_matrices.RDS"
if(!(file.exists(precision_matrices_fn))){
  source("spatiophylogenetic_modelling/analysis/simulations/make_precisionmatrices.R")}

precision_matrices = readRDS(precision_matrices_fn)
phylo_prec_mat = precision_matrices$phylogenetic_precision
spatial_prec_mat = precision_matrices$spatial_precision

#reading in AUTOTYP-area
if (!file.exists("output/non_GB_datasets/glottolog_AUTOTYP_areas.tsv")) { s
  source("unusualness/processing/assigning_AUTOTYP_areas.R") }		
autotyp_area <- read.delim("output/non_GB_datasets/glottolog_AUTOTYP_areas.tsv", sep = "\t") %>%
  dplyr::select(Language_ID, AUTOTYP_area_id_iid_model = AUTOTYP_area)

data <- gb %>% 
  left_join(autotyp_area, by = "Language_ID")

x <- assert_that(all(data$Language_ID == lgs_in_analysis$Language_ID), msg = "Data doesn't match!")

## Since we are using a sparse phylogenetic matrix, we need to math taxa to the correct
## rows in the matrix
data$phylo_id = match(data$Language_ID, rownames(phylo_prec_mat))
data$spatial_id = match(data$Language_ID, rownames(spatial_prec_mat))
data$obs_id = 1:nrow(data)

#INLA phylo only
#dual model
source("spatiophylogenetic_modelling/install_inla.R")

dual_model = INLA::inla(Surprisal ~
                    f(spatial_id,
                      model = "generic0",
                      Cmatrix = spatial_prec_mat,
                      hyper = pcprior) +
                    f(phylo_id,
                      model = "generic2",
                      Cmatrix = phylo_prec_mat,
                      hyper = pcprior),
                    control.compute = list(waic = TRUE, cpo = TRUE),
                  data = data#,
                  #control.family = list(hyper = list(prec = list(initial = log(1e+08), fixed = TRUE)))
                  )



dual_model_stripped <- strip_inla(dual_model)
