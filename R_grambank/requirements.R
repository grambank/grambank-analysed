# Please run this script first to make sure you have all the necessary packages 
# installed for running the rest of the scripts in this R project

#defining expected number of languages
n_total <- 2467
n_imputed <- 1509
n_overlap_imputed_and_jaeger_tree <- 1265


#installing packages
if (!suppressPackageStartupMessages(require("pacman"))) { install.packages("pacman") }

pacman::p_load(
	tidyverse,
	readr,
	fields,
	reshape2,
	broom,
#	plyr,
	broom.mixed, 
	naniar, 
	glue,
	forcats,
	magrittr,
	stringr,
	brms,
	purrr,
	rcompanion, 

	MASS,
	matrixStats,
	cluster, 

	# imputation
	missForest,
	
	#plotting graphs
	scales,
	RColorBrewer,
	ggpubr,
	ggplot2,
#	cowplot,
	ggrepel,
	gplots,
	ggridges,
	grid,
	gridExtra,
	scales,
	ggmap,
	nFactors,
	psych, #for scatterplot matrix
	viridis,
	rlang,
	devtools,
	patchwork,
	
	#making maps
	mapdata,
	maptools,
	mapproj,
	maps,
	
	#sensitivity testing spatial
#	rnaturalearth,
#	rnaturalearthdata,
	ggforce, 
	geosphere,
  sf,
	
	# phylogenetic packages
	ape,
	phytools,
	caper, 
  synchrony,

	# testing
	assertthat,
	testthat
)

# quiet down, tidyverse:
options(tidyverse.quiet = TRUE)
options(warn.conflicts = FALSE)
options(stringsAsFactors = FALSE)


GRAMBANK_LANGUAGES <- file.path("../..", "cldf", "languages.csv")
GRAMBANK_VALUES <- file.path("../..", "cldf", "values.csv")
GRAMBANK_PARAMETERS <- file.path("../..", "cldf", "parameters.csv")
GRAMBANK_CODES <- file.path("../..", "cldf", "codes.csv")

# The columns specifier for readr to parse ../cldf/values.csv
VALUES_COLSPEC <- c(
  ID = col_character(),
  Language_ID = col_character(),
  Parameter_ID = col_character(),
  Value = col_character(),
  Code_ID = col_character(),
  Comment = col_character(),
  Source = col_character()
)

LANGUAGES_COLSPEC = c(
  ID = col_character(),
  Name = col_character(),
  Macroarea = col_character(),
  Latitude = col_double(),
  Longitude = col_double(),
  Glottocode = col_character(),
  ISO639P3code = col_logical(),
  Coders = col_character(),
  provenance = col_character(),
  Family_name = col_character(),
  Family_level_ID = col_character(),
  Language_level_ID = col_character(),
  level = col_character(),
  lineage = col_character()
)

PARAMETERS_COLSPEC = c(
    ID = col_character(),
    Name = col_character(),
    Description = col_character(),
    patron = col_character(),
    name_in_french = col_character(),
    Grambank_ID_desc = col_character(),
    bound_morphology = col_character()
)

CODES_COLSPEC = c(
    ID = col_character(),
    Parameter_ID = col_character(),
    Name = col_character(),
    Description = col_character()
)

WIDE_COLSPEC = c(
    .default = col_integer(),
    Language_ID = col_character(),
    na_prop = col_double()
)

#parameters for spatial covariance matrix. values that are used across multiple scripts

#for detect_coderbias.R and spatiophylogenetic_jaegermodel.R
kappa = 2 # smoothness parameter as recommended by Dinnage et al. (2020)
sigma = c(1, 1.15) # Sigma parameter. First value is not used. 

#for waic_test.R and spatial_parameters.R
kappa_vec = c(2, 4, 1, 2, 2)
sigma_vec =  c(1.15, 2, 40, 10, 20)