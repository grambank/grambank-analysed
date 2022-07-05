# Please run this script first to make sure you have all the necessary packages
# installed for running the rest of the scripts in this R project


#installing packages
#if (!suppressPackageStartupMessages(require("pacman"))) { install.packages("pacman") }
source("fun_def_h_load.R")

h_load(verbose = F, pkg = c(
#pacman::p_load(
  "spam",
	"dplyr",
  "ggplot2",
  "tidyr",
  "forcats",
  "magrittr",
  "purrr",
  "stringr",
  "tibble",
  	"readr",
	"fields",
	"reshape2",
	"broom",
#	"plyr",
	"broom.mixed",
	"naniar",
	"glue",
	"forcats",
	"magrittr",
	"stringr",
	"brms",
	"purrr",
	"rcompanion",
  "qs",

	"MASS",
	"matrixStats",
	"cluster",
  "MCMCglmm",
"foreach",
"geiger",


	# "imputation
"randomForest",
	"missForest",

	#plotting graphs
  "ellipse",
	"scales",
	"RColorBrewer",
  "car",
	"ggpubr",
	"ggplot2",
#	cowplot",
	"ggrepel",
	"gplots",
	"ggridges",
	"grid",
	"gridExtra",
	"scales",
#	"ggmap",
	"nFactors",
	"psych", #for scatterplot matrix
	"viridis",
	"rlang",
	"devtools",
	"patchwork",

	#making maps
  "maps",
	"mapdata",
#	"maptools",
	"mapproj",

#unsualness
"e1071",
"mvtnorm",
"BayesLCA",
"factoextra",
"GGally",


	#sensitivity testing spatial
#	rnaturalearth",
#	rnaturalearthdata",
"ggforce",
	"geosphere",
"sf",

	# phylogenetic packages
"ape",
"phytools",
	"caper",
"synchrony",
  "geiger",
"synchrony",

	# testing
	"assertthat",
"testthat"
))

# quiet down, tidyverse:
options(tidyverse.quiet = TRUE)
options(warn.conflicts = FALSE)
options(stringsAsFactors = FALSE)

source("global_variables.R")
