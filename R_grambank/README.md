# Basic

This directory (R_grambank) contains the *R* scripts needed to generate the analysis, tables and figures for the Grambank paper: "Grambank reveals the importance of genealogical constraints on linguistic diversity and highlights the impact of language loss".

These scripts rely on data from git submodules. Read more about that [here](https://github.com/grambank/grambank-analysed#git-submodules).

This codebase is organised on the principle that we do not commit output of scripts to git that can be recreated by scripts here. For example, in order to run the scripts relating to PCA scatterplots you will first need to run the scripts to generate the PCA tables, which in turn rely on the make_wide scripts which bring in the Grambank data itself. The makefile rules "make data" and "make pca" run all the necessary scripts in the right order for this.

## Running code in the right order
The R-scripts can be run one by one, or called on from the Makefile in the directory R_grambank.

### setting up the data
In order to run all the scripts, we first need to just set-up the data-tables etc. This is the order to run those scripts in:

1.	make_glottolog-cldf_table.R                                                                 
2.	make_wide.R                                                                                 
3.	make_wide_binarized.R                                                                       
4.	impute_missing_values.R                                                                     
5.	make_theo_score_fusion.R                                                                    
6.	make_theo_scores.R                                                                          
7.	spatiophylogenetic_modelling/processing/pruning_EDGE_tree.R                                 
8.	unusualness/processing/assigning_AUTOTYP_areas.R                                            
### PCA
Once the data  is set up, you will ned to run PCA/PCA.R first in order to run any of the PCA scripts.                                                                   
### Functional Richness
The first script to run for the functional richness analysis is 

1. functional_richness/make_gower_mfd.R                                                        

### Unusualness 
The unsualness scripts to run first are:

1.   unusualness/analysis/get_unusualness_bayesLCA.R                                             
2.   unusualness/analysis/predict_surprisal.R                      

### Spatiophylogenetic modelling (SP)
The SP-scripts are the most complicated in this codebase, and involve running simulations and mutliple different values for the pcprior and spatial decay. The core script that the release paper relies most on is:

1.   spatiophylogenetic_modelling/analysis/INLA_multi_models.R     

## Running code in the right order: Makefile
`Makefile`: is a [GNU Make](https://www.gnu.org/software/make/) Makefile for compiling the analysis for the paper. We have organised the analysis in this repos into three categories: quick, medium and long. This represents the expected run time. The long analysis takes over 12h to compile, the quick just a few minutes and the medium under an hour on a typical personal computer.

The Makefile is organised into rules by type of analysis, for example for PCA, spatiophylogenetic modelling etc. These rules have been further summarised into:

1. quick_ones (wrangle data, PCA, coverage plots and Manhattan distances)
2. medium_length (compute unusualness and cultural fixation scores)
3. if_medium_are_finished (plotting  and cultural fixation scores)
4. long_ones (spatiophylogenetic modelling of real and simulated data, functional richness, predict unusualness and predict spatiophylogenetic effects, ancestral state reconstruction by spatiophylogenetic modelling)
5. if_long_ones_are_finished (plotting and summarising of various spatiophylogenetic output, functional richness etc in "long_ones")

## The data

Throughout these scripts "with question" refers to versions of the dataset where the datapoints that are coded as "?" are maintained as "?" and not coerced into missing values. These "?" values are instances where a coder has accessed descriptions and tried to make a coding judgment, but where the information was insufficient or non-existent.

Fully missing data is when no attempt has been made yet (often the case with inherited data).

The term "strict" refers to cases where "?" have indeed been coerced into missing data points, which is more appropriate for much of the analysis. These "with question" versions are mainly useful to calculate current coverage whereas "strict" should be used for actual analysis.

In short:

* "with question" = ? remain ?
* "strict" = ? are turned into missing values.

In these scripts, "na_prop" refers to the proportion of missing values. A high value means that many values are missing. The wide dataset both come with this information over each language in the second column. The imputed data does not come with a na_prop column, since there is not missing data because those gaps have been imputed.

### Git submodules
This Git repository contains git submodules. That means that this repository is linked to other git repositories in a principled way. In this instance this repository has git submodules for the following repostiroeies: grambank-cldf, AUTOTYP-data, glottolog-cldf and WALS.

If you want to run scripts in this repository on your machine, it is necessary not only to clone this repository but also after cloning to run:

`git submodule update --init`

This command will initialise and update the git submodules appropriately. Note that this includes data from grambank-cldf, so no script will run without initalising the git submodules.

You can read more about git submodules [here](https://git-scm.com/book/en/v2/Git-Tools-Submodules#_cloning_submodules).


## Requirements

To run these files you will need to be able to install packages, we use the [pacman](https://github.com/trinker/pacman) to handle this, and the script `requirements.R` will install these for you. Please make sure that *R* can install packages before running this. In addition, we also install a special package - INLA - which is needed for the spatiophylogenetic modelling. This package is installed with a special script: spatiophylogenetic_modelling/install_inla.R.

The package versions needed to run these scripts are:

* assertthat 0.2.1
* backports 1.1.4
* bitops 1.0-6
* broom 0.5.2
* callr 3.2.0
* cellranger 1.1.0
* cli 1.1.0
* codetools 0.2-16
* colorspace 1.4-1
* crayon 1.3.4
* desc 1.2.0
* devtools 2.0.2
* digest 0.6.18
* dplyr 0.8.0.1
* forcats 0.4.0
* foreach 1.4.4
* foreign 0.8-71
* fs 1.2.7
* generics 0.0.2
* gg3D 0.0.0.9000
* ggmap 3.0.0
* ggplot2 3.1.1
* ggpubr 0.2
* glue 1.3.1
* gridExtra 2.3
* gtable 0.3.0
* haven 2.1.0
* hms 0.4.2
* httr 1.4.0
* iterators 1.0.10
* itertools 0.1-3
* jpeg 0.1-8
* jsonlite 1.6
* lattice 0.20-38
* lazyeval 0.2.2
* lubridate 1.7.4
* magrittr 1.5
* mapdata 2.3.0
* maps 3.3.0
* maptools 0.9-5
* MASS 7.3-51.4
* memoise 1.1.0
* misc3d 0.8-4
* missForest 1.4
* modelr 0.1.4
* munsell 0.5.0
* nlme 3.1-139
* pacman 0.5.1
* pillar 1.3.1
* pkgbuild 1.0.3
* pkgconfig 2.0.2
* pkgload 1.0.2
* plot3D 1.1.1
* plyr 1.8.4
* png 0.1-7
* prettyunits 1.0.2
* processx 3.3.0
* ps 1.3.0
* purrr 0.3.2
* R6 2.4.0
* randomForest 4.6-14
* RColorBrewer 1.1-2
* Rcpp 1.0.1
* readr 1.3.1
* readxl 1.3.1
* remotes 2.0.4
* reshape2 1.4.3
* RgoogleMaps 1.4.3
* rjson 0.2.20
* rlang 0.3.4
* rprojroot 1.3-2
* rstudioapi 0.10
* rvest 0.3.3
* scales 1.0.0
* sessioninfo 1.1.1
* sp 1.3-1
* stringi 1.4.3
* stringr 1.4.0
* tibble 2.1.1
* tidyr 0.8.3
* tidyselect 0.2.5
* tidyverse 1.2.1
* usethis 1.5.0
* viridis 0.5.1
* viridisLite 0.3.0
* withr 2.1.2
* xml2 1.2.0
