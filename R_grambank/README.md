# R scripts for plotting Grambank Data.

This directory contains the *R* scripts needed to generate the figures for the Grambank paper.

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

## The data

Throughout these scripts "with question" refers to versions of the dataset where the datapoints that are coded as "?" are maintained as "?" and not coerced into missing values. These "?" values are instances where a coder has accessed descriptions and tried to make a coding judgment, but where the information was insufficient or non-existent.

Fully missing data is when no attempt has been made yet (often the case with inherited data).

The term "strict" refers to cases where "?" have indeed been coerced into missing data points, which is more appropriate for analysis. These "with question" versions are mainly useful to calculate current coverage whereas "strict" should be used for actual analysis.

In short:

* "with question" = ? remain ?
* "strict" = ? are turned into missing values.

In these scripts, "na_prop" refers to the proportion of missing values. A high value means that many values are missing. The wide dataset both come with this information over each language in the second column. The imputed data does not come with a na_prop column, since there is not missing data because those gaps have been imputed.

## Files

`Makefile`: is a [GNU Make](https://www.gnu.org/software/make/) Makefile for compiling and building all the plots. To build everything, run this from your terminal/console in the directory R_grambank:

1. `make requirements.log`: installs all the required R packages.
2. `make data`: renders wide data tables based on the CLDF-data. This step also generates the cropped and imputed dataframe as well as the glottolog-cldf table and pruning of the global Jäger-tree. This step also generates the data coverage comparison plot for WALS vs GB, which necessitates fetching WALS data from the git submodule.
3. `make pca`: runs PCA and generate all PCA plots. Also compares PCA to theoretical scores
4. `make get_unusualness`: calculates the unusualness scores per languages and generates plots.
5. `make phylo_signal_per_feature`: calculates the D-value per feature (in the binarised, cropped and imputed dataset)
6. `make endangerment_analysis`: calculates the dissimilarity scores per endangerment level

If you want to run all steps 1:6, you can also use the Makefile rule `make almost_all_fast`. This will run all the steps which are reasonable to run on a personal computer and which do not take up a lot of time. The next steps are more computationally expensive, which is why you may want to elect to run them when there is more time or more computational resources (cluster).

7. `make INLA`: runs the spatiophylogenetic modelling using the INLA-approach. This is feasible to run on a personal computer, but may take a few hours.
8. `make predict_unsualness`: runs a BRMS-analysis of the unusualness score calculated in step 3. This is preferable to run on a cluster instead of on a personal computer.

If you wish to run all analysis in one sweep, including step 6:7 and tests: `make all` will accomplish this.

9. `make clean`: deletes all new files and directories from previous steps, wiping the slate clean for running scripts anew.

### Git submodules
All scrits require git submodules to be initialised. Git submodules are a principled way to linked to other git repositories. If you want to run scripts in this repository on your machine, it is necessary not only to clone this repository but also after cloning to run:

`git submodules update --init`

This command will initialise and update the git submodules appropraitely. This includes the data from grambank-cldf. You can read more about git submodules [here](https://git-scm.com/book/en/v2/Git-Tools-Submodules#_cloning_submodules).
