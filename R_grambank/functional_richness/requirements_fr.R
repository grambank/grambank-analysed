source("requirements.R")

#installing and loading extra packages for Simon and Damian's code for functional richness.

h_load(verbose = F, pkg = c("mFD" , "ggforce", "concaveman", "cowplot", "fundiversity"))

OUTPUTDIR <- "output/functional_richness"

if(!dir.exists(OUTPUTDIR)){
  dir.create(OUTPUTDIR)
}

source("functional_richness/fun_def_make_group_matrix.R")