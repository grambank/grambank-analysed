source("requirements.R")

#installing and loading extra packages for Simon and Damian's code for functional richness.

h_load(verbose = F, pkg = c("fundiversity", "mFD"  ))

OUTPUTDIR <- "output/functional_richness"

if(!dir.exists(OUTPUTDIR)){
  dir.create(OUTPUTDIR)
}