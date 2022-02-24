source("requirements.R")

# 1. Install BiocManager using: 
pacman::p_load("BiocManager")

# 2. Install INLA dependencies with BiocManager using: 
BiocManager::install(c("graph","Rgraphviz",
                       "rgdal",
                       "rgl",
                       "spdep"))

# 3. Update foreach (although it unclear how vital this step was) using: 
#install.packages("foreach")

# 4. Install INLA using: 
# NOTE: This is a big download
install.packages("INLA", repos=c(getOption("repos"), 
                                 INLA="https://inla.r-inla-download.org/R/stable"), 
                 dep=TRUE)
