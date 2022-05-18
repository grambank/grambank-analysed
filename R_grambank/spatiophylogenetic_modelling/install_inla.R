R_version <- 4.1
testing <- "no"
experimental <- "yes"

source("fun_def_h_load.R")  

h_load("BiocManager")

  # Install INLA dependencies with BiocManager using: 
  BiocManager::install(c("graph","Rgraphviz",
                         "rgdal",
                         "rgl",
                         "spdep"))
  
  # Update foreach (although it unclear how vital this step was) using: 
 
 h_load("foreach")
  
  # 4. Install INLA using: 
  # NOTE: This is a big download

    if(testing != "yes"){
    
      if (!("INLA" %in% rownames(installed.packages()))) { 
      cat("INLA wasn't installed, installing now.\n") 
    
      install.packages("INLA", repos=c(getOption("repos"), 
                                   INLA="https://inla.r-inla-download.org/R/stable"), 
                   dep=TRUE)
  }}else{
  
  install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
  }
  

suppressPackageStartupMessages(
  library(INLA, quietly = T, warn.conflicts = F, verbose = F)
  )

if(experimental== "yes"){
inla.setOption(inla.mode="experimental")
}

cat(paste0("Loaded INLA version ", packageVersion("INLA"), ".\n"))
