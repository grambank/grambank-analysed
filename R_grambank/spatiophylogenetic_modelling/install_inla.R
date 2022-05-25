R_version <- 4.2
testing <- "yes"
experimental <- "yes"

source("fun_def_h_load.R")  

h_load("BiocManager")

if(!("spdep" %in% rownames(installed.packages()))){
  # Install INLA dependencies with BiocManager using: 
  BiocManager::install(c("graph","Rgraphviz",
                         "rgdal",
                         "rgl",
                         "spdep")) 
  }
  
  # Update foreach (although it unclear how vital this step was) using: 
 
 h_load("foreach")
  
  # 4. Install INLA using: 
  # NOTE: This is a big download

#If inla isn't installed, we don't want the testing version and R is of a version 4.2 or higher, do this
if (!("INLA" %in% rownames(installed.packages())) & testing != "yes" & R_version >=4.2) { 
      cat("INLA wasn't installed, installing stable version now.\n") 
    
      install.packages("INLA", repos=c( 
                                   INLA="https://inla.r-inla-download.org/R/stable"), 
                   dep=TRUE)
  }
  
 #If inla isn't installed, we do want the testing version and R is of a version 4.2 or higher, do this
 if (!("INLA" %in% rownames(installed.packages())) & testing == "yes" & R_version >= 4.2){
         install.packages("INLA",repos=c(INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
 }
 
 #If inla isn't installed, we do want the testing version and R isn't of a version 4.2 or higher, do this 
 if (!("INLA" %in% rownames(installed.packages())) & testing == "yes" & R_version < 4.2){
   h_load("remotes")
   remotes::install_version("INLA", version="22.05.03",repos=c(INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
 }
 
#whatever you do, now load INLA
suppressPackageStartupMessages(
  library(INLA, quietly = T, warn.conflicts = F, verbose = F)
  )



if(experimental== "yes"){
  INLA::inla.setOption(inla.mode="experimental")
}

cat(paste0("Loaded INLA version ", packageVersion("INLA"), ".\n"))


