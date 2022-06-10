verbose <- F
cluster = 0

if(verbose == T){
  cat("For the analysis in this paper, we are using the testing-version of INLA, beta if you will, version = INLA_22.06.03. If you have another version of INLA running, please set it to this for maximum reproducibility.")
  }


if(!("INLA" %in% rownames(installed.packages())) & cluster == 0){

source("fun_def_h_load.R")  

h_load(pkg = c("Deriv", "Ecdat", "HKprocess", "gsl", "mlogit", "orthopolynom", "rgdal", "rgeos", "sn", "splancs"))

install.packages("INLA",repos=c(INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)

}

if(!("INLA" %in% rownames(installed.packages())) & cluster == 1){
  
  if(!(dir.exists("../rlibs"))) {
mkdir("../rlibs")  
  }
  
install.packages("INLA",repos=c(INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE, lib = "../rlibs")
INLA::inla.binary.install(os = "Ubuntu-18.0", path = "../rlibs/INLA/bin/linux/")
}

library("INLA")
INLA::inla.setOption(inla.mode="experimental")



