verbose <- F

if(verbose == T){
  cat("For the analysis in this paper, we are using the testing-version of INLA, beta if you will, version = INLA_22.06.03. If you have another version of INLA running, please set it to this for maximum reproducibility.")
  }


if(!(INLA %in% rownames(installed.packages()))){

source("fun_def_h_load.R")  

h_load(pkg = c("Deriv", "Ecdat", "HKprocess", "gsl", "mlogit", "orthopolynom", "rgdal", "rgeos", "sn", "splancs"))

install.packages("INLA",repos=c(INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)

}

library("INLA")
INLA::inla.setOption(inla.mode="experimental")


