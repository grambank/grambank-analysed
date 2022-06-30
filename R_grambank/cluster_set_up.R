r  = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

options(Ncpus = 10)

cluster <- 1

if(!("nloptr" %in% rownames(installed.packages()))){
packageurl <- "https://cran.r-project.org/src/contrib/Archive/nloptr/nloptr_2.0.1.tar.gz"
install.packages(packageurl, type="source")
}

source("fun_def_h_load.R")

h_load(c("pacman", "car","lme4", "pbkrtest"))

source("global_variables.R")
