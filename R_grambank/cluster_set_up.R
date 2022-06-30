r  = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

options(Ncpus = 10)

packageurl <- "https://cran.r-project.org/src/contrib/Archive/nloptr/nloptr_2.0.1.tar.gz"
install.packages(packageurl, repos=NULL, type="source")