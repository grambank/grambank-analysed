r  = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

options(Ncpus = 10)

kappa = 2
sigma = c(1, 1.15)
