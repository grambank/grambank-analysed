
#function to check if a pkg is installed or not, if not it installs it and either way it's loaded.
#inspired by pacman::p_load()

h_load <- function(pkg, verbose = FALSE, version = NULL, repos = "http://cran.us.r-project.org"){
    for(p in pkg){
        #if no version is specified, check if it's installed and if not then go and install it as normal
        if (is.null(version) & (!(p %in% rownames(installed.packages())))) {
            install.packages(p, dependencies = T, repos = repos)
            if (verbose == TRUE) { cat(paste0("Installed ", p, ".\n")) }
        }
        if (!is.null(version) & (!(p %in% rownames(installed.packages())))) {

        #install devtools if it isn't already installed
        if (!"devtools" %in% rownames(installed.packages())){
            install.packages(devtools, dependencies = TRUE, repos = repos)
        }

        library(devtools, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE)
        devtools::install_version(p, version = version, dependencies = TRUE, repos = repos)

        if(verbose == TRUE){
            cat(paste0("Installed ", p, ", version ", version, ".\n"))}
        }

        if (verbose == TRUE) {
            library(p, character.only = TRUE, quietly = FALSE)
            cat(paste0("Loaded ", p, ", version ",packageVersion(p),".\n"))
        } else {
            suppressMessages(library(p, character.only = TRUE, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE))
        }
    }
}
