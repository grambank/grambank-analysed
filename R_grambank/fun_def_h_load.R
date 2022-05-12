
h_load <- function(pkg, verbose = T){
#  p <- "ggpubr"
#  pkg <- (geiger, ape)
    for(p in pkg){
      if(!(p %in% rownames(installed.packages()))){
        install.packages(p)

        if(verbose == T){
                cat(paste0("Installed", p, ". "))}
        
          }
      if(verbose == T){
        library(p, character.only = T, quietly = F)
              cat(paste0("Loaded ", p, ", version ",packageVersion(p)
,".\n "))
      }else{
        
        library(p, character.only = T, quietly = T, verbose = F)}
    
      }
  
}