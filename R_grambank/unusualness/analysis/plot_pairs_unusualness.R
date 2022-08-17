# Load pkgs
source("requirements.R")

surprisal_fn <- paste0(OUTPUTDIR_tables, "/surprisal.tsv")
if(!file.exists(surprisal_fn)){
  source("unusualness/analysis/get_unusualness_bayesLCA.R")
}
gb <- read_tsv(file = surprisal_fn, show_col_types = F) 

# Plot the pairwise relations between probabilities
pair_lower <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_smooth(method=loess, fill="deepskyblue", color="deepskyblue", ...)
  p
}

pair_diag<- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_density(fill="deepskyblue", color="deepskyblue", ...)
  p
}



pair_upp <- function(data, mapping, method, ...){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  corr <- cor(x, y, method=method, use='complete.obs')
  
  ggally_text(
    label = as.character(round(corr, 2)), 
    mapping = aes(),
    xP = 0.5, yP = 0.5,
    color = 'deepskyblue',
    size=6,
    ...
  )
}


gb %>%
  dplyr::select(c(Surprisal,Estimator,Language_ID)) %>%
  pivot_wider(id_cols = Language_ID,names_from = Estimator,values_from = Surprisal) %>%
  GGally::ggpairs(columns=c("Bayesian LCA","Kernel 1","Kernel 5","Kernel 10","Kernel 15","Kernel 20","Kernel 25","Kernel 30","Kernel 40"),
                  mapping = aes(alpha = 0.01),
                  lower = list(continuous=pair_lower),
                  diag=list(continuous=pair_diag),
                  upper = list(continuous=wrap(pair_upp,
                                               method="spearman")))+
  theme_minimal()+
  theme(panel.grid= element_blank(),
        axis.text = element_blank())