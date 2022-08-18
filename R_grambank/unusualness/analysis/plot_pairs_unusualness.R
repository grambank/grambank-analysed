# Load pkgs
source("requirements.R")

# Set working directory for output
#setup outpur dirs
OUTPUTDIR <- file.path("output/unusualness/")
if (!dir.exists(OUTPUTDIR)) { dir.create(OUTPUTDIR) }		

OUTPUTDIR_tables <- file.path("output/unusualness/tables")
if (!dir.exists(OUTPUTDIR_tables)) { dir.create(OUTPUTDIR_tables) }		

OUTPUTDIR_plots <- file.path("output/unusualness/plots")
if (!dir.exists(OUTPUTDIR_plots)) { dir.create(OUTPUTDIR_plots) }		

surprisal_fn <- paste0(OUTPUTDIR_tables, "/surprisal.tsv")
if(!file.exists(surprisal_fn)){
  source("unusualness/analysis/get_unusualness_bayesLCA.R")
}
gb <- read_tsv(file = surprisal_fn, show_col_types = F) 

# Plot the pairwise relations between probabilities
pair_lower <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_smooth(method=loess, fill="purple4", color="purple4", ...)
  p
}

pair_diag<- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_density(fill="purple4", color="purple4", ...)
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
    color = 'purple4',
    size=8,
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

ggsave(filename = paste0(OUTPUTDIR_plots, "/unsualenss_SLOM_compare_kernels.png"))
