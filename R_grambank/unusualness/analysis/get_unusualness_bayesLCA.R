# Analysis of Grambank unusualness
# Code by D. E. Blasi

#########################################
## (1) Setup
#########################################

# Set working directory for output
#setup output dirs
OUTPUTDIR <- file.path("output/unusualness/")
if (!dir.exists(OUTPUTDIR)) { dir.create(OUTPUTDIR) }		

OUTPUTDIR_tables <- file.path("output/unusualness/tables")
if (!dir.exists(OUTPUTDIR_tables)) { dir.create(OUTPUTDIR_tables) }		

OUTPUTDIR_plots <- file.path("output/unusualness/plots")
if (!dir.exists(OUTPUTDIR_plots)) { dir.create(OUTPUTDIR_plots) }		

# Load pkgs
source("requirements.R")

# Load imputed binarized GB data, Glottolog families, AUTOTYP areas, and AES endangerment
# Load imputed binarized GB data
if(!file.exists("output/GB_wide/GB_wide_imputed_binarized.tsv")){
  source("make_wide.R")
  source("make_wide_binarized.R")
  source("impute_missing_values.R")
}

gb<-read_tsv("output/GB_wide/GB_wide_imputed_binarized.tsv", col_types = cols())

# Load Glottolog families, AUTOTYP areas, and AES endangerment
# Add autotyp areas to language tibble
AUTOTYP_area_fn <- "output/non_GB_datasets/glottolog_AUTOTYP_areas.tsv"
if(!file.exists(AUTOTYP_area_fn)){
  source("unusualness/processing/assigning_AUTOTYP_areas.R")
}
autotyp_areas <- read_tsv(AUTOTYP_area_fn, col_types = cols()) %>%
  dplyr::select(Language_ID, AUTOTYP_area) 

glottolog_fn <- "output/non_GB_datasets/glottolog-cldf_wide_df.tsv"
if(!file.exists(glottolog_fn)) {
  source("make_glottolog-cldf_table.R")
}
glottolog_df <-read_tsv(glottolog_fn, col_types = cols()) %>%
  mutate(Language_level_ID = ifelse(is.na(Language_level_ID), Language_ID, Language_level_ID)) %>%
  mutate(Family_ID = ifelse(is.na(Family_ID), Language_level_ID, Family_ID)) %>%  # If no language family is assigned to a language then it's an isolate - in which case, use the name of the language as that of the family
  dplyr::select(Language_ID,Macroarea,Family_ID,aes,level,Name) 

gb <- gb %>% 
  left_join(autotyp_areas,  by = "Language_ID") %>% 
  left_join(glottolog_df,  by = "Language_ID")
  
# Check the coverage of the covariates
na_macroarea <- sum(!is.na(gb$Macroarea))/nrow(gb)
na_autotyp_area <- sum(!is.na(gb$AUTOTYP_area))/nrow(gb)
na_family <- sum(!is.na(gb$Family_ID))/nrow(gb)
na_aes <-sum(!is.na(gb$aes))/nrow(gb) 

cat(paste0(
  "There is \n",
  100 - (na_macroarea*100), "% missing data for macroarea\n",
  100 - (na_autotyp_area*100), "% missing data for AUTOTYP-area\n",
  100 - (na_family*100), "% missing data for Family memberhsip\n",
  round(100 - (na_aes*100), 2), "% missing data for aes-status.\n"
))


# Set a dummy variable that is all GB features
gb_features<-colnames(gb)[startsWith(colnames(gb),"GB")]

#########################################
## (2) Density estimation - kernel smoothing
#########################################

# Local kernel density estimation
# Idea: estimate the probability density of a given grammar based on the distribution of neighbors - 
# they contribute to the probability mass with an exponential weight of their distance to the grammar

# Investigate the distribution of Gower distances in the data
gb_dists<-cluster::daisy((gb[,gb_features]), metric = "gower", warnBin = F)

# Plot
gb_dists %>%
  hist()

# Define an exponential kernel function based on parameter k and distances ds
kernel_k<-function(ds,k){
  return(sum(exp(-k*ds))/length(ds))}

# Attach our probability estimates to a given df and a matrix of distances
add_kernel_probabilities<-function(df,distances,k){
  df[[paste0("prob_ker_",k)]]<-apply(distances,1,function(x) kernel_k(x,k))
  return(df)}

# Enrich GB data
gb<-gb %>%
  add_kernel_probabilities(gb_dists,1) %>%
  add_kernel_probabilities(gb_dists,5) %>%
  add_kernel_probabilities(gb_dists,10) %>%
  add_kernel_probabilities(gb_dists,20) %>%
  add_kernel_probabilities(gb_dists,30) %>%
  add_kernel_probabilities(gb_dists,40)
  

#########################################
## (3) Density estimation - Bayesian LCA
#########################################

# Bayesian LCA
# Idea: we organize GB features into correlated bundles and then we induce latent probabilistic classes on them

# First we sort all GB features into independent bundles such that: 
# I. they are highly correlated within them and 
# II. reasonably independent from features in other bundles

# Get gap statistic for each number of clusters
#test <- factoextra::fviz_nbclust(t(gb[,gb_features]), FUN = hcut, method = "gap_stat")

#test2 <- cluster::clusGap(t(gb[,gb_features]), FUN = hcut, method = "gap_stat", K.max = 30)

# 9 clusters selected as optimal according to this criterion
n_clusters<-9

# Obtain a hierarchical clustering of the features with that number of clusters
hier_gb<-hcut(t(gb[,gb_features]), k = n_clusters, stand = TRUE)
fviz_dend(hier_gb, rect = TRUE)

# Assign cluster number to each GB feature
hier_classes<-hier_gb$cluster

# Here we take the clusters from the previous section, and within each we induce an optimal 
# number of latent classes using a Bayesian implementation of LCA

# Set a range of latent clusters to explore for each cluster
k_range<-c(1:5)

# Function that yields the BIC of latent classes assignments for feature cluster i
k_choose<-function(i) {
  
    LCA_clusters<-lapply(k_range,
                     function(k) blca.em(gb[,names(hier_classes[hier_classes==i])],
                                         k, 
                                         restarts=250))

    # Get BICs
    LCA_BIC<-data.frame(BIC=sapply(LCA_clusters, function(x) x$BIC),
                    k=k_range,
                    cluster_n=i)
    
    return(LCA_clusters[[which(LCA_BIC$BIC==min(LCA_BIC$BIC))]])}

# Apply this to all clusters
LCA_clusters<-lapply(c(1:n_clusters),function(x) k_choose(x))

# Finally, and based on the latent classes obtained before, attach rarities to each observation given a cluster

# Function that attaches probabilities and surprisals according to this model *within each cluster*
estimate_unusualness_LCA<-function(lca,df) {
  
  p_item<-lca$itemprob # Get the latent probabilities
  
  class<-apply(Zscore(df,lca),
               1,
               function(x) which(x==max(x)))
  
  prob_assignment<-apply(Zscore(df,lca),
               1,
               function(x) max(x))
  
  prob<-sapply(1:nrow(df),
                  function(y) estimate_prob_LCA(df[y,],p_item[class[y],]))
  
  return(data.frame(prob=prob,
                    assignment=prob_assignment,
                    lg=1:length(prob)))
  }

# Auxuliary function that determines probability
estimate_prob_LCA<-function(obs,p){
  return(prod(obs*p+(1-obs)*(1-p)))}

# Apply
cluster_list<-c(1:n_clusters)
names(cluster_list)<-cluster_list
  
unusualness_df<-plyr::ldply(cluster_list,
                       function(x) estimate_unusualness_LCA(LCA_clusters[[x]],gb[,names(hier_classes[hier_classes==x])]),
                       .id="Cluster") %>%
  plyr::ddply("lg",function(x) data.frame(prob_lca=prod(x$prob)))


unusualness_df %>% write_tsv(file = paste0(OUTPUTDIR_tables, "unusualness.tsv"))

# Enrich GB data
gb$prob_lca<-unusualness_df$prob_lca

#########################################
## (4) Compare and visualize probabilities and surprisals
#########################################

# Standardize the resulting probabilities and compute surprisals
gb <- gb %>%
  pivot_longer(cols=c(prob_ker_1,prob_ker_5,prob_ker_10,prob_ker_20,prob_ker_30,prob_ker_40,prob_lca),names_to = "Estimator",values_to="Probability") %>%
  mutate(Surprisal=-log(Probability))

# Re-label for clarity
gb$Estimator<-fct_recode(gb$Estimator,
                         "Bayesian LCA"="prob_lca",
                         "Kernel 1"="prob_ker_1",
                         "Kernel 5"="prob_ker_5",
                         "Kernel 10"="prob_ker_10",
                         "Kernel 20"="prob_ker_20",
                         "Kernel 30"="prob_ker_30",
                         "Kernel 40"="prob_ker_40")

gb %>% write_tsv(file = paste0(OUTPUTDIR_tables, "/surprisal.tsv"))
