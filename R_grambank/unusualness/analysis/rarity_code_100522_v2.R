# Analysis of Grambank unusualness
# Code by D E Blasi

#########################################
## (1) Setup
#########################################

# Set working directory for output
#setup outpur dirs
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
k_range<-c(1:6)

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
estimate_rarity_LCA<-function(lca,df) {
  
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
  
rarity_df<-plyr::ldply(cluster_list,
                       function(x) estimate_rarity_LCA(LCA_clusters[[x]],gb[,names(hier_classes[hier_classes==x])]),
                       .id="Cluster") %>%
  plyr::ddply("lg",function(x) data.frame(prob_lca=prod(x$prob)))

# Enrich GB data
gb$prob_lca<-rarity_df$prob_lca
rm("rarity_df") # Remove unnecessary df

#########################################
## (4) Compare and visualize probabilities and surprisals
#########################################


# Standardize the resulting probabilities and compute surprisals
gb<-gb %>%
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

# Plot the pairwise relations between probabilities
gb %>%
  select(c(Surprisal,Estimator,Language_ID)) %>%
  pivot_wider(id_cols = Language_ID,names_from = Estimator,values_from = Surprisal) %>%
  GGally::ggpairs(columns=c("Bayesian LCA","Kernel 1","Kernel 5","Kernel 10","Kernel 20","Kernel 30","Kernel 40"),
                  mapping = aes(alpha = 0.1))+
  theme_minimal()

ggsave("comparison_surprisals.png",height=8,width = 10)

save.image("partial_workspace_rarity_160522.rdata")

### NEXT PART REQUIRES MATRICES ETC

# Zooming into the LCA and the kernel-20 approaches, and highlighting Macroarea
gb %>%
  select(c(Surprisal,Estimator,Name,Family,Macroarea)) %>%
  filter(Estimator %in% c("prob_lca","prob_ker_20")) %>%
  pivot_wider(id_cols = c(Name,Family,Macroarea),names_from = Estimator,values_from = Surprisal) %>%
  mutate(IE=ifelse(Family=="indo1319","IE","nIE")) %>%
  ggplot(aes(x=prob_lca,y=prob_ker_20,label=Name,color=Macroarea))+
  geom_text(alpha=0.7)+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x="Surprisal based on LCA approximation",y="Unnormalized surprisal based on 20-kernel")

ggsave("comparison_surprisals_pairwise_macroareas.png",width=25,height=18)

gb[gb$aes %in% c("threatened","not_endangered","nearly_extinct","moribund"),] %>%
  select(c(Surprisal,Estimator,Name,aes)) %>%
  filter(Estimator %in% c("prob_lca","prob_ker_20")) %>%
  pivot_wider(id_cols = c(Name,aes),names_from = Estimator,values_from = Surprisal) %>%
  mutate(Endangerment=ifelse(aes=="not_endangered","Safe","Not safe")) %>%
  ggplot(aes(x=prob_lca,y=prob_ker_20,label=Name,color=Endangerment))+
  geom_text(alpha=0.7)+
  theme_bw()+
  theme(legend.position = c(0.8,0.2))+
  labs(x="Surprisal based on LCA approximation",y="Unnormalized surprisal based on 20-kernel")

ggsave("comparison_surprisals_pairwise_endangerement.png",width=25,height=18)

# Plot individual surprisals
gb %>%
  filter(Estimator=="prob_lca") %>%
  ggplot(aes(x=Surprisal,y=..density..))+
  geom_density(color="dodgerblue3")+
  geom_histogram(bins = 50,alpha=0.7,fill="dodgerblue1")+
  labs(x="Surprisal of Grambank language",y="Density")+
  theme_minimal()+
  annotate("text",label="English",x=gb$Surprisal[gb$Name=="English"&gb$Estimator=="prob_lca"],y=0.04)+
  annotate("text",label="French",x=gb$Surprisal[gb$Name=="French"&gb$Estimator=="prob_lca"],y=0.041)+
  annotate("text",label="Mandarin",x=gb$Surprisal[gb$Name=="Mandarin Chinese"&gb$Estimator=="prob_lca"],y=0.04)+
  annotate("text",label="Turkish",x=gb$Surprisal[gb$Name=="Turkish"&gb$Estimator=="prob_lca"],y=0.041)+
  annotate("text",label="Japanese",x=gb$Surprisal[gb$Name=="Japanese"&gb$Estimator=="prob_lca"],y=0.041)+
  annotate("text",label="Igbo",x=gb$Surprisal[gb$Name=="Igbo"&gb$Estimator=="prob_lca"],y=0.042)+
  annotate("text",label="Mapudungun",x=gb$Surprisal[gb$Name=="Mapudungun"&gb$Estimator=="prob_lca"],y=0.042)+
  annotate("text",label="Yapese",x=gb$Surprisal[gb$Name=="Yapese"&gb$Estimator=="prob_lca"],y=0.041)+
  annotate("text",label="Yanomamö",x=gb$Surprisal[gb$Name=="Yanomamö"&gb$Estimator=="prob_lca"],y=0.039)+
  annotate("text",label="Tamil",x=gb$Surprisal[gb$Name=="Tamil"&gb$Estimator=="prob_lca"],y=0.042)+
  annotate("text",label="Standard Arabic",x=gb$Surprisal[gb$Name=="Standard Arabic"&gb$Estimator=="prob_lca"],y=0.039)+
  annotate("text",label="Russian",x=gb$Surprisal[gb$Name=="Russian"&gb$Estimator=="prob_lca"],y=0.039)+
  annotate("text",label="Hebrew",x=gb$Surprisal[gb$Name=="Hebrew"&gb$Estimator=="prob_lca"],y=0.039)+
  annotate("text",label="Standard Indonesian",x=gb$Surprisal[gb$Name=="Standard Indonesian"&gb$Estimator=="prob_lca"],y=0.039)+
  annotate("text",label="Hungarian",x=gb$Surprisal[gb$Name=="Hungarian"&gb$Estimator=="prob_lca"],y=0.043)
  
gb %>%
  filter(Estimator=="prob_ker_20") %>%
  ggplot(aes(x=Surprisal,y=..density..))+
  geom_density(color="dodgerblue3")+
  geom_histogram(bins = 50,alpha=0.7,fill="dodgerblue1")+
  labs(x="Surprisal of Grambank language",y="Density")+
  theme_minimal()+
  annotate("text",label="English",x=gb$Surprisal[gb$Name=="English"&gb$Estimator=="prob_ker_20"],y=0.74)+
  annotate("text",label="French",x=gb$Surprisal[gb$Name=="French"&gb$Estimator=="prob_ker_20"],y=0.72)+
  annotate("text",label="Mandarin",x=gb$Surprisal[gb$Name=="Mandarin Chinese"&gb$Estimator=="prob_ker_20"],y=0.7)+
  annotate("text",label="Turkish",x=gb$Surprisal[gb$Name=="Turkish"&gb$Estimator=="prob_ker_20"],y=0.78)+
  annotate("text",label="Japanese",x=gb$Surprisal[gb$Name=="Japanese"&gb$Estimator=="prob_ker_20"],y=0.74)+
  annotate("text",label="Igbo",x=gb$Surprisal[gb$Name=="Igbo"&gb$Estimator=="prob_ker_20"],y=0.78)+
  annotate("text",label="Mapudungun",x=gb$Surprisal[gb$Name=="Mapudungun"&gb$Estimator=="prob_ker_20"],y=0.64)+
  annotate("text",label="Yapese",x=gb$Surprisal[gb$Name=="Yapese"&gb$Estimator=="prob_ker_20"],y=0.77)+
  annotate("text",label="Yanomamö",x=gb$Surprisal[gb$Name=="Yanomamö"&gb$Estimator=="prob_ker_20"],y=0.78)+
  annotate("text",label="Tamil",x=gb$Surprisal[gb$Name=="Tamil"&gb$Estimator=="prob_ker_20"],y=0.742)+
  annotate("text",label="Standard Arabic",x=gb$Surprisal[gb$Name=="Standard Arabic"&gb$Estimator=="prob_ker_20"],y=0.739)+
  annotate("text",label="Russian",x=gb$Surprisal[gb$Name=="Russian"&gb$Estimator=="prob_ker_20"],y=0.65)+
  annotate("text",label="Hebrew",x=gb$Surprisal[gb$Name=="Hebrew"&gb$Estimator=="prob_ker_20"],y=0.7)+
  annotate("text",label="Standard Indonesian",x=gb$Surprisal[gb$Name=="Standard Indonesian"&gb$Estimator=="prob_ker_20"],y=0.76)+
  annotate("text",label="Hungarian",x=gb$Surprisal[gb$Name=="Hungarian"&gb$Estimator=="prob_ker_20"],y=0.77)

#########################################
## (5) Model rarity in terms of genealogical, areal covariates, and endangerement status
#########################################

# Load libraries
require(brms)

# Recode endangerment
gb<-gb %>%
  mutate(Endangerement=ifelse(aes %in% c("threatened","moribund","nearly_extinct"),"endangered",aes))

# Function that obtains a predictive model of surprisal
model_surprisal<-function(df) {
  m<-brm(Surprisal~
          (1|Macroarea)+
          (1|AUTOTYP_area)+
          (1|Family)+
           Endangerement,
        data=df,
        chains = 4,
        iter = 12000,
        warmup = 5000,
        cores = 4,
        control = list(adapt_delta=0.99),
        backend="cmdstanr")
  return(m)}

# From now on we centered our analyses on two estimates: LCA and kernel-20
gb_redux<-gb %>%
  filter(Estimator %in% c("prob_lca","prob_ker_20"))

# Get models
models_surprisal<-plyr::dlply(gb_redux,"Estimator",model_surprisal)

# Check summary
summary(model_surprisal_lca)
conditional_effects(model_surprisal_lca)

# Produce predictions and attach them to the empirical surprisals
gb_lca<-gb[gb$Estimator=="prob_lca"&!is.na(gb$Endangerement),]
predicted_surprisal_lca<-as.data.frame(predict(model_surprisal_lca))[,c("Estimate","Est.Error")]
gb_lca$surprisal_pred<-predicted_surprisal_lca$Estimate
gb_lca$surprisal_er<-predicted_surprisal_lca$'Est.Error'
gb_lca<-gb_lca %>%
  mutate(surprisal_z=(Surprisal-surprisal_pred)/surprisal_er)

# Plot this
gb_lca %>%
  ggplot(aes(x=surprisal_z,y=..density..,fill=stat(abs(x))>2))+
  geom_histogram(bins=50)+
  labs(x="Unexpected surprisal of Grambank language (z-score)",y="Density")+
  theme_minimal()+
  scale_fill_manual(values = c("deepskyblue","tomato"))+
  theme(legend.position="none")+
  theme(plot.background = element_rect(fill="white"))

ggsave("unexpected_surprisal.png",height=4,width=7)

# Check outliers
gb_lca %>%
  filter(abs(gb_lca$surprisal_z)>2) %>%
  select(Name,Endangerement,AUTOTYP_area,Macroarea,Family,surprisal_z) %>%
  arrange(desc(surprisal_z))

# Plot outliers
require(ggrepel)
rarity_ext %>%
  filter(abs(rarity_ext$surprisal_z)>2.3) %>%
  ggplot(aes(x=1,y=surprisal_z,label=Name,color=surprisal_z))+
  geom_text_repel(max.overlaps = 400,force_pull = 0.5)+
  theme_minimal()+
  theme(plot.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        axis.title=element_blank(),
        axis.text = element_blank())+
  scale_color_gradient2()


ggsave("unexpected_surprisal_lgs.png",height=7,width=6)

## Test
plyr::ddply(gb[gb$Estimator=="prob_lca"&!is.na(gb$aes),],"AUTOTYP_area",function(x) data.frame(E=sum(x$Endangerement=="endangered")/sum(x$Endangerement %in% c("not_endangered","endangered")),
                                                               U=mean(x$Surprisal),
                                                               N=nrow(x))) %>%
  ggplot(aes(x=U,y=E,label=AUTOTYP_area))+
  annotate("rect",xmin = 52.5,xmax = 60,ymin=0.75,ymax=1.05,fill="tomato",alpha=0.1)+
  geom_point(alpha=0.2,aes(size=5*N))+
  geom_text_repel(size=5)+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x="Mean unexpectedness in region",y="Proportion of threatened/moribund/near extinct lgs in region")

ggsave("endangerment_and_area.png")


plyr::ddply(gb[gb$Estimator=="prob_lca"&!is.na(gb$aes),],"Macroarea",function(x) data.frame(E=sum(x$Endangerement=="endangered")/sum(x$Endangerement %in% c("not_endangered","endangered")),
                                                                                               U=mean(x$Surprisal),
                                                                                               N=nrow(x))) %>%
  ggplot(aes(x=U,y=E,label=Macroarea))+
  geom_point(aes(size=log(N)),alpha=0.5)+
  geom_text()+
  theme_bw()+
  theme(legend.position = "none")

#########################################
## (6) Rarity and endangerment
#########################################

# Plot rarity
gb %>%
  filter(Estimator=="prob_lca") %>%
  mutate(Endagerement=fct_reorder(aes,Surprisal)) %>%
  ggplot(aes(y=Surprisal,x=Endagerement))+
  geom_boxplot()+
  labs(x="Unexpected surprisal of Grambank language",y="Density")+
  theme_minimal()+
  scale_fill_manual(values = c("deepskyblue","tomato"))+
  theme(legend.position="none")+
  coord_flip()

# Plot z-rarity
rarity_ext %>%
  mutate(Endagerement=fct_reorder(aes,surprisal_z)) %>%
  ggplot(aes(y=surprisal_z,x=Endagerement))+
  geom_boxplot()+
  labs(x="Surprisal of Grambank language",y="Density")+
  theme_minimal()+
  scale_fill_manual(values = c("deepskyblue","tomato"))+
  theme(legend.position="none")+
  coord_flip()

