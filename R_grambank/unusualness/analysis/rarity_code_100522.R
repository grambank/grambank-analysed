# Analysis of Grambank rarity
# Code by D E Blasi

#########################################
## (1) Setup
#########################################

# Load pkgs
source("requirements.R")

# Load imputed binarized GB data
gb<-read_tsv("output/GB_wide/GB_wide_imputed_binarized.tsv", col_types = cols())

# Load Glottolog families, AUTOTYP areas, and AES endangerment
# Add autotyp areas to language tibble
autotyp_areas <- read_tsv(file.path("output", "non_GB_datasets", "glottolog_AUTOTYP_areas.tsv"), col_types = cols()) %>%
  dplyr::select(Language_ID, AUTOTYP_area) 

gb_ext<-read_tsv("output/non_GB_datasets/glottolog-cldf_wide_df.tsv", col_types = cols()) %>%
  dplyr::select(Language_ID,Macroarea,Family_ID,aes,level,Name) %>%
  full_join(autotyp_areas, by = "Language_ID"  ) %>% 
  right_join(gb[,"Language_ID"], by = "Language_ID"  )

# If no language family is assigned to a language then it's an isolate - in which case, use the name of the language as that of the family
gb_ext<-gb_ext %>%
  mutate(Family=ifelse(is.na(Family_ID)&level=="language",Language_ID,Family_ID)) %>%
  select(-Family_ID)

# Check the coverage of the covariates
sum(!is.na(gb_ext$Macroarea))/nrow(gb_ext)
sum(!is.na(gb_ext$AUTOTYP_area))/nrow(gb_ext)
sum(!is.na(gb_ext$Family))/nrow(gb_ext)
sum(!is.na(gb_ext$aes))/nrow(gb_ext)

#########################################
## (2) Determine clusters 
#########################################

# In this section we sort all GB features into independent bundles such that: 
# I. they are highly correlated within them and 
# II. reasonably independent from features in other bundles

# Produce a hierarchical clustering of the features
gower_gb<-cluster::daisy(t(gb[,-1]), metric = "gower")
hier_gb<-hclust(gower_gb, method = "complete")

# Plot the resulting clusters and highlight the four largest clusters
plot(hier_gb)
rect.hclust(hier_gb , k = 4, border = 2:6)

# 
n_clusters<-4 # Clearly three large groups
hier_classes <- cutree(hier_gb, 
                       k = n_clusters)

#########################################
## (3) Induce a hierarchical probability distribution with latent classes
#########################################

# Here we take the clusters from the previous section, and within each we induce an optimal number of latent classes using
# a Bayesian implementation of LCA

# Load library
require(BayesLCA)

# Set a range of latent clusters to explore for each cluster
k_range<-c(2:4)

# Function that yields the BIC of latent classes assignments for feature cluster i
k_choose<-function(i) {
  
    LCA_clusters<-lapply(k_range,
                     function(k) blca.em(gb[,names(hier_classes[hier_classes==i])],
                                         k, 
                                         restarts=150))

    # Get BICs
    LCA_BIC<-data.frame(BIC=sapply(LCA_clusters, function(x) x$BIC),
                    k=k_range,
                    cluster_n=i)
    
    return(LCA_clusters[[which(LCA_BIC$BIC==min(LCA_BIC$BIC))]])}

# Apply this to all clusters
LCA_clusters<-lapply(c(1:n_clusters),function(x) k_choose(x))

#########################################
## (4) Assign a probability to each GB entry
#########################################

# Based on the latent classes obtained before, attach rarities to each observation given a cluster
estimate_rarity_LCA<-function(lca,df) {
  
  p_item<-lca$itemprob
  
  class<-apply(lca$Z,
               1,
               function(x) which(x==max(x)))
  
  prob_assignment<-apply(lca$Z,
               1,
               function(x) max(x))
  
  prob<-sapply(1:nrow(df),
                  function(y) estimate_prob_LCA(df[y,],p_item[class[y],]))
  
  return(data.frame(prob=prob,
                    surprisal=(-log(prob)),
                    assignment=prob_assignment,
                    lg=1:length(prob)))
  }

# Attach probabilities to each observation given a cluster
estimate_prob_LCA<-function(obs,p){
  return(prod(obs*p+(1-obs)*(1-p)))}

# Apply
cluster_list<-c(1:n_clusters)
names(cluster_list)<-cluster_list
  
rarity_df<-plyr::ldply(cluster_list,
                       function(x) estimate_rarity_LCA(LCA_clusters[[x]],gb[,names(hier_classes[hier_classes==x])]),
                       .id="Cluster") %>%
  plyr::ddply("lg",function(x) data.frame(prob=prod(x$prob),
                                    surprisal=sum(x$surprisal)))

# Reattach glottocodes
rarity_df$Language_ID<-gb$Language_ID
rarity_df<-rarity_df %>%
  select(-lg)
  
# Plot
ggplot(rarity_df,aes(x=surprisal,y=..density..))+
  geom_density(color="dodgerblue3")+
  geom_histogram(bins = 50,alpha=0.7,fill="dodgerblue1")+
  labs(x="Surprisal of Grambank language",y="Density")+
  theme_minimal()

#########################################
## (5) Model rarity in terms of genealogical and areal covariates
#########################################

# Load libraries
require(brms)

# Enrich GB 
rarity_ext<-rarity_df %>%
  left_join(gb_ext)

# Set up a simple model
model_rarity<-brm(surprisal~
                    (1|Macroarea)+
                    (1|AUTOTYP_area)+
                    (1|Family),
                  data=rarity_ext,
                  chains = 4,
                  iter = 12000,
                  warmup = 5000,
                  cores = 4)

# Check summary
summary(model_rarity)

# Produce predictions and attach them to the empirical surprisals
predicted_surprisal<-as.data.frame(predict(model_rarity))[,c("Estimate","Est.Error")]
rarity_ext$surprisal_pred<-predicted_surprisal$Estimate
rarity_ext$surprisal_er<-predicted_surprisal$'Est.Error'
rarity_ext<-rarity_ext %>%
  mutate(surprisal_z=(surprisal-surprisal_pred)/surprisal_er)

# Plot this
rarity_ext %>%
  ggplot(aes(x=surprisal_z,y=..density..,fill=stat(abs(x))>2))+
  geom_histogram(bins=50)+
  labs(x="Unexpected surprisal of Grambank language (z-score)",y="Density")+
  theme_minimal()+
  scale_fill_manual(values = c("deepskyblue","tomato"))+
  theme(legend.position="none")+
  theme(plot.background = element_rect(fill="white"))

ggsave("unexpected_surprisal.png",height=4,width=7)

# Check outliers
rarity_ext %>%
  filter(abs(rarity_ext$surprisal_z)>2) %>%
  select(Name,aes,AUTOTYP_area,Macroarea,Family,surprisal_z) %>%
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


#########################################
## (6) Rarity and endangerment
#########################################

# Plot rarity
rarity_ext %>%
  mutate(Endagerement=fct_reorder(aes,surprisal)) %>%
  ggplot(aes(y=surprisal,x=Endagerement))+
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

