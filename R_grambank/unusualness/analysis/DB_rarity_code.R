source("requirements.R")
p_load("BayesLCA")

if(!file.exists("output/GB_wide/GB_wide_imputed_binarized.tsv")){
source("make_wide.R")
source("make_wide_binarized.R")
source("impute_missing_values.R")
}

if(!dir.exists("output/unusualness/")){
  dir.create("unusualness/")
}
if(!dir.exists("output/unusualness/plots")){
  dir.create("unusualness/plots")
}

# Grambank unusualness

# Load imputed binarized data
gb<-read.delim("output/GB_wide/GB_wide_imputed_binarized.tsv", sep = "\t")

# Analyze clusters given a fixed k (number of clusters)
k_range<-c(2:8)
LCA_clusters<-lapply(k_range,
                     function(k) blca.em(gb[,-1],
                                         k, 
                                         restarts=50))

# Get BICs
LCA_BIC<-data.frame(BIC=sapply(LCA_clusters, function(x) x$BIC),
                    k=k_range)

# Plot them
ggplot(LCA_BIC,aes(x=k,y=BIC))+
  geom_point()

# Attach rarities to each observation given a cluster
estimate_rarity_LCA<-function(lca,df) {
  p_item<-lca$itemprob
  class<-apply(lca$Z,
               1,
               function(x) which(x==max(x)))
  df$prob<-sapply(1:nrow(df),
                  function(y) estimate_prob_LCA(df[y,],p_item[class[y],]))
  df$surprisal<-(-log(df$prob))
  return(df)
}

# Attach probabilities to each observation given a cluster
estimate_prob_LCA<-function(obs,p){
  return(prod(obs*p+(1-obs)*(1-p)))}

# Apply
rarity_df<-estimate_rarity_LCA(LCA_clusters[[1]],gb[,-1])

rarity_df$Language_ID <- gb$Language_ID

write_tsv(rarity_df, "output/unusualness/tables/DB_rarity.tsv", na = "")

# Plot

db_rarity_plot <- ggplot(rarity_df,aes(x=surprisal,y=..density..))+
  geom_density(color="dodgerblue3")+
  geom_histogram(bins = 50,alpha=0.7,fill="dodgerblue1")+
  labs(x="Surprisal",y="Density")+
  theme_minimal()

png("output/unusualness/plots/DB_rarity_plot.png")
plot(db_rarity_plot)
x <- dev.off()

# Randomization test

