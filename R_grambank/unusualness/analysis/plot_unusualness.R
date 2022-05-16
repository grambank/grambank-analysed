#This script requires both get_unusualness_bayesLCA and predict surprisal to have be run prior.

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

surprisal_fn <- paste0(OUTPUTDIR_tables, "surprisal.tsv")
if(!file.exists(surprisal_fn)){
  source("unusualness/analysis/get_unusualness_bayesLCA.R")
}
gb <- read_tsv(file = surprisal_fn)

# Plot the pairwise relations between probabilities
gb %>%
  select(c(Surprisal,Estimator,Language_ID)) %>%
  pivot_wider(id_cols = Language_ID,names_from = Estimator,values_from = Surprisal) %>%
  GGally::ggpairs(columns=c("Bayesian LCA","Kernel 1","Kernel 5","Kernel 10","Kernel 20","Kernel 30","Kernel 40"),
                  mapping = aes(alpha = 0.1))+
  theme_minimal()

ggsave("comparison_surprisals.png",height=8,width = 10)



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
unusualness_ext %>%
  filter(abs(unusualness_ext$surprisal_z)>2.3) %>%
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
## (6) unusualness and endangerment
#########################################

# Plot unusualness
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

# Plot z-unusualness
unusualness_ext %>%
  mutate(Endagerement=fct_reorder(aes,surprisal_z)) %>%
  ggplot(aes(y=surprisal_z,x=Endagerement))+
  geom_boxplot()+
  labs(x="Surprisal of Grambank language",y="Density")+
  theme_minimal()+
  scale_fill_manual(values = c("deepskyblue","tomato"))+
  theme(legend.position="none")+
  coord_flip()

