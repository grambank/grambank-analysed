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

surprisal_fn <- paste0(OUTPUTDIR_tables, "/surprisal.tsv")
if(!file.exists(surprisal_fn)){
  source("unusualness/analysis/get_unusualness_bayesLCA.R")
}
gb <- read_tsv(file = surprisal_fn, show_col_types = F) 

# Plot the pairwise relations between probabilities
gb %>%
  dplyr::select(c(Surprisal,Estimator,Language_ID)) %>%
  pivot_wider(id_cols = Language_ID,names_from = Estimator,values_from = Surprisal) %>% 
  GGally::ggpairs(columns=c("Bayesian LCA","Kernel 1","Kernel 5","Kernel 15", "Kernel 10","Kernel 20","Kernel 30","Kernel 40"),
                  mapping = aes(alpha = 0.1), 
                  upper = list(continuous = wrap("cor", method = "spearman"))
                  )+
  theme_classic()

ggsave(file.path(OUTPUTDIR_plots,"comparison_surprisals.png"),height=8,width = 10)

# Zooming into the LCA and the kernel-15 approaches, and highlighting Macroarea
gb %>%
  dplyr::select(c(Surprisal,Estimator,Name,Family_ID,Macroarea)) %>% 
  filter(Estimator == "Bayesian LCA" |
           Estimator == "Kernel 15") %>% 
  pivot_wider(id_cols = c(Name,Family_ID,Macroarea),names_from = Estimator,values_from = Surprisal) %>% 
  mutate(IE=ifelse(Family_ID=="indo1319","IE","nIE")) %>%
  ggplot(aes(x=`Bayesian LCA`,y=`Kernel 15`,label=Name,color=Macroarea))+
  geom_text(alpha=0.7)+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x="Surprisal based on LCA approximation",y="Unnormalized surprisal based on 15-kernel")

ggsave(filename = file.path(OUTPUTDIR_plots, "comparison_surprisals_pairwise_macroareas.png"),width=25,height=18)

gb[gb$aes %in% c("threatened","not_endangered","nearly_extinct","moribund"),] %>%
  dplyr::select(c(Surprisal,Estimator,Name,aes)) %>%
  filter(Estimator %in% c("Bayesian LCA","Kernel 15")) %>%
  pivot_wider(id_cols = c(Name,aes),names_from = Estimator,values_from = Surprisal) %>%
  mutate(Endangerment=ifelse(aes=="not_endangered","Safe","Not safe")) %>%
  ggplot(aes(x=`Bayesian LCA`,y=`Kernel 15`,label=Name,color=Endangerment))+
  geom_text(alpha=0.7)+
  theme_bw()+
  theme(legend.position = c(0.8,0.2))+
  labs(x="Surprisal based on LCA approximation",y="Unnormalized surprisal based on 15-kernel")

ggsave(file.path(OUTPUTDIR_plots,"comparison_surprisals_pairwise_endangerement.png"),width=25,height=18)

# Plot individual surprisals
gb %>%
  filter(Estimator=="Bayesian LCA") %>%
  ggplot(aes(x=Surprisal,y=..density..))+
  geom_density(color="dodgerblue3")+
  geom_histogram(bins = 50,alpha=0.7,fill="dodgerblue1")+
  labs(x="Surprisal of Grambank language",y="Density")+
  theme_minimal()+
  annotate("text",label="English",x=gb$Surprisal[gb$Name=="English"&gb$Estimator=="Bayesian LCA"],y=0.04)+
  annotate("text",label="French",x=gb$Surprisal[gb$Name=="French"&gb$Estimator=="Bayesian LCA"],y=0.041)+
  annotate("text",label="Mandarin",x=gb$Surprisal[gb$Name=="Mandarin Chinese"&gb$Estimator=="Bayesian LCA"],y=0.04)+
  annotate("text",label="Turkish",x=gb$Surprisal[gb$Name=="Turkish"&gb$Estimator=="Bayesian LCA"],y=0.041)+
  annotate("text",label="Japanese",x=gb$Surprisal[gb$Name=="Japanese"&gb$Estimator=="Bayesian LCA"],y=0.041)+
  annotate("text",label="Igbo",x=gb$Surprisal[gb$Name=="Igbo"&gb$Estimator=="Bayesian LCA"],y=0.042)+
  annotate("text",label="Mapudungun",x=gb$Surprisal[gb$Name=="Mapudungun"&gb$Estimator=="Bayesian LCA"],y=0.042)+
  annotate("text",label="Yapese",x=gb$Surprisal[gb$Name=="Yapese"&gb$Estimator=="Bayesian LCA"],y=0.041)+
  annotate("text",label="Yanomamö",x=gb$Surprisal[gb$Name=="Yanomamö"&gb$Estimator=="Bayesian LCA"],y=0.039)+
  annotate("text",label="Tamil",x=gb$Surprisal[gb$Name=="Tamil"&gb$Estimator=="Bayesian LCA"],y=0.042)+
  annotate("text",label="Standard Arabic",x=gb$Surprisal[gb$Name=="Standard Arabic"&gb$Estimator=="Bayesian LCA"],y=0.039)+
  annotate("text",label="Russian",x=gb$Surprisal[gb$Name=="Russian"&gb$Estimator=="Bayesian LCA"],y=0.039)+
  annotate("text",label="Hebrew",x=gb$Surprisal[gb$Name=="Hebrew"&gb$Estimator=="Bayesian LCA"],y=0.039)+
  annotate("text",label="Standard Indonesian",x=gb$Surprisal[gb$Name=="Standard Indonesian"&gb$Estimator=="Bayesian LCA"],y=0.039)+
  annotate("text",label="Hungarian",x=gb$Surprisal[gb$Name=="Hungarian"&gb$Estimator=="Bayesian LCA"],y=0.043)

gb %>%
  filter(Estimator=="Kernel 15") %>%
  ggplot(aes(x=Surprisal,y=..density..))+
  geom_density(color="dodgerblue3")+
  geom_histogram(bins = 50,alpha=0.7,fill="dodgerblue1")+
  labs(x="Surprisal of Grambank language",y="Density")+
  theme_minimal()+
  annotate("text",label="English",x=gb$Surprisal[gb$Name=="English"&gb$Estimator=="Kernel 15"],y=0.74)+
  annotate("text",label="French",x=gb$Surprisal[gb$Name=="French"&gb$Estimator=="Kernel 15"],y=0.72)+
  annotate("text",label="Mandarin",x=gb$Surprisal[gb$Name=="Mandarin Chinese"&gb$Estimator=="Kernel 15"],y=0.7)+
  annotate("text",label="Turkish",x=gb$Surprisal[gb$Name=="Turkish"&gb$Estimator=="Kernel 15"],y=0.78)+
  annotate("text",label="Japanese",x=gb$Surprisal[gb$Name=="Japanese"&gb$Estimator=="Kernel 15"],y=0.74)+
  annotate("text",label="Igbo",x=gb$Surprisal[gb$Name=="Igbo"&gb$Estimator=="Kernel 15"],y=0.78)+
  annotate("text",label="Mapudungun",x=gb$Surprisal[gb$Name=="Mapudungun"&gb$Estimator=="Kernel 15"],y=0.64)+
  annotate("text",label="Yapese",x=gb$Surprisal[gb$Name=="Yapese"&gb$Estimator=="Kernel 15"],y=0.77)+
  annotate("text",label="Yanomamö",x=gb$Surprisal[gb$Name=="Yanomamö"&gb$Estimator=="Kernel 15"],y=0.78)+
  annotate("text",label="Tamil",x=gb$Surprisal[gb$Name=="Tamil"&gb$Estimator=="Kernel 15"],y=0.742)+
  annotate("text",label="Standard Arabic",x=gb$Surprisal[gb$Name=="Standard Arabic"&gb$Estimator=="Kernel 15"],y=0.739)+
  annotate("text",label="Russian",x=gb$Surprisal[gb$Name=="Russian"&gb$Estimator=="Kernel 15"],y=0.65)+
  annotate("text",label="Hebrew",x=gb$Surprisal[gb$Name=="Hebrew"&gb$Estimator=="Kernel 15"],y=0.7)+
  annotate("text",label="Standard Indonesian",x=gb$Surprisal[gb$Name=="Standard Indonesian"&gb$Estimator=="Kernel 15"],y=0.76)+
  annotate("text",label="Hungarian",x=gb$Surprisal[gb$Name=="Hungarian"&gb$Estimator=="Kernel 15"],y=0.77)

# Produce predictions and attach them to the empirical surprisals
gb_lca<-gb[gb$Estimator=="Bayesian LCA"&!is.na(gb$aes),]
#predicted_surprisal_lca<-as.data.frame(predict(model_surprisal_lca))[,c("Estimate","Est.Error")]
#gb_lca$surprisal_pred<-predicted_surprisal_lca$Estimate
#gb_lca$surprisal_er<-predicted_surprisal_lca$'Est.Error'
#gb_lca<-gb_lca %>%
#  mutate(surprisal_z=(Surprisal-surprisal_pred)/surprisal_er)

# Plot this
# gb_lca %>%
#   ggplot(aes(x=surprisal_z,y=..density..,fill=stat(abs(x))>2))+
#   geom_histogram(bins=50)+
#   labs(x="Unexpected surprisal of Grambank language (z-score)",y="Density")+
#   theme_minimal()+
#   scale_fill_manual(values = c("deepskyblue","tomato"))+
#   theme(legend.position="none")+
#   theme(plot.background = element_rect(fill="white"))
# 
# ggsave(file.path(OUTPUTDIR_plots,"unexpected_surprisal.png"),height=4,width=7)

# # Check outliers
# gb_lca %>%
#   filter(abs(gb_lca$surprisal_z)>2) %>%
#   dplyr::select(Name,Endangerement,AUTOTYP_area,Macroarea,Family,surprisal_z) %>%
#   arrange(desc(surprisal_z))

# # Plot outliers
# unusualness_ext %>%
#   filter(abs(unusualness_ext$surprisal_z)>2.3) %>%
#   ggplot(aes(x=1,y=surprisal_z,label=Name,color=surprisal_z))+
#   geom_text_repel(max.overlaps = 400,force_pull = 0.5)+
#   theme_minimal()+
#   theme(plot.background = element_rect(fill="white"),
#         panel.grid.major = element_blank(),
#         axis.title=element_blank(),
#         axis.text = element_blank())+
#   scale_color_gradient2()
# 
# 
# ggsave(file.path(OUTPUTDIR_plots,"unexpected_surprisal_lgs.png"),height=7,width=6)

## Test
#plot_df <-  plyr::ddply(gb[gb$Estimator=="Bayesian LCA"&!is.na(gb$aes),],"AUTOTYP_area",function(x) data.frame(E=sum(x$Endangerement=="endangered")/sum(x$Endangerement %in% c("not_endangered","endangered")),
 #                                                                                              U=mean(x$Surprisal),
  #                                                                                             N=nrow(x))) 
  plot_df <- gb %>%
    dplyr::select(Language_ID, Estimator, aes, AUTOTYP_area, Surprisal) %>% 
  filter(Estimator == "Bayesian LCA") %>% 
  filter(!is.na(aes)) %>% 
  mutate(Endangerment=ifelse(aes=="not_endangered",0,1)) %>% 
  group_by(AUTOTYP_area) %>% 
  summarise(x = mean(Surprisal), 
            N = mean(Endangerment), 
            n_lgs = n()) 
  
  plot_df %>%
  ggplot(aes(y=N, x=x,label=AUTOTYP_area, size = n_lgs))+
  annotate("rect",ymin = 0.7,ymax = max(plot_df$N) + 0.1,xmin=60,xmax=max(plot_df$x)+2,fill="tomato",alpha=0.1)+
  geom_point(alpha = 0.6) +
  geom_label_repel(size=6, box.padding = 0.3, fill = "white",  
                   min.segment.length	= 0.1)+
  theme_classic(base_size = 13)+
  theme(legend.position = "none", 
        axis.title = element_text(size = 18))+
  labs(x="Mean unexpectedness in region",y="Proportion of threatened/moribund/near extinct lgs in region")

ggsave(file.path(OUTPUTDIR_plots,"endangerment_and_area.png"), width = 10, height = 10)


plyr::ddply(gb[gb$Estimator=="Bayesian LCA"&!is.na(gb$aes),],"Macroarea",function(x) data.frame(E=sum(x$Endangerement=="endangered")/sum(x$Endangerement %in% c("not_endangered","endangered")),
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
  filter(Estimator=="Bayesian LCA") %>%
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

