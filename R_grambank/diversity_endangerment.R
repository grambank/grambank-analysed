#script written by Damián Blasi
source("requirements.R")

#output directory
OUTPUTDIR <- file.path('.', "output", 'diversity_endangerment')
if (!dir.exists(OUTPUTDIR)){dir.create(OUTPUTDIR)}

#read in the endangerment status df
status_df<-read_tsv("output/non_GB_datasets/glottolog-cldf_wide_df.tsv", col_types = cols()) %>% 
  dplyr::select(Language_ID, endangerment_status = aes)

# Load the imputed and binarized dataset
gb <-read_tsv("output/GB_wide/GB_wide_imputed_binarized.tsv", col_types = cols()) %>% 
    mutate_at(vars(-Language_ID), as.factor) #make all GB feature cols be factors to ensure that they are understood as binary instead of interval scaled by cluster::daisy()

x <- assert_that(nrow(gb) == n_imputed, msg = "The number of languages has changed.")

# Get the names of the GB features
gb_features<-colnames(gb[,-1])

# Add endangerment status (and filter those for which we cannot determine it)
gb<-left_join(gb,status_df, by = "Language_ID") %>% 
  filter(!is.na(endangerment_status)) %>% 
  filter(endangerment_status != "moribund") %>% 
  filter(endangerment_status != "extinct") %>% 
  filter(endangerment_status != "nearly_extinct") 

# Collapse multiple endangerment statuses
gb$status_collapsed<-sapply(gb$endangerment_status,function(x)
  ifelse(is.na(x),NA,
         ifelse(x=="not_endangered","Not Endangered",
                ifelse(x=="shifting","Shifting", "Threatened"))))

# Define function that computes structural diversity 
compute_l<-function(df) {
    gb <- df %>% 
      dplyr::select(all_of(gb_features))
    return(mean(cluster::daisy(gb,metric="gower", warnBin = F)))}

# Define function that produces random samples out of the not endangered category (NEnd) and compute L and H 
rsampling<-function(DF,nsim)  {
  
  N<-plyr::ddply(DF[DF$status_collapsed!="Not Endangered",],"status_collapsed",nrow)
  ntypes<-nrow(N)

  df<-DF[DF$status_collapsed=="Not Endangered",]
  nrw<-nrow(df)
  output<-data.frame(sample_size= c(sapply(N$V1,function(x) rep(x,nsim)),N$V1),
                     l=rep(0,ntypes*(nsim+1)),
                     end_group= c(sapply(N$status_collapsed, function(x) rep(x,nsim)),N$status_collapsed),
                     type="random",
                     stringsAsFactors = FALSE)
  
  dummy<-0

  for(n in N$V1) {
    
    for(s in 1:nsim)  {
      
      sdf<-dplyr::sample_n(df,n)
      dummy<-dummy+1
      output[dummy,"l"]<-compute_l(sdf)
    }}
  
  for(type in N$status){
      dummy<-dummy+1
      edf<-DF[DF$status_collapsed==type,]
      output[dummy,"l"]<-compute_l(edf)
      output[dummy,"type"]<-"empirical"
  }
  
  return(output)}

# Produce reference distances
p1<-compute_l(gb[gb$Language_ID %in% c("yuec1235","mand1415"),gb_features])
# Egyptian Arabic / Standard Arabic
p2<-compute_l(gb[gb$Language_ID %in% c("egyp1253","stan1318"),gb_features])
# Estonian / Finnish
p4<-compute_l(gb[gb$Language_ID %in% c("esto1258","finn1318"),gb_features])

# Diverse sample: Latin, M Chinese, S Arabic, Maori, Evenki, Basque, Korean, Chippewa, K'iche, Zuni, Tupinamba, Hixkaryana, Tsez
q<-compute_l(gb[gb$Language_ID %in% c("lati1261","mand1415","stan1318","maor1246","even1259","basq1248","kore1280","chip1241","kich1262","zuni1245","tupi1273","hixk1239","dido1241"),gb_features])
ref_points<-c(p1,p2,p4,q)

cat("Running the sampling, this takes a while.\n")
# Run the analysis and produce plots
random_samples_rna<-rsampling(gb,1000)

# Add the reference comparisons
refs<-data.frame(l=ref_points,
                 case=c("Yue vs Mandarin Chinese",
                        "Standard vs Egyptian Arabic",
                        "Estonian vs Finnish",
                        "Highly diverse sample"),
                 type=c("similar","similar","similar","extreme"))

# Transform all distances into percentiles of the distribution of distances in all of GB
gb_dists<-cluster::daisy(gb[,gb_features],metric="gower", warnBin = F)

# Define function that turns distances into cumulative proportions
dist_to_cump<-function(x){sum(gb_dists<=x)/length(gb_dists)}

# Apply it to the distances
refs$l<-sapply(refs$l,dist_to_cump)
random_samples_rna$l<-sapply(random_samples_rna$l,dist_to_cump)

# Setup plot to display comparisons
new_plot_df<-random_samples_rna[random_samples_rna$type=="empirical",c("l","end_group","type")]
colnames(new_plot_df)[2]<-"case"
new_plot_df<-rbind(new_plot_df,refs)
new_plot_df<-rbind(new_plot_df,data.frame(l=mean(random_samples_rna[(random_samples_rna$type=="random")&(random_samples_rna$end_group=="Threatened"),"l"]),case="Not Endangered (control)",type="control"))
new_plot_df$x<-0

new_plot_df <- new_plot_df %>% 
  mutate(shape = ifelse(str_detect(case, "vs"), 17, 15))

#isolating specific values to variables that can be called in the plot easily
Shifting_value <- new_plot_df %>% 
  filter(case == "Shifting") %>% 
  dplyr::select(l) %>% 
  .[1,1]

Threatened_value <- new_plot_df %>% 
  filter(case == "Threatened") %>% 
  dplyr::select(l) %>% 
  .[1,1]

high_diverse_value <- new_plot_df %>% 
  filter(case== "Highly diverse sample")%>% 
  dplyr::select(l) %>% 
  .[1,1]

nend_value <- new_plot_df %>% 
  filter(case== "Not Endangered (control)")%>% 
  dplyr::select(l) %>% 
  .[1,1]

# Yue Chinese / Mandarin Chinese
chinese_pair_value <- new_plot_df %>% 
  filter(case== "Yue vs Mandarin Chinese")%>% 
  dplyr::select(l) %>% 
  .[1,1]

#standard and egyptian arabic
arabic_pair_value <- new_plot_df %>% 
  filter(case== "Standard vs Egyptian Arabic")%>% 
  dplyr::select(l) %>% 
  .[1,1]

#estonian and finnish
uralic_pair_value <-   new_plot_df %>% 
  filter(case== "Estonian vs Finnish")%>% 
  dplyr::select(l) %>% 
  .[1,1]

new_plot_df %>% 
  write_tsv(file = file.path(OUTPUTDIR, "diversity_endangerment_table.tsv"))

plot <- ggplot(data=new_plot_df,
       aes(y=l,x=x, color=factor(type),
           label=case)) +
  geom_vline(xintercept = 0,size=2,color="gray95")+
      geom_hline(yintercept = 0,size=2,color="gray95")+
  geom_segment(inherit.aes=FALSE,
               data=data.frame(y1=quantile(random_samples_rna[(random_samples_rna$type=="random")&(random_samples_rna$end_group=="Threatened"),"l"],probs=0.975),
                               y2=quantile(random_samples_rna[(random_samples_rna$type=="random")&(random_samples_rna$end_group=="Threatened"),"l"],probs=0.025)),
               aes(x=0,xend=0,y=y1,yend=y2),color="#1a8cff",size=10)+
  geom_point(size=5, shape =new_plot_df$shape)+
  annotate("text",label="Shifting",x=0.27,y=0.55,hjust=0,size=10,color="tomato2")+
  annotate("segment",x=0,xend=0.25,y=Shifting_value,yend=0.55,color="tomato2",size=1)+
  annotate("text",label="Highly diverse sample",x=0.27,y=high_diverse_value,hjust=0,size=10,color="maroon")+
  annotate("segment",x=0,xend=0.25,y=high_diverse_value,yend=high_diverse_value,color="maroon",size=1)+
  annotate("text",label="Threatened",x=0.27,y=0.51,hjust=0,size=10,color="tomato2")+
  annotate("segment",x=0,xend=0.25,y=Threatened_value,yend=0.51,color="tomato2",size=1)+
  annotate("text",label="Not Endangered",x=0.27,y=nend_value,hjust=0,size=10,color="#1a8cff")+
  annotate("segment",x=0,xend=0.25,y=nend_value,yend=nend_value ,color="#1a8cff",size=1)+
  annotate("text",label="Yue vs Mandarin Chinese",x=0.27,y= 0.16,hjust=0,size=10,color="#2eb873")+
  annotate("segment",x=0,xend=0.25,y= chinese_pair_value  ,yend=0.16,color="#2eb873",size=1)+
   annotate("text",label="Standard vs Egyptian Arabic",x=0.27,y= 0.12 ,hjust=0,size=10,color="#2eb873")+
  annotate("segment",x=0,xend=0.25,y=arabic_pair_value,yend= 0.12,color="#2eb873",size=1)+
   annotate("text",label="Estonian vs Finnish",x=0.27,y=0.08,hjust=0,size=10,color="#2eb873")+
  annotate("segment",x=0,xend=0.25,y=uralic_pair_value,yend=0.08 ,color="#2eb873",size=1)+
  labs(y="Cumulative structural diversity in Grambank")+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(size=19),
        axis.ticks.length.y =  unit(.5,"cm"),
        axis.line.y = element_blank(),
        axis.title.y = element_text(size=22, vjust = 2))+
  scale_x_continuous(limits=c(0,1))+
    scale_y_continuous(limits=c(0,0.65))+
    scale_colour_manual(values = c("control" = "#1a8cff", "empirical" = "tomato2","similar"="#2eb873","extreme"="maroon"))

# Save this
png(filename = file.path("OUTPUTDIR", "gb_diversity_endangerment.png"), height = 13, width = 10)
plot(plot)
x <- dev.off()

tiff(ilename = file.path("OUTPUTDIR", "gb_diversity_endangerment.tiff"), height = 13, width = 10)
plot(plot)
x <- dev.off()

cat("Done with the endangerment dissimilarity analysis.\n")