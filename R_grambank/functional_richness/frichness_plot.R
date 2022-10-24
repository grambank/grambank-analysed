# Functional richness using PCoA

#script written by Simon Greenhill and Damián Blasi.

# (0) Read libraries
source("functional_richness/requirements_fr.R")

# (1) Read data and auxiliary functions (taken almost verbatim from Simon)

# see https://github.com/glottolog/glottolog/blob/master/config/aes_status.ini

#assign numbers to aes-status, splitting them into two categories: endangered and not endangered.
aes2numbers <- data.frame(
  aesid = c(2, 1, 1, 1, 1, 1),
  aes = c("not_endangered", "threatened", "shifting", "moribund", "nearly_extinct", "extinct")
)


gb <-read_tsv("output/GB_wide/GB_wide_imputed_binarized.tsv", show_col_types = F) 
  
languages <- read_tsv("output/non_GB_datasets/glottolog-cldf_wide_df.tsv", show_col_types = F)
languages <- languages[languages$Language_ID %in% gb$Language_ID, ]

# which languages have no endangerment?
languages[is.na(languages$aes), c("Name", "Glottocode", "Countries")]
# assume these are all safe, non-endangered languages
languages[is.na(languages$aes), ]$aes <- 'not_endangered'

languages <- languages %>% left_join(aes2numbers, by='aes')

table(languages$aes)

gb <- languages %>% dplyr::select('Language_ID', 'Family_ID') %>% left_join(gb, by="Language_ID")

autotyp <- read.delim('output/non_GB_datasets/glottolog_AUTOTYP_areas.tsv', na.strings="NA", header=TRUE) %>%
  dplyr::select(c("Language_ID", "AUTOTYP_area"))

# merge in regions
gb <- gb %>% left_join(autotyp, by = "Language_ID")

#gb.gower.mfd_fn <- "output/functional_richness/gb.gower.mfd.RDS"
#if(!file.exists(gb.gower.mfd_fn)) {
#  source("functional_richness/make_gower_mfd.R")
#}

#gb.dist <- readRDS(gb.gower.mfd_fn)

# drop languages we don't have
#gb <- gb[gb$Language_ID %in% attr(gb.dist, 'Labels'),]

# (2) Check AES endangerment categories
languages<-languages %>%
  filter(Language_ID %in% gb$Language_ID)

# Check distribution of endangerment
table(languages$aesid)
#1   2 
#944 560 

gb_dist_gower <- gb %>% 
  dplyr::select(-Language_ID,-Family_ID,-AUTOTYP_area) %>%
  cluster::daisy(metric = "gower", warnBin = F)

# (3) Run MDS on this data
mds_gb<-gb_dist_gower %>% 
  cmdscale(eig=TRUE, k=2)

# Add the MDS projections to the languages data frame
languages$mds1<-mds_gb$points[,1]
languages$mds2<-mds_gb$points[,2]

# Plot this

mds_plots<-languages %>%
  left_join(gb,by="Language_ID") %>%
  ggplot(aes(mds1,mds2,label=AUTOTYP_area))+
  geom_mark_hull(concavity = 10,
                 color="NA",
                 fill="#cae6d3",
                 alpha=1)+
  geom_point(color=ifelse(languages$aesid==1,"gray","black"),
             alpha=0.6)+
  facet_wrap(~AUTOTYP_area,ncol=4)+
  theme_bw() +
  theme(strip.text = element_text(size = 12),
    legend.position = "none",
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color="black",fill="NA"))

# Create subsets of language data
langsNEndangered<-languages$Language_ID[languages$aesid==2]

languagesNEndangered<-languages %>%
  filter(aesid==2) %>%
  dplyr::select(mds1,mds2)  %>% 
  as.data.frame() #making a tibble into a df only makes the next line not complain that you're setting rownames on a tibble

rownames(languagesNEndangered)<-langsNEndangered

languagesAll<-languages %>%
  dplyr::select(mds1,mds2) %>% 
  as.data.frame() #making a tibble into a df only makes the next line not complain that you're setting rownames on a tibble

rownames(languagesAll)<-languages$Language_ID

# Create required matrices for the functional richness analysis
areas_matrixNEndangered<-gb %>%
  filter(Language_ID %in% langsNEndangered) %>%
  make_group_matrix('AUTOTYP_area')

areas_matrixAll<-gb %>%
  make_group_matrix('AUTOTYP_area')

areas_matrixGlobalNEndangered<-gb %>%
  filter(Language_ID %in% langsNEndangered) %>%
  mutate(dummy="Global") %>%
  make_group_matrix('dummy')

areas_matrixGlobalAll<-gb %>%
  mutate(dummy="Global") %>%
  make_group_matrix('dummy')


# (3) Run Functional Richness with this binary distinction
frichNEndangered<-fundiversity::fd_fric(languagesNEndangered,
                              areas_matrixNEndangered,
                              stand=FALSE) %>%
  add_row(fundiversity::fd_fric(languagesNEndangered,
                                areas_matrixGlobalNEndangered,
                                stand=FALSE))

frichAll<-fundiversity::fd_fric(languagesAll,
                              areas_matrixAll,
                              stand=FALSE) %>%
  add_row(fundiversity::fd_fric(languagesAll,
                                areas_matrixGlobalAll,
                                stand=FALSE))

# (4) Create df and plots
frich_df<-frichNEndangered %>%
  full_join(frichAll,by="site")

colnames(frich_df)<-c("Area","NEndangered","All")
frich_df<-frich_df %>%
  pivot_longer(cols = c("All","NEndangered"),names_to = "Languages",values_to = "FRichness") %>%
  mutate(FRichness=ifelse(is.na(FRichness),0,FRichness))

# Transformations relevant to the plot
#frich_df<-frich_df %>%
#  mutate(Languages=fct_rev(Languages))

# Normalization
frich_df<-frich_df %>%
  mutate(FRichness=FRichness/max(frich_df$FRichness))

frich_df<-frich_df %>%
  mutate(Area=fct_reorder(Area,FRichness,max))

frich_df<-frich_df %>%
  group_by(Area) %>%
  mutate(Difference=FRichness-min(FRichness))

frich_df<-frich_df %>%
  group_by(Area) %>%
  mutate(FRichness=ifelse(Languages=="All",
                          Difference,
                          FRichness))

frich_plot<-frich_df %>%
  ggplot(aes(x=Area,y=FRichness,fill=Languages))+
  geom_bar(stat="identity",
           width=0.7,
           position="stack")+
  coord_flip()+
  labs(x="",y="Functional richness")+
  scale_fill_manual(values=c("#cae6d3","#71a674")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        legend.position = "none",
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        #        axis.text.y = element_text(margin = ggplot2::margin(r=-130),
        #                                   hjust=0,vjust=-0.7)
  )
  

mds_grid <- plot_grid(frich_plot,mds_plots,nrow=2, rel_heights = c(0.7, 1))

save_plot("output/functional_richness/frichness_paper.png",mds_grid,base_height = 13,base_width = 9)
