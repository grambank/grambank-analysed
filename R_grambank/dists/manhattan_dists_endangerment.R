source("requirements.R")

OUTPUTDIR <- "output/dists/"
if (!dir.exists(OUTPUTDIR)) {dir.create(OUTPUTDIR)}

#reading in GB
GB <- read.delim("output/GB_wide/GB_cropped_for_missing.tsv", sep ="\t")

GB_matrix <- GB %>%
  column_to_rownames("Language_ID") %>%
  as.matrix()

#calculating manhattan distances %>% 
GB_dist <- GB_matrix %>% 
  cluster::daisy(metric = "manhattan", warnBin = F) %>% 
  as.matrix()

#insert back the names
rownames(GB_dist) <- rownames(GB_matrix)
colnames(GB_dist) <- rownames(GB_matrix)

GB_dist[upper.tri(GB_dist, diag = T)] <- NA

#all
GB_dist_list <- GB_dist %>% 
  reshape2::melt() %>%
  filter(!is.na(value)) 

GB_dist_list$value %>% mean()

GB_dist_list %>% 
  group_by(value) %>% 
  summarise(n = n()) %>%
  ggplot() +
  geom_point(aes(x = value, y = n, col = value)) +
  theme_classic() +
  theme(axis.title = element_blank(), 
        legend.position = "None") 

ggsave("output/dists/plot_manhattan_dists.png", height = 4, width = 5)
