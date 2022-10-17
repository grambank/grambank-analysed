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
        legend.position = "None") +
    xlab("Manhattan distance") +
    ylab("n")

GB_dist_list_binned <- GB_dist_list %>% 
  transform(bin = cut(value, 40)) %>% #grouping the mean scores into 20 bins		
  group_by(bin) %>% 		
  dplyr::summarize(n = n())	%>% 
  separate(bin, into = c("order", "second"), sep = ",", remove = F) %>% 
  mutate(order = str_replace_all(order, "\\(", "") %>% as.numeric() %>% round(digits = 0)) %>% 
  mutate(second = str_replace_all(second, "[\\(|\\[|\\]]", "") %>% as.numeric() %>% round(digits = 0)) %>% 
  unite("order", "second",col = "bin_display", sep = " - ", remove = F, )
 

GB_dist_list_binned$bin_display <- fct_reorder(GB_dist_list_binned$bin_display, GB_dist_list_binned$order)

GB_dist_list_binned %>% 
  ggplot() +
  geom_bar(aes(x = bin_display, y = n, fill = order), stat = "identity") +
  theme_classic() +
  theme(legend.position = "None") +
  xlab("Manhattan distance") +
  ylab("n") +
  viridis::scale_fill_viridis(discrete = F) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1)) 

ggsave("output/dists/plot_manhattan_dists.png", height = 4, width = 5)