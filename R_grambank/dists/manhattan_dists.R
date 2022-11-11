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
  unite("order", "second",col = "bin_display", sep = "-", remove = F, )
 
#list of languages that have zero dist

glottolog_df <- read_tsv("output/non_GB_datasets/glottolog-cldf_wide_df.tsv", show_col_types = F) %>% 
  dplyr::select(Language_ID, Name, Macroarea, Family_ID, level)

glottolog_fam<- glottolog_df %>% 
  filter(level == "family") %>% 
  dplyr::select(Family_ID = Language_ID, Family_name = Name) 

glottolog_df <- glottolog_df %>% 
  left_join(glottolog_fam, by = "Family_ID")
  
zeros <- GB_dist_list %>% 
  filter(value ==0) %>% 
  dplyr::select(Var1) 

zeros_2 <- GB_dist_list %>% 
  filter(value ==0) %>% 
  dplyr::select(Var2) %>% 
  rename(Var1 = Var2)

zeros_df <- rbind(zeros, zeros_2) %>% 
  as.data.frame() %>%
  rename(Language_ID = "Var1") %>% 
  left_join(glottolog_df, by = "Language_ID") %>% 
  dplyr::select(Name, Language_ID, Macroarea, Family_name)

zeros_df %>% 
  write_tsv("output/dists/manhattan_dist_0.tsv", na = "")

GB_dist_list_binned$bin_display <- fct_reorder(GB_dist_list_binned$bin_display, GB_dist_list_binned$order)

GB_dist_list_binned %>% 
  ggplot() +
  geom_bar(aes(x = bin_display, y = n, fill = order, color = order), stat = "identity") +
  theme_classic() +
  theme(legend.position = "None") +
  theme(plot.margin = unit(c(0.5,0.5,0,0), "cm")) +
  xlab("Manhattan distances") +
  ylab("n") +
  viridis::scale_fill_viridis(discrete = F) +
  viridis::scale_color_viridis(discrete = F) +
  theme(axis.text.x = element_text(angle = 0)) +
  scale_x_discrete(breaks = GB_dist_list_binned$bin_display[c(1, 10, 20, 30, 40)])

ggsave("output/dists/plot_manhattan_dists.png", height = 4, width = 5)
