source("requirements.R")

Language_meta_data <-  read_tsv("output/non_GB_datasets/glottolog-cldf_wide_df.tsv", show_col_types = F) %>% 
  dplyr::select(Language_ID, level, Family_ID, Name, Macroarea) %>% 
  mutate(Family_ID = ifelse(is.na(Family_ID), "Isolate", Family_ID))

#sometimes languages in an autotyp-area are not classified all as the same macroarea. we take the macroarea per autotyp-area that is the most common. for example, N Coast Asia has 58 languages glottolog puts in Eurasia and 2 in North Americas. We assign N Coast Asia to Eurasia macroarea

autotyp_area <- read_tsv("output/non_GB_datasets/glottolog_AUTOTYP_areas.tsv", col_types = cols()) %>%
  dplyr::select(Language_ID, AUTOTYP_area)

macroarea_per_autotyp_area <- Language_meta_data %>% 
  left_join(autotyp_area, by = "Language_ID") %>% 
  group_by(AUTOTYP_area, Macroarea) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  full_join(autotyp_area, by = "AUTOTYP_area") %>% 
  group_by(Language_ID, AUTOTYP_area) %>% 
  slice_max(n) %>% 
  dplyr::select(Language_ID, AUTOTYP_area, Macroarea)
  
Language_meta_data <- Language_meta_data  %>% 
  dplyr::select(-Macroarea) %>% 
  left_join(macroarea_per_autotyp_area, by = "Language_ID") %>% 
  mutate(americas = ifelse(str_detect(Macroarea, "merica"), "americas", "not americas"))

left <- Language_meta_data %>% 
  dplyr::select(Var1 = AUTOTYP_area, americas_var1 = americas)

right <- Language_meta_data %>% 
  dplyr::select(Var2 = AUTOTYP_area, americas_var2 = americas)

df_long <- read_tsv("output/dists/cfx_AUTOTYP_area_cut_off_0_list.tsv", show_col_types = F) %>% 
  separate("Vars", into = c("Var1", "Var2"), sep = " - ") 

#recover symmetric distance matrix
h_load("igraph")
g <- igraph::graph.data.frame(df_long, directed=FALSE)

cfx_dist_matrix_sym <- igraph::get.adjacency(g, attr="Value_cfx", sparse=FALSE)

#calculating modularity

membership_vec <- g[8] %>%
  names() %>% 
  as.data.frame() %>% 
  rename("AUTOTYP_area" = ".") %>% 
  left_join(distinct(Language_meta_data, `AUTOTYP_area`, `americas`), by = "AUTOTYP_area") %>% 
  mutate(americas_num = ifelse(americas == "americas", 1, 2)) %>% 
  dplyr::select(americas_num) %>% 
  as.matrix() %>% 
  as.vector()

igraph::modularity(g, membership = membership_vec, weights = df_long$Value_cfx)

random_membership <- sample(1:2, size = 24, replace = T)

igraph::modularity(g, membership = random_membership, weights = df_long$Value_cfx)


#make network graph
h_load("qgraph")

cfx_dist_matrix_sym_inverse <- (1 - cfx_dist_matrix_sym)^2

labels <- colnames(cfx_dist_matrix_sym) %>% str_replace_all(" ", "\n")  %>% str_replace_all("N\n", "N ")%>% str_replace_all("S\n", "S ")
groups<-  colnames(cfx_dist_matrix_sym) %>% 
  as.data.frame() %>% 
  rename("AUTOTYP_area" = ".") %>% 
  left_join(distinct(Language_meta_data, `AUTOTYP_area`, `americas`), by = "AUTOTYP_area") %>% 
  mutate(americas_color = ifelse(americas == "americas", "orange", "turquoise3")) %>% 
  dplyr::select(americas_color) %>% 
  as.matrix() %>% 
  as.vector()

qgraph_plot <- qgraph(cfx_dist_matrix_sym_inverse, 
       layout='spring', 
       labels = labels, 
       label.scale = F, 
       shape = "rectangle", 
       node.width = 1.5, 
       node.height= 0.9, 
       maximum = 0.5,
       vTrans= 230, 
       color = groups,
       border.color = groups,
       edge.color = "purple4")

plot(qgraph_plot) 

png(filename = "output/dists/cfx_network_plot.png", width = 10, height = 10, units = "in", res = 400)
plot(qgraph_plot) 
x <- dev.off()


#h_load("dbscan")

#hdbscan_obj <- dbscan::hdbscan(cfx_dist_matrix_sym, minPts = 2)

#colnames(cfx_dist_matrix_sym) %>% 
#  as.data.frame() %>% 
#  mutate(HDBSCAN_cluster <- hdbscan_obj$cluster) %>% View()

#melting and testing dists


df_long_sided <- cfx_dist_matrix_sym %>% 
  reshape2::melt() %>% 
  left_join(left) %>% 
  left_join(right) %>% 
  distinct() %>% 
  filter(Var1 != Var2)

df_long_sided %>% 
  group_by(americas_var1, americas_var2) %>%
  summarise(mean= mean(value)) %>% 
  reshape2::dcast(americas_var1 ~ americas_var2) %>%
  column_to_rownames("americas_var1") %>% 
  as.matrix() %>% 
  round(2)

#MDS PLOT
mds <- MASS::isoMDS(cfx_dist_matrix_sym)

mds_df <- mds$points %>% 
  as.data.frame() %>% 
  rownames_to_column("Label") %>% 
  left_join(Language_meta_data, by = c("Label" = "AUTOTYP_area")) %>% 
  distinct()

mds_df %>% 
write_tsv("output/dists/cfx_mds_autotyp.tsv")


# The palette with black:
cbbPalette <- c( "#E69F00", "#009E73", "#F0E442", "#D55E00", "#CC79A7","#000000")

mds_plot <- mds_df %>% 
  distinct(Label, Macroarea, V1, V2, americas) %>% 
  ggplot(aes(x = V1, y = V2, color = americas,  group = americas)) +
  geom_point() +
  geom_density2d(bins = 12) +
  geom_label_repel(aes(label = Label)) +
  theme_minimal() +
  coord_fixed() +
  theme(legend.title = element_blank())

mds_plot

png("output/dists/MDS_plot_autotyp_area.png", width = 12, height = 12, units = "in", res = 300)
plot(mds_plot)
x <- dev.off()