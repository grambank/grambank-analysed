source("requirements.R")

OUTPUTDIR <- "output/dists/"
if (!dir.exists(OUTPUTDIR)) {dir.create(OUTPUTDIR)}

cfx_macroarea_list <-   read_tsv("output/dists/cfx_Macroarea_cut_off_0_list.tsv", show_col_types = F)

cfx_macroarea_list$Vars <- fct_reorder(cfx_macroarea_list$Vars, cfx_macroarea_list$Value_cfx)

cfx_macroarea_list %>% 
  ggplot(aes(x = Vars, y = Value_cfx)) +
  geom_bar(aes(fill = 1 - Value_cfx), stat = "identity") +
  theme_classic() +
  theme(text = element_text(angle = 70, hjust = 1, size = 20), 
        legend.position = "None", 
        axis.title.x = element_blank()) +
  scale_fill_viridis()+
  ylab("Cultural Fixation Score")

mean(cfx_macroarea_list$Value_cfx)

ggsave(filename = file.path(OUTPUTDIR, "cfx_barplot_macroarea.png"), height =  7.89, width =  8.61)

#heatmap

#macroarea
cfx_macroarea_list_separated <- cfx_macroarea_list %>%
  separate(Vars, into = c("Var1", "Var2"), sep = " - ")

#the dist list isn't symmetric, which isn't an issue really because we only want to display one half of the symmetric matrix, but currently it's not displaying one triangle consistently. So, we make it symmetric again and expicitly chop of one half
h_load("igraph")
g <- igraph::graph.data.frame(cfx_macroarea_list_separated, directed=FALSE)

cfx_dist_matrix_sym <- igraph::get.adjacency(g, attr="Value_cfx", sparse=FALSE)

cfx_dist_matrix_sym[upper.tri(cfx_dist_matrix_sym, diag = T)] <- NA

#ggplot pipe for heatmap
cfx_dist_matrix_sym %>% 
  reshape2::melt() %>% 
  filter(!is.na(value)) %>% 
  ggplot(aes(Var1, Var2)) +
  theme_classic() +
  geom_tile(aes(fill = value), color='white', linewidth = 0.6) +
  scale_fill_viridis(direction = -1) +
  geom_text(aes(label = round(value, 2)),color = "white", size = 8) +
  theme(axis.title = element_blank(),
        legend.position = "none", 
        axis.text = element_text(size = 14)) +
  scale_x_discrete(limits=rev)

ggsave(file.path(OUTPUTDIR, "cfx_heatmap_macroarea.png"), height =  7.89, width =  8.61)

#heatmap for autotyp-areas
cfx_autotyp_list <-   read_tsv("output/dists/cfx_AUTOTYP_area_cut_off_0_list.tsv", show_col_types = F) %>%
  separate(Vars, into = c("Var1", "Var2"), sep = " - ")

g <- igraph::graph.data.frame(cfx_autotyp_list, directed=FALSE)

cfx_dist_matrix_sym <- igraph::get.adjacency(g, attr="Value_cfx", sparse=FALSE)

cfx_dist_matrix_sym[upper.tri(cfx_dist_matrix_sym, diag = T)] <- NA

cfx_dist_matrix_sym %>% 
reshape2::melt() %>% 
  filter(!is.na(value)) %>% 
  ggplot(aes(Var1, Var2)) +
  theme_classic() +
  geom_tile(aes(fill = value), color='white', size = 0.6) +
  scale_fill_viridis(direction = -1) +
  geom_text(aes(label = round(value, 2)),color = "white", size = 3.5) +
  theme(axis.title = element_blank(),
        legend.position = "none", 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 50, hjust =1)) +
  scale_x_discrete(limits=rev)

ggsave(file.path(OUTPUTDIR, "cfx_heatmap_autotyp.png"), height =  10, width =  11)
