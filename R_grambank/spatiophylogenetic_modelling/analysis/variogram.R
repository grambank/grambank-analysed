source('requirements.R')

#### Read Data ####

pca_filename = 'PCA/PCA_language_values.tsv'
pca_components = read_tsv(pca_filename, col_types = cols()) %>% 
  mutate(PC1 = scale(PC1), PC2 = scale(PC2), PC3 = scale(PC3)) 

languages <- read_csv(GRAMBANK_LANGUAGES, col_types=LANGUAGES_COLSPEC) %>%		
  dplyr::select(Language_ID = Language_level_ID, Family_name, Name, Longitude, Latitude, Macroarea) %>% 
  distinct(Language_ID, .keep_all = T) %>% 
  inner_join(pca_components, by = "Language_ID")

#### Spatial variogram ####
vario.pc1 = vario(n.bins = 15,
                   data =
                    languages[,c("Latitude", "Longitude", "PC1")],
                   type = "moran")

vario.pc2 = vario(n.bins = 15,
                  data =
                    languages[,c("Latitude", "Longitude", "PC2")],
                  type = "moran")

vario.pc3 = vario(n.bins = 15,
                  data =
                    languages[,c("Latitude", "Longitude", "PC3")],
                  type = "moran")

#### phylogenetic variogram ####
# Organise Tree & Data
tree_filename = 'spatiophylogenetic_modelling/processed_data/jaeger_pruned.tree'
phylogenetic_tree = read.tree(tree_filename)

# Subset PCA and languages to Jaeger set
pca_components = pca_components[pca_components$Language_ID %in% phylogenetic_tree$tip.label,]
languages = languages[languages$Language_ID %in% pca_components$Language_ID,]

# prune tree to dataset
taxa = pca_components$Language_ID
phylogenetic_tree = keep.tip(phylogenetic_tree, tip = taxa)

# Get Language levels
glottolog = read_tsv("non_GB_datasets/glottolog-cldf_wide_df.tsv", col_types = cols())

classification = str_split(glottolog$classification, "/")
classification = lapply(classification, function(x) x[c(1:3, length(x))]) %>%
  do.call(rbind, .) %>% 
  as.data.frame()

classification = classification[complete.cases(classification),]
colnames(classification) = c("top", "mid", "bottom", "taxa")

languages = left_join(languages, classification, 
                      by = c("Language_ID" = "taxa"))

languages$all = "all"
languages_complete = languages[!is.na(languages$top),]
languages_complete$PC1 = c(languages_complete$PC1)
languages_complete$all = factor("all")
languages_complete$top = factor(languages_complete$top)
languages_complete$mid = factor(languages_complete$mid)
languages_complete$bottom = factor(languages_complete$bottom)

languages_complete = as.data.frame(languages_complete)


### Calculate Moran's I for all distances
pc1.form <- PC1 ~ all/top/mid/bottom/Language_ID
pc1.phyvar <- correlogram.formula(pc1.form, 
                                  data = languages_complete) 

pc2.form <- PC2 ~ all/top/mid/bottom/Language_ID
pc2.phyvar <- correlogram.formula(pc2.form, 
                                  data = languages_complete) 

pc3.form <- PC3 ~ all/top/mid/bottom/Language_ID
pc3.phyvar <- correlogram.formula(pc3.form, 
                                  data = languages_complete) 


#### Make Plots ####
col_vector <- c("orange", "purple4", "turquoise3")

# Spatial Plot
plot_df = data.frame(pc1 = vario.pc1$vario,
                     pc2 = vario.pc2$vario,
                     pc3 = vario.pc3$vario)
plot_df$levels = vario.pc1$mean.bin.dist

plot_df = pivot_longer(plot_df, cols = c("pc1", "pc2", "pc3"))
plot_df$type = "spatial"

plot_df$name = factor(plot_df$name, levels = c("pc1", "pc2", "pc3"), 
                      labels = c("PC1", "PC2", "PC3"))

spatial_vario = ggplot(plot_df, aes(y = value, x = levels, col = name)) +
  geom_line() + 
  ylab("") + 
  xlab("Distance") + 
  ylim(c(-0.4, 1)) + 
  scale_colour_manual(values = col_vector) + 
  theme_classic() + 
    theme(legend.position = 'bottom', 
        legend.title = element_blank()) 

# Language family plot
phy_df = data.frame(
  pc1 = pc1.phyvar$obs,
  pc2 = pc2.phyvar$obs,
  pc3 = pc3.phyvar$obs,
  levels = pc1.phyvar$labels
)
phy_df$levels = factor(phy_df$levels, 
                       levels = c("Language_ID", 
                                  "bottom", "mid", 
                                  "top", "all"),
                       labels = c("Taxa", "LF - 2",
                                  "LF - 1", "Language Family", "All Languages"))

phy_long = pivot_longer(phy_df, 
                        cols = c("pc1", "pc2", "pc3"))
phy_long$type = "phylogeny"

phy_long$name = factor(phy_long$name, levels = c("pc1", "pc2", "pc3"), 
                      labels = c("PC1", "PC2", "PC3"))



phylo_vario = ggplot(phy_long, 
       aes(y = value, 
           x = levels, 
           group = name, 
           col = name)) +
  geom_line() + 
  ylim(c(-0.4, 1)) + 
  ylab("Moran's I") + 
  xlab("Levels") + 
  scale_colour_manual(values = col_vector) + 
  theme_classic() + 
  theme(legend.position = 'bottom', 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1)) 

# combine plots using patchwork
sidebyside = phylo_vario + spatial_vario + 
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom', legend.title = element_blank()) 

ggsave(sidebyside, 
       filename = 'spatiophylogenetic_modelling/figures/variogram_sidebyside.png')
