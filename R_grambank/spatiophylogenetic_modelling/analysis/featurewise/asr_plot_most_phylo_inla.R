## This script visualises the output from a dual-process model for each Grambank feature
source('requirements.R')

source("spatiophylogenetic_modelling/analysis/INLA_parameters.R")

color_vector <-c("#593d9cff", "#f68f46ff", "#F9F9F9")

#### Format Posterior Data ####
## Feature Metadata
parameters_binary <- read_csv("feature_grouping_for_analysis.csv", show_col_types = F)

## Read in model posteriors
model_output_files = list.files(path = "output/spatiophylogenetic_modelling/featurewise/",
                                pattern = "*.qs",
                                full.names = TRUE)

model_output = lapply(model_output_files, qread)
names(model_output) = basename(model_output_files) %>%
  tools::file_path_sans_ext(.)

## Extract the posterior distrubtions for each dual_process model from the hypersample
dual_posterior = lapply(model_output, function(m) {
  dd = m[[4]]$hyper_sample
  binomial_error = pi^2 / 3
  # Calculate h^2
  posterior = (1 / dd) / (rowSums(1 / dd) + 1 + binomial_error)
  
  posterior
})
# Convert this to a dataframe
dual_posterior = map_df(dual_posterior, ~as.data.frame(.x), .id="id")

colnames(dual_posterior) = c(
  "Feature_ID",
  "spatial",
  "phylogenetic"
)

# join feature metadata to posterior
dual_posterior = left_join(dual_posterior, parameters_binary, by ="Feature_ID")

# Summarise the posterior distributions
dual_summary = dual_posterior %>%
  group_by(Feature_ID) %>%
  summarise(mean_phylogenetic = mean(phylogenetic),
            mean_spatial = mean(spatial),
            error_phylogenetic = sd(phylogenetic),
            error_spatial = sd(spatial),
            domain = first(Main_domain),
            cor = list(cor(cbind(phylogenetic, spatial))))

five_most_phylo_features <- dual_summary %>% 
  top_n(n = 5, wt = mean_phylogenetic) %>% 
  arrange(desc(mean_phylogenetic)) %>% 
  dplyr::select(Feature_ID) %>% 
  as.matrix() %>% 
  as.vector()

#reading in tree
tree_fn <- "output/spatiophylogenetic_modelling/processed_data/EDGE_pruned_tree.tree"
if(!(file.exists(tree_fn))){
  source("spatiophylogenetic_modelling/processing/pruning_EDGE_tree.R")}
tree = read.tree(tree_fn)

#double check that subset to lgs in GB cropped dataset
tree <- ape::keep.tip(tree, lgs_in_analysis$Language_ID)

glottolog_df <- read_tsv("output/non_GB_datasets/glottolog-cldf_wide_df.tsv") %>% 
  mutate(Family_ID = ifelse(is.na(Family_ID), "Isolate", Family_ID))

h_load("randomcoloR")

  
tree_all_df <- tree$tip.label %>% 
  as.data.frame() %>% 
  rename("Language_ID" = ".") %>% 
    left_join(glottolog_df)

n <- length(unique(tree_all_df$Family_ID)) #counting how many distinct colors we need. Note that "NA" is also counted, and represents Isolates.

color_vector_families <- randomcoloR::distinctColorPalette(n)

tree_all_df$tip.color_vec <- color_vector_families[as.factor(tree_all_df$Family_ID)]

png(file = "output/coverage_plots/EDGE_tree_full.png", width = 15.27, height = 15.69, units = "in", res = 600)

plot.phylo(ladderize(tree  , right = F), 
           col="grey", 
           tip.color = tree_all_df$tip.color, 
           type = "fan", 
           cex = 0.3,
           label.offset = 0.02)

x <- dev.off()


#reading in GB
GB_fn <- "output/GB_wide/GB_cropped_for_missing.tsv"
if(!file.exists(GB_fn)){
  cat(paste0("Making GB data wide, binarising, cropping etc..\n"))
  source("make_wide.R")
  source("make_wide_binarized.R")
  source("impute_missing_values.R")
}  
GB_df <- readr::read_tsv(file =   GB_fn,show_col_types = F) %>% 
  inner_join(lgs_in_analysis)

GB_id_desc <- readr::read_tsv("output/GB_wide/parameters_binary.tsv", show_col_types = F) %>% 
  dplyr::select(Feature_ID = ID, Grambank_ID_desc)

index <- 0

for(feature in five_most_phylo_features) {
  
#  feature <- five_most_phylo_features[2]
  index <- index + 1

  feature_df <-GB_df %>% 
    dplyr::select(Language_ID, all_of(feature)) %>% 
    rename(Feature = 2) %>% 
    mutate(tip.color = as.character(Feature)) %>% 
    mutate(tip.color = str_replace_all(tip.color, "1", color_vector[1])) %>% 
    mutate(tip.color = str_replace_all(tip.color, "0", color_vector[2])) 
  
  #removing missing data
  missing <- which(is.na(feature_df[,2]))
  
  missing_language_IDs <- feature_df[which(is.na(feature_df[,2])), 1][[1]]
  tree_feature <- ape::drop.tip(tree, missing_language_IDs)
  
  feature_df <- feature_df[-missing,]
  x <- feature_df[,2][[1]]

  #generating plot title  
  plot_title <- GB_id_desc %>% 
    filter(Feature_ID == feature) %>% 
    dplyr::select(Grambank_ID_desc) %>% 
    as.matrix() %>% 
    as.vector() 
  
    filename <- paste("output/spatiophylogenetic_modelling/most_signal_", "_" , str_replace(plot_title, " ", "_"), ".tiff", sep = "")
  filename_png <- paste("output/spatiophylogenetic_modelling/most_signal_", "_" , str_replace(plot_title, " ", "_"), ".png", sep = "")
  
#running the contrasting algorithm reconstruction. Note: for the analysis we are using the tree with the original branch lengths even if we're visualizing using the imputed branch lengths.
  asr_most_signal<- ape::ace(x = x, phy = tree_feature, method = "ML", type = "discrete", model = "ARD")

  asr_most_signal %>% 
    qs::qsave(paste0("output/spatiophylogenetic_modelling/most_signal_", "_" , str_replace(plot_title, " ", "_"), ".qs"))
    
  tiff(file = filename, width = 15.27, height = 15.69, units = "in", res = 600)
  
  plot.phylo(ladderize(tree_feature  , right = F), 
             col="grey", 
             tip.color = feature_df$tip.color, 
             type = "fan", 
             cex = 0.25,
             label.offset = 0.02,main = plot_title)
  
  lastPP<-get("last_plot.phylo",env=.PlotPhyloEnv)
  ss<-unique(x) %>% sort(decreasing = T)
  par(fg="black")
  colors<-setNames(color_vector[1:length(ss)],ss)
#  add.simmap.legend(colors=colors,
#                    vertical=T,
#                    x=--1,
#                    y=-1,
#                    prompt=F)
  
  nodelabels(node=1:tree_feature$Nnode+Ntip(tree_feature),
             pie=asr_most_signal$lik.anc,
             piecol=color_vector, cex = 0.3)
  
  x <- dev.off()
  
  
  png(file = filename_png, width = 15.27, height = 15.69, units = "in", res = 400)
  
  plot.phylo(ladderize(tree_feature  , right = F), 
             col="grey", 
             tip.color = feature_df$tip.color, 
             type = "fan", 
             cex = 0.25,
             label.offset = 0.02,main = plot_title)
  
  lastPP<-get("last_plot.phylo",env=.PlotPhyloEnv)
  ss<-unique(x) %>% sort(decreasing = T)
  par(fg="black")
  colors<-setNames(color_vector[1:length(ss)],ss)

  #  add.simmap.legend(colors=colors,
  #                  vertical=T,
  #                  x=--1,
  #                  y=-1,
  #                  prompt=F)
  
  nodelabels(node=1:tree_feature$Nnode+Ntip(tree_feature),
             pie=asr_most_signal$lik.anc,
             piecol=color_vector, cex = 0.3)
  
  x <- dev.off()
  
  
}

