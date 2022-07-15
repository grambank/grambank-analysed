## This script visualises the output from a dual-process model for each Grambank feature
source('requirements.R')

OUTPUTDIR <- "output/spatiophylogenetic_modelling/spec_feature_plots/"
if(!dir.exists(OUTPUTDIR)){
  dir.create(OUTPUTDIR)
}

beep = 0

#description of featurs
GB_id_desc <- readr::read_tsv("output/GB_wide/parameters_binary.tsv", show_col_types = F) %>%
  dplyr::select(Feature_ID = ID, Grambank_ID_desc, Name)

source("spatiophylogenetic_modelling/analysis/INLA_parameters.R")

cat("\n###\nLoading covariance matrices...\n")

precision_matrices_fn <- "output/spatiophylogenetic_modelling/processed_data/precision_matrices_kappa_2_sigma_1.15.RDS"
if(!(file.exists(precision_matrices_fn))){
  source("spatiophylogenetic_modelling/analysis/make_precisionmatrices.R")}

precision_matrices = readRDS(precision_matrices_fn)
phylo_prec_mat = precision_matrices$phylogenetic_precision
spatial_prec_mat = precision_matrices$spatial_precision

color_vector <-c("#593d9cff", "#f68f46ff", "#d6d4d4")

#### Format Posterior Data ####
## Feature Metadata
parameters_binary <- read_csv("feature_grouping_for_analysis.csv", show_col_types = F)

## Read in model posteriors
model_output_files = list.files(path = "output/spatiophylogenetic_modelling/featurewise/",
                                pattern = "*.qs",
                                full.names = TRUE)

model_output_files <- model_output_files[!str_detect(model_output_files, "kapp")]

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

dual_hyperpar = lapply(model_output, function(m) {
  dd = apply(m[[4]]$hyper_sample, 2, function(x) median(log(x)))
  dd
})

dual_hyperpar <- do.call(rbind, dual_hyperpar) %>%
  as.data.frame()
dual_hyperpar$Feature_ID <- rownames(dual_hyperpar)

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
            cor = list(cor(cbind(phylogenetic, spatial))),
            median_phylogenetic = median(phylogenetic),
            median_spatial = median(spatial)) %>%
  arrange(desc(mean_phylogenetic)) %>%
  left_join(GB_id_desc, by = "Feature_ID")

three_most_phylo_features <- dual_summary %>%
  top_n(n = 3, wt = mean_phylogenetic) %>%
  arrange(desc(mean_phylogenetic)) %>%
  dplyr::select(Feature_ID) %>%
  as.matrix() %>%
  as.vector()

two_least_phylo_features <- dual_summary %>%
  top_n(n = -2, wt = mean_phylogenetic) %>%
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

glottolog_df <- read_tsv("output/non_GB_datasets/glottolog-cldf_wide_df.tsv", show_col_types = F) %>%
  mutate(Family_ID = ifelse(is.na(Family_ID), "Isolate", Family_ID))

#reading in GB
GB_fn <- "output/GB_wide/GB_cropped_for_missing.tsv"
if(!file.exists(GB_fn)){
  cat(paste0("Making GB data wide, binarising, cropping etc..\n"))
  source("make_wide.R")
  source("make_wide_binarized.R")
  source("impute_missing_values.R")
}
GB_df <- readr::read_tsv(file =   GB_fn,show_col_types = F) %>%
  inner_join(lgs_in_analysis, by = "Language_ID")

#values for clade labels
dat.offset = 1.02
ln.offset = 1.06
lab.offset = 1.11
fsize = 3
cex = 1.8

index <- 0

for(feature in c(three_most_phylo_features)) {

#  feature <- three_most_phylo_features[1]
  index <- index + 1
  cat(paste0("I'm at feature ", feature, " which is index ", index, ".\n"))

  feature_df <-GB_df %>%
    dplyr::select(Language_ID, all_of(feature)) #%>%
  #   rename(Feature = 2) %>%
  #   mutate(tip.color = as.character(Feature)) %>%
  #   mutate(tip.color = str_replace_all(tip.color, "1", color_vector[1])) %>%
  #   mutate(tip.color = str_replace_all(tip.color, "0", color_vector[2]))

  data <- GB_df %>%
    dplyr::select(Language_ID, all_of(feature))
  data$phylo_id = match(data$Language_ID, rownames(phylo_prec_mat))
  data$spatial_id = match(data$Language_ID, rownames(spatial_prec_mat))
  data$obs_id = 1:nrow(data)

  all_nodes <- rownames(phylo_prec_mat)
  all_node_df <- tibble(Language_ID = all_nodes,
                        spatial_id = NA,
                        obs_id = NA)

  all_node_df$phylo_id = match(all_node_df$Language_ID, rownames(phylo_prec_mat))
  all_node_df$what <- "predictions"
  data$what <- "observed_data"

  data <- bind_rows(data, all_node_df)

  formula <-  eval(substitute(expr = this_feature ~
                                f(spatial_id,
                                  model = "generic0",
                                  Cmatrix = spatial_prec_mat,
                                  hyper = pcprior) +
                                f(phylo_id,
                                  model = "generic0",
                                  Cmatrix = phylo_prec_mat,
                                  hyper = pcprior)  +
                                f(obs_id, model = "iid", hyper = obs_hyper),
                              list(this_feature=as.name(feature))))

  start_pars <- dual_hyperpar %>%
    filter(Feature_ID == feature) %>%
    dplyr::select(-Feature_ID) %>%
    unlist()

  dual_model <- INLA::inla(formula = formula,
                          family = "binomial",
                          control.predictor=list(link = 1, compute = TRUE),
                          control.mode = list(theta = start_pars, restart = FALSE),
                          data = data,
                          num.threads = 6)

  preds <- dual_model$summary.fitted.values$`0.5quant`[data$what == "predictions"]
  pred_df <- tibble(Language_ID = data$Language_ID[data$what == "predictions"], pred = preds)

  link_to_nodes <- makeNodeLabel(tree)
  link_to_nodes <- tibble(Language_ID = c(link_to_nodes$tip.label, link_to_nodes$node.label),
                          node_num = seq_len(Ntip(tree) + Nnode(tree)))

  pred_df <- pred_df %>%
    left_join(link_to_nodes) %>%
    arrange(node_num)

  #removing missing data
  #missing <- which(is.na(feature_df[,2]))

  #missing_language_IDs <- feature_df[which(is.na(feature_df[,2])), 1][[1]]
  #tree_feature <- ape::drop.tip(tree, missing_language_IDs)

  #feature_df <- feature_df[-missing,]
  # x <- feature_df[,2][[1]]
  # names(x) <- feature_df$Language_ID

  #generating plot title
  plot_title <- GB_id_desc %>%
    filter(Feature_ID == feature) %>%
    dplyr::select(Name) %>%
    as.matrix() %>%
    as.vector()

  plot_title <- paste0(feature, " ", plot_title)
  plot_title <- gsub('(.{1,60})(\\s|$)', '\\1\n', plot_title)

  filename <- paste(OUTPUTDIR, "/most_signal_phylo_",index, "_" , feature, ".tiff", sep = "")
  filename_png <- paste(OUTPUTDIR, "/most_signal_phylo_",index, "_" , feature, ".png", sep = "")

#running the contrasting algorithm reconstruction. Note: for the analysis we are using the tree with the original branch lengths even if we're visualizing using the imputed branch lengths.
  # asr_most_signal<- ape::ace(x = x, phy = tree_feature, method = "ML", type = "discrete", model = "ARD")

  pred_df %>%
    qs::qsave(paste0(OUTPUTDIR, "asr_object_",index , "_" , feature, ".qs"))

  cols <- colour_ramp(viridis(500))

  root <- plogis(dual_model$summary.fixed$mean[1])

  tips <- pred_df$pred[1:length(tree$tip.label)]
  names(tips) <- pred_df$Language_ID[1:length(tree$tip.label)]
  nodes <- c(root, pred_df$pred[-1:-length(tree$tip.label)])

  tiff(file = filename, width = 15.27, height = 15.69, units = "in", res = 400)

    par(mar = c(1,6,1,1),
        xpd = NA) ## xpd=NA prevents clipping in the margins (just turns off clipping entirely)

    p <- contMap(ladderize(tree), tips,
            anc.states = nodes,
            method = "user",
            type = "fan", ftype = "off",
            plot = FALSE)

    p <- setMap(p, viridis(500))
    plot(p, type = "fan", ftype = "off",
         mar = c(4, 4, 6, 4), legend = FALSE)

    add.color.bar(75, p$cols, title = "Predicted Probability\n",
                  prompt = FALSE, x = 0, y = 10, lwd = 12,
                  fsize = 1.5)

    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

    tip_points <- cbind(lastPP$xx[1:Ntip(tree)], lastPP$yy[1:Ntip(tree)],
                        feature_df[match(tree$tip.label, feature_df$Language_ID), feature])
    ## plot 'yes' values, blank is assumed 'no'
    points(tip_points[tip_points[ , 3] == 1, 1:2] * dat.offset, col = "black", pch = 19)
    ## plot missing data
    points(tip_points[is.na(tip_points[ , 3]), 1:2] * dat.offset, col = "grey", pch = 19)

    arc.cladelabels(text="Austronesian",node =   getMRCA(tree, tip = c("kana1286", "samo1305")), mark.node=FALSE,
                    ln.offset = ln.offset , lab.offset = lab.offset,fsize=fsize,orientation="curved", cex = cex)

    arc.cladelabels(text="Otomanguean",node =   getMRCA(tree, tip = c("mali1285", "yatz1235")), mark.node=FALSE,
                    ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

    arc.cladelabels(text="Uto-Aztecan",node =   getMRCA(tree, tip = c("hopi1249", "isth1240")), mark.node=FALSE,
                    ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

    arc.cladelabels(text="Nuclear Trans New Guinea",node =   getMRCA(tree, tip = c("gira1247", "domm1246")), mark.node=FALSE,
                    ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

    arc.cladelabels(text="Afro-Asiatic",node =   getMRCA(tree, tip = c("glav1244", "xamt1239")), mark.node=FALSE,
                    ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

    arc.cladelabels(text="Indo-European",node =   getMRCA(tree, tip = c("port1283", "mode1248")), mark.node=FALSE,
                    ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

    arc.cladelabels(text="Atlantic-Congo",node =   getMRCA(tree, tip = c("tswa1255", "noon1242")), mark.node=FALSE,
                    ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

    arc.cladelabels(text="Tibeto-Burman",node =   getMRCA(tree, tip = c("koir1240", "cent2004")), mark.node=FALSE,
                    ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

    arc.cladelabels(text="Austroasiatic",node =   getMRCA(tree, tip = c("seme1247", "aheu1239")), mark.node=FALSE,
                    ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

    arc.cladelabels(text="Uralic",node =   getMRCA(tree, tip = c("livv1243", "lule1254")), mark.node=FALSE,
                    ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

    title(plot_title, cex.main = cex)

  dev.off()

  png(file = filename_png, width = 15.27, height = 15.69, units = "in", res = 400)


    par(mar = c(1,6,1,1),
        xpd = NA) ## xpd=NA prevents clipping in the margins (just turns off clipping entirely)

    p <- contMap(ladderize(tree), tips,
            anc.states = nodes,
            method = "user",
            type = "fan", ftype = "off",
            plot = FALSE)

    p <- setMap(p, viridis(500))
    plot(p, type = "fan", ftype = "off",
         mar = c(4, 4, 6, 4),
         legend = FALSE)

    add.color.bar(75, p$cols, title = "Predicted Probability\n",
                  prompt = FALSE, x = 0, y = 10, lwd = 12,
                  fsize = 1.5)

    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

    tip_points <- cbind(lastPP$xx[1:Ntip(tree)], lastPP$yy[1:Ntip(tree)],
                        feature_df[match(tree$tip.label, feature_df$Language_ID), feature])
    ## plot 'yes' values, blank is assumed 'no'
    points(tip_points[tip_points[ , 3] == 1, 1:2] * dat.offset, col = "black", pch = 19)
    ## plot missing data
    points(tip_points[is.na(tip_points[ , 3]), 1:2] * dat.offset, col = "grey", pch = 19)

    arc.cladelabels(text="Austronesian",node =   getMRCA(tree, tip = c("kana1286", "samo1305")), mark.node=FALSE,
                    ln.offset = ln.offset , lab.offset = lab.offset,fsize=fsize,orientation="curved", cex = cex)

    arc.cladelabels(text="Otomanguean",node =   getMRCA(tree, tip = c("mali1285", "yatz1235")), mark.node=FALSE,
                    ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

    arc.cladelabels(text="Uto-Aztecan",node =   getMRCA(tree, tip = c("hopi1249", "isth1240")), mark.node=FALSE,
                    ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

    arc.cladelabels(text="Nuclear Trans New Guinea",node =   getMRCA(tree, tip = c("gira1247", "domm1246")), mark.node=FALSE,
                    ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

    arc.cladelabels(text="Afro-Asiatic",node =   getMRCA(tree, tip = c("glav1244", "xamt1239")), mark.node=FALSE,
                    ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

    arc.cladelabels(text="Indo-European",node =   getMRCA(tree, tip = c("port1283", "mode1248")), mark.node=FALSE,
                    ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

    arc.cladelabels(text="Atlantic-Congo",node =   getMRCA(tree, tip = c("tswa1255", "noon1242")), mark.node=FALSE,
                    ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

    arc.cladelabels(text="Tibeto-Burman",node =   getMRCA(tree, tip = c("koir1240", "cent2004")), mark.node=FALSE,
                    ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

    arc.cladelabels(text="Austroasiatic",node =   getMRCA(tree, tip = c("seme1247", "aheu1239")), mark.node=FALSE,
                    ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

    arc.cladelabels(text="Uralic",node =   getMRCA(tree, tip = c("livv1243", "lule1254")), mark.node=FALSE,
                    ln.offset= ln.offset,lab.offset = lab.offset,fsize = fsize,orientation="curved", cex = cex)

    title(plot_title, cex.main = cex)

  dev.off()


}

if(beep == 1){
  h_load("beepr")
  beep(2)
  }


#making map plots to show the feature that has the most spatial effect in the model

five_most_spatial_features <- dual_summary %>%
  top_n(n = 5, wt = mean_spatial) %>%
  arrange(desc(mean_spatial)) %>%
  dplyr::select(Feature_ID) %>%
  as.matrix() %>%
  as.vector()

#read in meta data
Language_meta_data_fn <- "output/non_GB_datasets/glottolog-cldf_wide_df.tsv"
if(!file.exists(Language_meta_data_fn)){
  source("make_glottolog-cldf_table.R")
}
Language_meta_data  <- read_tsv(Language_meta_data_fn, show_col_types = F) %>%
  dplyr::select(Language_ID, Longitude, Latitude)


df_for_maps <- GB_df %>%
  inner_join(Language_meta_data, by = "Language_ID") %>%
  mutate(Longitude = if_else(Longitude <= -25, Longitude + 360, Longitude)) #shift for pacific

#basemap
#rendering a worldmap that is pacific centered
world <- map_data('world', wrap=c(-25,335), ylim=c(-56,80), margin=T)
lakes <- map_data("lakes", wrap=c(-25,335), col="white", border="gray", ylim=c(-55,65), margin=T)

basemap <- ggplot(df_for_maps) +
  geom_polygon(data=world, aes(x=long, y=lat, group=group),
               colour="gray87",
               fill="gray87", size = 0.5) +
  geom_polygon(data=lakes, aes(x=long, y=lat, group=group),
               colour="gray87",
               fill="white", size = 0.3) +
  theme(
    # all of these lines are just removing default things like grid lines, axes etc
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  coord_map(projection = "vandergrinten", ylim=c(-56,67)) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

index <- 0

for(feature in five_most_spatial_features){

#  feature <- five_most_spatial_features[1]

  index <- index + 1
  #generating plot title
  plot_title <- GB_id_desc %>%
    filter(Feature_ID == feature) %>%
    dplyr::select(Name) %>%
    as.matrix() %>%
    as.vector()

  plot_title <- paste0(feature, " ", plot_title)
  plot_title <- gsub('(.{1,60})(\\s|$)', '\\1\n', plot_title)

  #filename
  filename <- paste(OUTPUTDIR, "/most_signal_spatial_",index, "_" , feature, ".png", sep = "")

  feature_df <-df_for_maps %>%
    filter(!is.na("Longitude")) %>%
    dplyr::select(Language_ID, all_of(feature), Longitude, Latitude) %>%
    rename(Feature = 2) %>%
    mutate(point.color = as.character(Feature)) %>%
    mutate(point.color = str_replace_all(point.color, "1", color_vector[1])) %>%
    mutate(point.color = str_replace_all(point.color, "0", color_vector[2])) %>%
    mutate(point.color = ifelse(is.na(point.color), color_vector[3], point.color)) %>%
    arrange(Feature)

  basemap +
    geom_jitter(data = feature_df, mapping = aes(x = Longitude, y = Latitude),  color = feature_df$point.color, alpha = 0.6, width = 3) +
    ggtitle(plot_title)

  ggsave(filename = filename, width = 9.3, height = 9.2, units = "in")

}

#comparing phylo and spatial effects with prop of 1s

comp_df <- GB_df %>%
  reshape2::melt(id.vars = "Language_ID") %>%
  group_by(variable) %>%
  summarise(one_prop = mean(value, na.rm = T)) %>%
  rename(Feature_ID = variable) %>%
  full_join(dual_summary, by = "Feature_ID") %>%
  dplyr::select(Feature_ID, one_prop, mean_phylogenetic, mean_spatial)

png("output/spatiophylogenetic_modelling/effect_plots/splom_dual.png")

psych::pairs.panels(comp_df[,-1],
method = "pearson", # correlation method
hist.col = "#a3afd1",# "#a9d1a3","",""),
density = TRUE,  # show density plots
ellipses = F, # show correlation ellipses
cex.labels= 1,
#           smoother= T,
cor=T,
lm=T,
ci = T, cex.cor = 0.9,stars = T)

x <- dev.off()

