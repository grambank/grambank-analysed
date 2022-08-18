source("spatiophylogenetic_modelling/analysis/featurewise/INLA_results_most_prep.R")

three_most_phylo_features <- dual_summary %>%
  top_n(n = 3, wt = mean_phylogenetic) %>%
  arrange(desc(mean_phylogenetic)) %>%
  dplyr::select(Feature_ID) %>%
  as.matrix() %>%
  as.vector()


index <- 0

for(feature in c(three_most_phylo_features)) {

#  feature <- three_most_phylo_features[1]
  index <- index + 1
  cat(paste0("I'm at feature ", feature, " which is index ", index, ".\n"))

  feature_df <-GB_df %>%
    dplyr::select(Language_ID, all_of(feature)) 
  
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
    dplyr::select(`Precision for spatial_id`, `Precision for phylo_id`) %>% 
    unlist()

  dual_model <- INLA::inla(formula = formula,
                          family = "binomial",
                          control.predictor=list(link = 1, compute = TRUE),
                          control.mode = list(theta = start_pars, restart = FALSE),
                          data = data,
                          num.threads = 6)
  
  qs::qsave(x = dual_model, file = paste0(OUTPUTDIR, "/INLA_obj_",index , "_" , feature, ".qs"))

  preds <- dual_model$summary.fitted.values$`0.5quant`[data$what == "predictions"]
  pred_df <- tibble(Language_ID = data$Language_ID[data$what == "predictions"], pred = preds)

  link_to_nodes <- makeNodeLabel(tree)
  link_to_nodes <- tibble(Language_ID = c(link_to_nodes$tip.label, link_to_nodes$node.label),
                          node_num = seq_len(Ntip(tree) + Nnode(tree)))

  pred_df <- pred_df %>%
    left_join(link_to_nodes, by = "Language_ID") %>%
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
    write_tsv(paste0(OUTPUTDIR, "/INLA_ASR_pred_df_",index , "_" , feature, ".tsv"))

output <- list(pred_df = pred_df, plot_title = plot_title, filename = filename, filename_png = filename_png, feature_df = feature_df, dual_model = dual_model)

qs::qsave(x = output, file = paste0(OUTPUTDIR, "/INLA_ASR_OBJ_", feature, ".qs"))  
}

if(beep == 1){
  h_load("beepr")
  beep(2)
  }



