source("requirements.R")
p_load("spdep")

#get the phylo.d table   
if (!file.exists("phylosig/pylosig_table.tsv")) { source("phylo_signal_features_d_stat.R") }		
phylo_d_table <- read_tsv("phylosig/pylosig_table.tsv", col_types = cols()) %>% 
  dplyr::rename(Feature_ID = Feature)

GB <- read_tsv("GB_wide/GB_wide_imputed_binarized.tsv", col_types = cols())


parameters <- read_csv("../../grambank_grambank/grambank/docs/feature_groupings/feature_grouping_for_analysis.csv", show_col_types = F, col_types = cols())


Language_meta_data <-  read_csv(GRAMBANK_LANGUAGES, col_types=LANGUAGES_COLSPEC) %>%		
  dplyr::select(Language_ID = Language_level_ID, Latitude, Longitude, Family_name) %>% 
  distinct(Language_ID, .keep_all = T) %>% 
  mutate(Family_name = ifelse(is.na(Family_name), "Isolate", Family_name))

GB_locations <- GB %>% 
  left_join(Language_meta_data, by = "Language_ID") %>% 
  dplyr::select(Language_ID,Longitude, Latitude) %>% 
  column_to_rownames("Language_ID") %>% 
  as.matrix()

GB_knn <- knearneigh(GB_locations, k = 5, longlat = T)

GB_nb <- knn2nb(GB_knn)

features <- GB %>% 
  dplyr::select(-Language_ID) %>% 
  colnames()

#empty df to bind to
df <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(df) <- c("Feature_ID","join_count_statistic", "p.value", "estimate_same_color","estimate_expectation", "estimate_variance", "alternative") 
df$Feature_ID <- as.character(df$Feature_ID)
df$join_count_statistic <- as.numeric(df$join_count_statistic)
df$p.value <- as.numeric(df$p.value)
df$estimate_same_color <- as.numeric(df$estimate_same_color)
df$estimate_expectation <- as.numeric(df$estimate_expectation)
df$estimate_variance <- as.numeric(df$estimate_variance)
df$alternative <- as.character(df$alternative)

index <- 0
for(feature in features){
  
  #feature <- features[1]
  
  cat(paste0("# Caclualting the join count statistic for ", feature, ". That means I'm ", round(index/length(features) * 100, 2), "% done.\n"))
  index <- index + 1 
  
gb_spec <- GB[[feature]] %>% as.factor()
names(gb_spec) <- GB$Language_ID
join_count_object <- spdep::joincount.test(gb_spec, listw = nb2listw(GB_nb, style = "B"))

df_spec <-  data.frame(
Feature_ID = feature,
join_count_statistic = join_count_object[[1]]$statistic[[1]],
p.value = join_count_object[[1]]$p.value[[1]],
estimate_same_color =join_count_object[[1]]$estimate[[1]],
estimate_expectation =join_count_object[[1]]$estimate[[2]],
estimate_variance =join_count_object[[1]]$estimate[[3]],
alternative = join_count_object[[1]]$alternative)

df <- df %>% 
  full_join(df_spec, by = c("Feature_ID", "join_count_statistic", "p.value", "estimate_same_color", "estimate_expectation", "estimate_variance", "alternative"))
}


joined_df <- df %>% 
  full_join(phylo_d_table, by = "Feature_ID")

png("temp_scripts_for_meetings/phylo_d_vs_join_count.png")
joined_df %>% 
  left_join(parameters, by = "Feature_ID") %>% 
  ggplot() +
  geom_point(aes(x = join_count_statistic, y = `D-estimate`, col = Main_domain)) +
  theme_classic()

x <- dev.off()



  


