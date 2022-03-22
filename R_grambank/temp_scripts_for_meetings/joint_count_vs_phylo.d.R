source("requirements.R")
p_load("spdep")

#get the phylo.d table   
if (!file.exists("phylosig/pylosig_table.tsv")) { source("phylo_signal_features_d_stat.R") }		
phylo_d_table <- read.delim("phylosig/pylosig_table.tsv", sep = "\t") %>% 
  dplyr::rename(Feature_ID = Feature)

GB <- read.delim("GB_wide/GB_wide_imputed_binarized.tsv", sep = "\t")

parameters <- read.delim("feature_grouping_for_analysis.csv", sep = ",")


Language_meta_data <-  read.delim(GRAMBANK_LANGUAGES, sep = ",") %>%		
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
df <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(df) <- c("Feature_ID","jtot_z_value") 
df$Feature_ID <- as.character(df$Feature_ID)
df$jtot_z_value <- as.numeric(df$jtot_z_value)

index <- 0
for(feature in features){
  
  #feature <- features[1]
  
  cat(paste0("# Caclualting the join count statistic for ", feature, ". That means I'm ", round(index/length(features) * 100, 2), "% done.\n"))
  index <- index + 1 
  
gb_spec <- GB[[feature]] %>% as.factor()
names(gb_spec) <- GB$Language_ID
join_count_object <- spdep::joincount.multi(fx = gb_spec, listw = nb2listw(GB_nb, style = "B"))

jtot_z_value<- join_count_object["Jtot", "z-value"] 

df_spec <-  data.frame(
Feature_ID = feature,
jtot_z_value = jtot_z_value)

df <- df %>% 
  full_join(df_spec, by = c("Feature_ID", "jtot_z_value"))
}

cat("I'm a 100% done! Go me!")


joined_df <- df %>% 
  full_join(phylo_d_table, by = "Feature_ID") %>% 
  mutate(jtot_z_value = scale(jtot_z_value))

png("temp_scripts_for_meetings/phylo_d_vs_join_count.png")
joined_df %>% 
  left_join(parameters, by = "Feature_ID") %>% 
  ggplot(aes(x = jtot_z_value, y = `D-estimate`)) +
  geom_point(aes(col = Main_domain)) +
  theme_classic() +
  geom_smooth()

x <- dev.off()



  


