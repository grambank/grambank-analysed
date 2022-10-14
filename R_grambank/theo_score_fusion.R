source("requirements.R")		

#script written by Hedvig Skirgård

OUTPUTDIR <- file.path('.', "output", "fusion_score")		
# create output dir if it does not exist.		
if (!dir.exists(OUTPUTDIR)) { dir.create(OUTPUTDIR) }		

GB_wide <- read_tsv(file.path("output", "GB_wide", "GB_wide_strict.tsv"), col_types=WIDE_COLSPEC)	

Language_meta_data <-read_csv(GRAMBANK_LANGUAGES, col_types=LANGUAGES_COLSPEC) %>%		
  dplyr::select(Language_ID = Language_level_ID, Family_name, Name, Longitude, Latitude, Macroarea) %>% 
  distinct(Language_ID, .keep_all = T)

#read in sheet with scores for whether a feature denotes fusion
GB_fusion_points <- data.table::fread(GRAMBANK_PARAMETERS,
                                                      encoding = 'UTF-8', 
                                                      quote = "\"", header = TRUE, 
                                                      sep = ",") %>% 
  dplyr::select(Parameter_ID = ID, Fusion = boundness,informativity) %>% 		
  mutate(Fusion = as.numeric(Fusion)) 

#remove features for which there is only a feature for the free or bound kind of marking, only keep those where there is one for each type of marking

# GB_fusion_points_only_with_alternatives <- GB_fusion_points %>% 
#   filter(!is.na(Fusion))	%>% 
#   filter(!is.na(informativity))	%>% 
#   group_by(informativity) %>% 
#   dplyr::summarise(count_informativity_categories = n()) %>% 
#   filter(count_informativity_categories > 1) %>% 
#   inner_join(GB_fusion_points,  by = "informativity") %>% 
#   dplyr::select(Parameter_ID, Fusion)

df_morph_count <- GB_wide %>%
  filter(na_prop <= 0.25 ) %>% #exclude languages with more than 25% missing data
  dplyr::select(-na_prop) %>% 
  reshape2::melt(id.vars = "Language_ID") %>% 
  dplyr::rename(Parameter_ID = variable) %>% 
  inner_join(GB_fusion_points, by = "Parameter_ID") %>% 
  filter(Fusion == 1) %>% 
  filter(!is.na(value)) %>% 
  mutate(value_weighted = if_else(Fusion == 0.5 & value == 1, 0.5, value )) %>% #replacing all instances of 1 for a feature that is weighted to 0.5 bound morph points to 0.5
#  mutate(value_weighted = if_else(Fusion == 0, abs(value-1), value_weighted)) %>% # reversing the values of the features that refer to free-standing markers 
  group_by(Language_ID) %>% 
  dplyr::summarise(mean_morph = mean(value_weighted))

df_morph_count  %>% 		
  write_tsv(file.path(OUTPUTDIR, "fusion_score.tsv"))		

GB_wide_morph_counts_summary <- df_morph_count %>% 		
  transform(bin = cut(mean_morph, 20)) %>% #grouping the mean scores into 20 bins		
  group_by(bin) %>% 		
  dplyr::summarize(n = n())		

# Tidy up the tick labels for the bin columns in the graph		
# as the plain output from cut() was a little bit messy.		
GB_wide_morph_counts_summary$bin_display <- GB_wide_morph_counts_summary$bin %>% 		
  str_replace_all("[\\(|\\[|\\]]", "") 		

GB_wide_morph_counts_summary <- GB_wide_morph_counts_summary %>% 		
  separate(col = bin_display, sep = ",", into = c("low", "high")) %>% 		
  mutate(low = round(as.numeric(low), digits = 2)) %>% 		
  mutate(high = round(as.numeric(high), digits = 2)) %>% 		
  unite(low, high, col = bin_display, sep = "-") 		

GB_morph_counts_plot <- GB_wide_morph_counts_summary %>% 		
  ggplot() + geom_bar(aes(x = factor(bin_display), y = n, fill = factor(bin_display)), stat = "identity", colour = "grey25", size = 1) +		
  theme_classic() +		
  theme(axis.title = element_text(size=18), 		
        axis.text.x = element_text(size = 18, angle = 90, hjust=0.95), 		
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.3, size = 16), 		
        legend.position = "None") +		
  labs(title="",		
       x ="Fusion scores", y = "Number of languages") + 		
  scale_fill_viridis(discrete = T)  		
plot(GB_morph_counts_plot)

ggsave(file.path(OUTPUTDIR, "GB_morph_score_mean_binned.tiff"), width = 7, height = 7)
