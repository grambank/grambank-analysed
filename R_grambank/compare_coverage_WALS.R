source("requirements.R")

#script written by Hedvig Skirgård

#this scripts compares the coverage of features per language in WALS and Grambank and produces a plot. The data is standardized such that dialects are aggregated into languages for both datasets to make a more accurate comparison. WALS languoids which lack glottocodes are not included in the comparison.

cat("Make plot to compare WALS and GB coverage.\n")

#getting table from Glottolog which contains information on the Language-level parent of dialects
glottolog_cldf_df <- read_tsv("non_GB_datasets/glottolog-cldf_wide_df.tsv",col_types = cols()) %>% 
  mutate(Language_level_ID = ifelse(is.na(Language_level_ID), Language_ID, Language_level_ID)) #making language-level entities their own parent, so that we can use this column for aggregation easier.

##grambank import and aggregation to language-level (i.e. merge dialects)
GB_wide <-read_tsv(file.path("GB_wide", "GB_wide_strict.tsv"), col_types = cols()) %>% 
  dplyr::select(Language_ID, GB_na_prop = na_prop) %>% 
  mutate(GB_na_prop_reverse = 1 - GB_na_prop) 

#creating column which indexes the languages in order of coverage, to be used for x-axis in plot. The most well-covered langauge in GB and WALS will differ, we are only interested for this plot in comparing the rank - not the specific languages themselves

GB_wide$x_ticks <- seq(1:nrow(GB_wide))

##WALS
#wals data needs to be extracted from a zipped file which is created temporarily and discarded once the necessary data is extracted.

wals_fn <- "../wals/cldf/"

cat("Fetching WALS data from ", wals_fn, ".\n")

WALS_data <- read_csv(file.path(wals_fn, "values.csv"),col_types = cols()) %>% 
  dplyr::select(ID = Language_ID, Value, Parameter_ID)

WALS_lgs_codes <- read_csv(file.path(wals_fn, "languages.csv"), col_types = cols()) 

#Important: languoids which do not have a glottocode in the WALS dataset are ignored since we can't confidently matche them to dialect parents or tell them apart from duplicates.

#counting how many languoids in WALS don't have glottocodes
wals_lgs_without_glottocodes <- WALS_data %>% 
  left_join(WALS_lgs_codes, by = "ID") %>% 
  distinct(ID, .keep_all = T) %>% 
  filter(is.na(Glottocode)) %>% 
  nrow()

cat(wals_lgs_without_glottocodes, "languoids in WALS did not have an associated Glottocode and where therefore ignored for further comaprison.\n")
  
#making WALS data wide instead of long
wals_wide <- WALS_data %>% 
  left_join(WALS_lgs_codes, by = "ID") %>% 
  filter(!is.na(Glottocode)) %>% 
  filter(!is.na(Value)) %>% 
  dplyr::select(Language_ID = Glottocode, Parameter_ID, Value) %>% 
  distinct(Language_ID, Parameter_ID, .keep_all = T) %>% 
  arrange(Parameter_ID) %>% 
  spread(key = Parameter_ID, value = Value, drop = FALSE) %>%  
  column_to_rownames("Language_ID")

#counting missing values in wals data
apply(wals_wide, 1, function(x) mean(is.na(x))) -> wals_wide$wals_na_prop

#wranging wals data into a similar shape to the GB data, i.e. if more than one languoid is coded for the same language (as defined by having the same language-level glottocode) the one with the least amount of missing data is kept. This means that dialects are removed and the dialect with the leats amount of missing data gets assigned the language-level glottocode. 
wals_wide_language_level <- wals_wide %>% 
  rownames_to_column("Language_ID") %>% 
  left_join(glottolog_cldf_df, by = "Language_ID") %>% 
  arrange(-wals_na_prop) %>% #sorting the dataframe my amount of missing data
  dplyr::select(-Language_ID) %>% 
  group_by(Language_level_ID) %>% 
  top_n(wt = -wals_na_prop, n = 1) %>% 
  dplyr::select(Language_ID = Language_level_ID, wals_na_prop) %>% 
  mutate(wals_na_prop_reverse = 1 - wals_na_prop) %>% 
  filter(!is.na(Language_ID)) %>% 
  arrange(wals_na_prop)

wals_wide_language_level$x_ticks <- seq(1:nrow(wals_wide_language_level))

#coverage plot with both distributions overlayed each other
WALS_GB_coverage_overlay <- ggplot() +
  geom_area(data = GB_wide, aes(x = x_ticks, y = GB_na_prop_reverse), stat = "identity", fill = "turquoise3", alpha = 7) + 
  geom_area(data = wals_wide_language_level, aes(x = x_ticks, y = wals_na_prop_reverse), stat = "identity", fill = "purple4", alpha = 0.8) +
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 25), 
        axis.text.y = element_text(hjust = 2, size = 25, angle = 30),
        axis.title.x = element_text(size = 35),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 35), 
        axis.line = element_blank()) +
  xlab("Languages") +
  ylab("Feature coverage") +
  theme(plot.title = element_text(size = 50, face = "bold", color = "black")) +
  scale_y_continuous(labels = percent)  +
  scale_x_continuous(breaks = c( 500, 1000, 1500, 2000)) +
  annotation_custom(grid.text("WALS", x=0.17,  y=0.17, gp=gpar(col = "white", fontsize=50, fontface="bold")))  +
  annotation_custom(grid.text("Grambank", x=0.40,  y=0.40, gp=gpar(col = "white", fontsize=50, fontface="bold")))  

tiff(file.path("coverage_plots", "WALS_GB_coverage_overlay.tiff"), height = 1100, width = 1600,    units = "px",  bg = "white")
plot(WALS_GB_coverage_overlay)
x <- dev.off()

png(file.path("coverage_plots", "WALS_GB_coverage_overlay.png"), height = 1100, width = 1600,    units = "px",  bg = "white")
plot(WALS_GB_coverage_overlay)
x <- dev.off()

cat("WALS-GB-comparison plot made.\n")
