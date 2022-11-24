#This is a script for imputing values with missForest in GB.

#script written by Hedvig Skirgård
source("requirements.R")
source('global_variables.R')

source("set_random_seed.R")

OUTPUTFILE <- file.path("output", "GB_wide", "GB_wide_imputed_binarized.tsv")

#CUTOFFS

# remove languages with less than LG_NA_PROP_CUTOFF data
# 0 = good (no missing data), 1 = bad (entirely missing)
LG_NA_PROP_CUTOFF = 0.25

# remove features present in less than FEAT_NA_PROP_CUTOFF languages
FEAT_NA_PROP_CUTOFF = 0.25

# For the imputation, the values necessarily need to be character strings, not numeric.
# The pipe below replaces the numeric datapoints with strings that are unable to be coerced
# into numeric (sometimes functions coerce into numeric by default).

#using the biniarised version
GB_wide_binarised <- read.delim(file.path("output" , "GB_wide", "GB_wide_binarized.tsv"), sep = "\t")

Binary_feat_n <- GB_wide_binarised %>% 
  dplyr::select(-na_prop, -Language_ID) %>% 
  ncol()

#In order to make sure the imputation is doing categorically and not for continuous integers the values are replaces by strings.
GB_wide_binarised_values_prepped_for_imputation <- GB_wide_binarised %>% 
  dplyr::select(-na_prop) %>% 
  reshape2::melt(id.vars = "Language_ID") %>% 
  mutate(value = str_replace_all(value, "0", "0 - absent")) %>%
  mutate(value = str_replace_all(value, "1", "1 - present")) %>%
  mutate(value = str_replace_all(value, "2", "2 - multistate")) %>%
  mutate(value = str_replace_all(value, "3", "3 - multistate")) %>%
  mutate(value = str_replace_all(value, "4", "4 - multistate")) %>%
  mutate(value = as.factor(value)) %>%
  dplyr::select(Language_ID, variable, value) %>%
  spread(key = variable, value, drop = FALSE) %>%
  left_join(dplyr::select(GB_wide_binarised, Language_ID, na_prop), by = "Language_ID") #adding back in the na_prop value per language

missing_value_percentage_before_imputation <- GB_wide_binarised_values_prepped_for_imputation  %>%
  dplyr::select(-na_prop, Language_ID) %>%
  is.na() %>%
  mean()

cat(sprintf("%0.2f%% of the data read was missing before imputation.\n",
  missing_value_percentage_before_imputation * 100
))

# In order to count the number of missing values per feature, not just per
# language, we need to transpose the dataframe.
#
# Note that Parameter_ID are row names here and therefore not counted.
GB_wide_transposed <- GB_wide_binarised_values_prepped_for_imputation %>%
  dplyr::select(-na_prop) %>%
  column_to_rownames("Language_ID") %>%
  t(.) %>%
  data.frame(stringsAsFactors = FALSE) 

#calculating the number of missing values per features as opposed to per language
GB_wide_transposed$na_prop_feat <- apply(GB_wide_transposed, 1, function(x) mean(is.na(x)))

# Subsetting the full dataset to contain only the languages that have 
# 0.76 percentage missing values and parameters that have less than 500
# missing languages
Features_to_keep <- GB_wide_transposed %>%
  rownames_to_column("Parameter_ID") %>%
  dplyr::select(Parameter_ID, na_prop_feat, everything()) %>%
  filter(na_prop_feat <=FEAT_NA_PROP_CUTOFF) %>%
  dplyr::select(Parameter_ID) %>%
  as.matrix() %>%
  as.vector()

Features_to_keep_n <- Features_to_keep %>% length()

GB_wide_cropped_for_imputation <- GB_wide_binarised_values_prepped_for_imputation %>%
  filter(na_prop <= LG_NA_PROP_CUTOFF ) %>%
  dplyr::select(Language_ID, all_of(Features_to_keep)) 

missing_value_percentage_after_cropping <- mean(is.na(GB_wide_cropped_for_imputation))
cat(sprintf(
  "The data has now been cropped, and there remains %0.2f%% missing values.\n",
  missing_value_percentage_after_cropping * 100
))
cat("These are the values that will be imputed by random forests.\n")

#Saving the language IDs separately, it's better to leave them out of the imputation workflow
GB_wide_cropped_lgs <- GB_wide_cropped_for_imputation$Language_ID

#for the MDS we want to use the cropped dataset, but not imputed. So we're writing that to file here before the imputation starts

GB_wide_cropped_for_imputation %>% 
  mutate_all(as.character) %>% 
  reshape2::melt(id = "Language_ID")  %>% #melting to make it easier to turn the values back to single numbers instead of strings
  mutate(value =  str_replace_all(value, "0 - absent", "0")) %>%
  mutate(value =  str_replace_all(value,  "1 - present", "1")) %>%
  dplyr::select(Language_ID, Parameter_ID = variable, value) %>%
  spread(key = Parameter_ID, value = value, drop = FALSE) %>% 
  write_tsv(file = file.path("output", "GB_wide", "GB_cropped_for_missing.tsv"))

#actual imputation
imputed_data <- GB_wide_cropped_for_imputation %>%
  dplyr::select(-Language_ID) %>%
  as.matrix() %>%
  data.frame() %>%
  mutate_all(as.factor) %>% 
  missForest() 

imputed_data$OOBerror

GB_imputed <- imputed_data$ximp
GB_imputed$Language_ID <- GB_wide_cropped_lgs

GB_imputed_num <- GB_imputed %>%
  mutate_all(as.character) %>% 
  reshape2::melt(id = "Language_ID")  %>%
  mutate(value =  str_replace_all(value, "0 - absent", "0")) %>%
  mutate(value =  str_replace_all(value,  "1 - present", "1")) %>%
  dplyr::select(Language_ID, Parameter_ID = variable, value) %>%
  spread(key = Parameter_ID, value = value, drop = FALSE)

GB_imputed_num %>%
  dplyr::select(Language_ID, everything()) %>%
  write_tsv(OUTPUTFILE)

cat("Wrote ", OUTPUTFILE, "\n")

cat("The imputation of missing values (excluding question marks) is done. The imputation was based on a binarised version of the dataset. ", round(mean(is.na(GB_wide_cropped_for_imputation))*100,2), "% of the data was imputed. 
    Note: we did not impute all missing values in the entire set. Before imputation the dataset was cropped to a subset where languages with more than ",   LG_NA_PROP_CUTOFF *100 ,"% missing values were removed and features with more than",  FEAT_NA_PROP_CUTOFF*100,"% missing languages were removed. This left us with", GB_imputed_num %>% nrow(), "languages and", Features_to_keep_n,"features (i.e. we cropped away",GB_wide_binarised %>% nrow - GB_imputed_num %>% nrow(), "languages and", Binary_feat_n - Features_to_keep_n, "features)

In the full dataset, there was",round(missing_value_percentage_before_imputation, 2)*100, "% of the data missing.

After this cropping the dataset contained ", round(mean(is.na(GB_wide_cropped_for_imputation))*100,2), "% missing data, this is what was imputed. The Out of Bag error is", round(imputed_data$OOBerror, 2), "

After imputation, the dataset contains",round(mean(is.na(GB_imputed_num)*100,2)), "% missing values and is ready for analysis that requires a complete dataset (PCA or otherwise).")


