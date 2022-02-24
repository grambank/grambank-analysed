source("requirements.R")

#Script written by Hedvig Skirgård

cat("Matching all languages in Grambank to an AUTOTYP-area.\n")

#combining the tables languages and values from glottolog_df-cldf into one wide dataframe.
#this can be replaced with any list of Language_IDs, long and lat
glottolog_df <- read_tsv("non_GB_datasets/glottolog-cldf_wide_df.tsv",col_types = cols()) %>% 
dplyr::select(Language_ID, Longitude, Latitude)

##Adding in areas of linguistic contact from AUTOTYP

AUTOTYP_FN <- "https://raw.githubusercontent.com/autotyp/autotyp-data/master/data/csv/Register.csv"

cat("Fetching AUTOTYP data from", AUTOTYP_FN, ".\n")

AUTOTYP <- read_csv(AUTOTYP_FN  ,col_types = cols()) %>% 
  dplyr::select(Language_ID = Glottocode, Area, Longitude, Latitude) %>% 
  group_by(Language_ID) %>% 
  sample_n(1) #when a language is assigned to more than one area, pick randomly.

#This next bit where we find the autotyp areas of languages was written by Seán Roberts
# We know the autotyp-area of langauges in autotyp and their long lat. We don't know the autotyp area of languages in Glottolog. We also can't be sure that the long lat of languoids with the same glottoids in autotyp and glottolog_df have the exact identical long lat. First let's make two datasets, one for autotyp languages (hence lgs where we know the area) and those that we wish to know about, the Glottolog ones.

lgs_with_known_area <- as.matrix(AUTOTYP[!is.na(AUTOTYP$Area),c("Longitude","Latitude")])
rownames(lgs_with_known_area) <- AUTOTYP[!is.na(AUTOTYP$Area),]$Language_ID

known_areas <- AUTOTYP %>% 
  dplyr::filter(!is.na(Area)) %>% 
  dplyr::select(Language_ID, Area) %>% 
  distinct() %>% 
  dplyr::select(AUTOTYP_Language_ID = Language_ID, everything())

rm(AUTOTYP)

lgs_with_unknown_area <- as.matrix(glottolog_df[,c("Longitude","Latitude")])
rownames(lgs_with_unknown_area) <- glottolog_df$Language_ID

# For missing, find area of closest language

cat("Calculating the geographical distance between languages with known AUTOTYP-areas and those without a matched AUTOTYP-area.\n")
atDist <- rdist.earth(lgs_with_known_area,lgs_with_unknown_area, miles = F)

rm(lgs_with_known_area, lgs_with_unknown_area)

df_matched_up <- as.data.frame(unlist(apply(atDist, 2, function(x){names(which.min(x))})), stringsAsFactors = F) %>% 
  dplyr::rename(AUTOTYP_Language_ID = `unlist(apply(atDist, 2, function(x) {     names(which.min(x)) }))`)

cat("Matching languages without known AUTOTYP-area to the AUTOTYP-area of its closest neighbour with has a known AUTOTYP-area.\n")

glottolog_df_with_AUTOTYP <- df_matched_up %>% 
  tibble::rownames_to_column("Language_ID") %>%
  full_join(known_areas, by = "AUTOTYP_Language_ID") %>% 
  right_join(glottolog_df,  by = "Language_ID") %>% 
  dplyr::select(-AUTOTYP_Language_ID) %>% 
  dplyr::rename(AUTOTYP_area = Area) 

glottolog_df_with_AUTOTYP %>% 
  write_tsv("non_GB_datasets/glottolog_AUTOTYP_areas.tsv")