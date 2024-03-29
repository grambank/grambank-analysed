source("fun_def_h_load.R")
h_load("tidyverse")
source("global_variables.R")


#This script takes the values and languages tables from the glottolog cldf-release and combines then and transforms them to a wide data format from a long. It does not take into account the parameter or code tables.

OUTPUTDIR <- "output/non_GB_datasets/"
if (!dir.exists(OUTPUTDIR)){dir.create(OUTPUTDIR)}

#script written by Hedvig Skirgård

#finding the filenames for the two tables we are interested in, the language and value tables. The specific filenames can vary, so instead of identifying them via the filename we should check which of the tables conform to particular CLDF-standards and then take the filenames for the tables that conform to those standards from the meta-data json.

cldf_github_folder <- "../glottolog-cldf/cldf/"

cldf_json <- jsonlite::read_json(paste0(cldf_github_folder, "cldf-metadata.json"))

cat("Fetched glottolog-cldf data from", cldf_github_folder ,". \n")

#finding the fileanme for the relevant tables by checking which of the tables entered into the json meta data file conforms to a given cldf-standard and pulling the filename from there

index <- 0

#going over each table in the json and checking which one conforms and saving the index of that one to a separate variable. First: "values"
for (table in cldf_json$tables ) {
  
  index <- index +1
  
  if("dc:conformsTo" %in% names(table) & !is.null(table$`dc:conformsTo`)) { #not every table in a cldf dataset has this attribute, or it can be set to "null". So we gotta check that this is even a thing in this table before we proceed
    if(table$`dc:conformsTo` == "http://cldf.clld.org/v1.0/terms.rdf#ValueTable") {index_ValueTable  <- index}
  }}

#using the index we derived above, pulling out the filename for that table
values_fn_name <- cldf_json$tables[index_ValueTable][[1]]$url #not sure why this has the name "url" when it is just the filename but that is the way
values_csv_fn <- paste0(cldf_github_folder, values_fn_name) #creating the fn path

#doing the same thing for the second table we are interested in, "languages"

index <- 0

for (table in cldf_json$tables ) {
  
  index <- index +1
  
  if("dc:conformsTo" %in% names(table) & !is.null(table$`dc:conformsTo`)) { #not every table in a cldf dataset has this attribute, or it can be set to "null". So we gotta check that this is even a thing in this table before we proceed
    if(table$`dc:conformsTo` == "http://cldf.clld.org/v1.0/terms.rdf#LanguageTable") {index_LangaugeTable <- index}
  }}

#using the index we derived above, pulling out the filename for that table
language_fn_name <- cldf_json$tables[index_LangaugeTable][[1]]$url
languages_csv_fn <- paste0(cldf_github_folder, language_fn_name) #creating the fn path

#reading in data and making it wide
values <- readr::read_csv(values_csv_fn, na = c("","<NA>"), col_types = cols()) %>% 
  mutate(Value = ifelse(Parameter_ID == "aes", Code_ID, Value)) %>% 
  mutate(Value = ifelse(Parameter_ID == "med", Code_ID, Value)) %>% 
  reshape2::dcast(Language_ID ~ Parameter_ID, value.var = "Value") %>%  #making long data wide %>% 
  mutate(med = str_replace(med, "med-", "")) %>% 
  mutate(aes = str_replace(aes, "aes-", "")) 

languages <- readr::read_csv(languages_csv_fn, na = c("","<NA>"), col_types = cols())

#The languages-table from glottolog-cldf contains a paramter called "Language_ID" which is NOT the same as the parameter "Language_ID" in the values tables. This parameter is in fact the language leveled parent of a dialect. In order to avoid confusion, let's rename the parameter in the languages tables to the more transparent "Language_level_ID". This set-up first test if this is indeed a problem (i.e. if this is glottolog-cldf) and only does the renaming then.

if(str_detect(cldf_github_folder, "glottolog")) {
  languages <- languages %>% 
    dplyr::rename(Language_level_ID = Language_ID) %>% 
    dplyr::rename(Language_ID = ID)
} else{  languages <- languages %>% 
  dplyr::rename(Language_ID = ID)}

cldf_wide_df <- dplyr::full_join(values,languages, by = "Language_ID") 

write_tsv(cldf_wide_df, file = paste0(OUTPUTDIR ,"glottolog-cldf_wide_df.tsv"))

cat("Wrote output/non_GB_datasets/glottolog-cldf_wide_df.tsv\n", sep = "")