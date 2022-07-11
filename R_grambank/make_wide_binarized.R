source("requirements.R")
#GB contains a small set of multistate features. They can be binarised, but they need to be done so in a particular way. This code renders a appropriately binarised version of the dataset.

GB_wide_strict <- read.delim(file.path("output", "GB_wide", "GB_wide_strict.tsv"), sep = "\t")

multistate_features <- read_csv(GRAMBANK_CODES, col_types=CODES_COLSPEC)  %>% 
  unite(possible_values_split, Name,  Description, sep = ": ") %>% 
  group_by(Parameter_ID) %>% 
  mutate(possible_values = paste0(possible_values_split, collapse = ", ")) %>% 
  distinct(Parameter_ID, possible_values) %>% 
  filter(str_detect(possible_values, "2")) %>% 
  dplyr::select(Parameter_ID) %>% 
  as.matrix() %>% 
  as.vector()

GB_wide_strict_multi <- GB_wide_strict %>% 
  column_to_rownames("Language_ID") %>% 
  dplyr::select(contains(multistate_features))

#GB024 multistate 1; Num-N; 2: N-Num; 3: both.
if("GB024" %in% colnames(GB_wide_strict_multi)){
  GB_wide_strict_multi$GB024a <- if_else(GB_wide_strict_multi$GB024 == "1"|GB_wide_strict_multi$GB024 == "3", "1", ifelse(GB_wide_strict_multi$GB024 == "2", "0", NA)) 
  
  GB_wide_strict_multi$GB024b <- if_else(GB_wide_strict_multi$GB024 == "2"|GB_wide_strict_multi$GB024 == "3", "1", ifelse(GB_wide_strict_multi$GB024 == "1", "0", NA)) 
}

#GB025 multistate 1: Dem-N; 2: N-Dem; 3: both.
if("GB025" %in% colnames(GB_wide_strict_multi)){
  GB_wide_strict_multi$GB025a <- if_else(GB_wide_strict_multi$GB025 == "1"|GB_wide_strict_multi$GB025 == "3", "1", ifelse(GB_wide_strict_multi$GB025 == "2", "0", NA)) 
  
  GB_wide_strict_multi$GB025b <- ifelse(GB_wide_strict_multi$GB025 == "2"|GB_wide_strict_multi$GB025 == "3", "1", ifelse(GB_wide_strict_multi$GB025 == "1", "0", NA)) 
}

#GB065 multistate 1:Possessor-Possessed; 2:Possessed-Possessor; 3: both
if("GB065" %in% colnames(GB_wide_strict_multi)){
  GB_wide_strict_multi$GB065a <- if_else(GB_wide_strict_multi$GB065 == "1"|GB_wide_strict_multi$GB065 == "3", "1", ifelse(GB_wide_strict_multi$GB065 == "2", "0", NA)) 
  
  GB_wide_strict_multi$GB065b <- if_else(GB_wide_strict_multi$GB065 == "2"|GB_wide_strict_multi$GB065 == "3", "1", ifelse(GB_wide_strict_multi$GB065 == "1", "0", NA)) 
}

#GB130 multistate 1: SV; 2: VS; 3: both
if("GB130" %in% colnames(GB_wide_strict_multi)){
  GB_wide_strict_multi$GB130a <- if_else(GB_wide_strict_multi$GB130 == "1"|GB_wide_strict_multi$GB130 == "3", "1", ifelse(GB_wide_strict_multi$GB130 == "2", "0", NA)) 
  
  GB_wide_strict_multi$GB130b <- if_else(GB_wide_strict_multi$GB130 == "2"|GB_wide_strict_multi$GB130 == "3", "1", ifelse(GB_wide_strict_multi$GB130 == "1", "0", NA)) 
}

#GB193 multistate 0: they cannot be used attributively, 1: ANM-N; 2: N-ANM; 3: both.
if("GB193" %in% colnames(GB_wide_strict_multi)){
  GB_wide_strict_multi$GB193a <- if_else(GB_wide_strict_multi$GB193 == "1"|GB_wide_strict_multi$GB193 == "3", "1", ifelse(GB_wide_strict_multi$GB193 == "2"|GB_wide_strict_multi$GB193 == "0", "0", NA)) 

  GB_wide_strict_multi$GB193b <- if_else(GB_wide_strict_multi$GB193 == "2"|GB_wide_strict_multi$GB193 == "3", "1", ifelse(GB_wide_strict_multi$GB193 == "1"|GB_wide_strict_multi$GB193 == "0", "0", NA)) 
}
#GB203 multistate 0: no UQ, 1: UQ-N; 2: N-UQ; 3: both.
if("GB203" %in% colnames(GB_wide_strict_multi)){
  GB_wide_strict_multi$GB203a <- if_else(GB_wide_strict_multi$GB203 == "1"|GB_wide_strict_multi$GB203 == "3", "1", ifelse(GB_wide_strict_multi$GB203 == "2"|GB_wide_strict_multi$GB203 == "0", "0", NA)) 

  GB_wide_strict_multi$GB203b <- if_else(GB_wide_strict_multi$GB203 == "2"|GB_wide_strict_multi$GB203 == "3", "1", ifelse(GB_wide_strict_multi$GB203 == "1"|GB_wide_strict_multi$GB203 == "0", "0", NA)) 
}

GB_wide_strict_multi_only_binarized <- GB_wide_strict_multi %>% 
  dplyr::select(-all_of(multistate_features)) %>% 
  rownames_to_column("Language_ID")

stopifnot(all(!multistate_features %in% colnames(GB_wide_strict_multi_only_binarized)))

output_path <- file.path("output", "GB_wide", "GB_wide_binarized.tsv")

cat("Wrote", output_path, "\n")

#making parameters table_binary 
GB_wide_strict %>% 
  dplyr::select(-all_of(multistate_features)) %>% 
  full_join(GB_wide_strict_multi_only_binarized,by = "Language_ID") %>% 
  dplyr::select(Language_ID, na_prop, everything()) %>% # reordering columns for inspection convenience
  write_tsv(output_path)

parameters_multi <- data.table::fread(GRAMBANK_PARAMETERS ,
                                      encoding = 'UTF-8', 
                                      quote = "\"", header = TRUE, 
                                      sep = ",") 

Parameter_desc_binary <- tibble(
  ID = c(
      "GB024a", "GB024b",
      "GB025a", "GB025b",
      "GB065a", "GB065b",
      "GB130a","GB130b",
      "GB193a","GB193b",
      "GB203a", "GB203b"
  ),
  Grambank_ID_desc = c(
    "GB024a NUMOrder_Num-N",
    "GB024b NUMOrder_N-Num",  
    "GB025a DEMOrder_Dem-N",
    "GB025b DEMOrder_N-Dem",
    "GB065a POSSOrder_PSR-PSD",
    "GB065b POSSOrder_PSD-PSR",
    "GB130a IntransOrder_SV",
    "GB130b IntransOrder_VS",
    "GB193a ANMOrder_ANM-N",
    "GB193b ANMOrder_N-ANM",
    "GB203a UQOrder_UQ-N",
    "GB203b UQOrder_N-UQ"
),
Name = c("Is the order of the numeral and noun Num-N?", 
        "Is the order of the numeral and noun N-Num?",
        "Is the order of the adnominal demonstrative and noun Dem-N?",
        "Is the order of the adnominal demonstrative and noun N-Dem?",
        "Is the pragmatically unmarked order of adnominal possessor noun and possessed noun PSR-PSD?",
        "Is the pragmatically unmarked order of adnominal possessor noun and possessed noun PSD-PSR?",
        "Is the pragmatically unmarked order of S and V in intransitive clauses S-V?",
        "Is the pragmatically unmarked order of S and V in intransitive clauses V-S?",
        "Is the order of the adnominal property word (ANM) and noun ANM-N?",
        "Is the order of the adnominal property word (ANM) and noun N-ANM?",
        "Is the order of the adnominal collective universal quantifier (UQ) and noun UQ-N?",
        "Is the order of the adnominal collective universal quantifier (UQ) and noun N-QU?" ),

`OV vs VO types (excl affixes)`= c(
"OV",
  "VO",
  "OV",
  "VO",
  "OV",
  "VO",
  NA,
  NA,
  "OV",
  "VO",
  "OV",
  "VO"),
`OV VO score for counting`= c(
  0,
  1,
  0,
  1,
  0,
  1,
  NA,
  NA,
  0,
  1,
  0,
  1
), 
Binary_Multistate = c("binarised","binarised","binarised","binarised","binarised","binarised","binarised","binarised","binarised","binarised","binarised","binarised")) %>% 
  
  full_join(parameters_multi, by = c("ID", "Grambank_ID_desc"))

Parameter_desc_binary %>% 
  mutate(Binary_Multistate= ifelse(ID %in% multistate_features, "Multi", Binary_Multistate)) %>% 
  mutate(Binary_Multistate = ifelse(is.na(Binary_Multistate), "Binary", Binary_Multistate)) %>% 
  dplyr::select(-Description) %>% 
  write_tsv(file.path("output/GB_wide/parameters_binary.tsv"))

