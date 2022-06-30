#defining expected number of languages
n_total <- 2467
n_imputed <- 1509
n_overlap_imputed_and_jaeger_tree <- 1265
n_overlap_imputed_and_EDGE_tree <- 1403

GRAMBANK_LANGUAGES <- file.path("../grambank", "cldf", "languages.csv")
GRAMBANK_VALUES <- file.path("../grambank", "cldf", "values.csv")
GRAMBANK_PARAMETERS <- file.path("../grambank", "cldf", "parameters.csv")
GRAMBANK_CODES <- file.path("../grambank", "cldf", "codes.csv")

source("fun_def_h_load.R")
h_load(pkg = c("readr"))

# The columns specifier for readr to parse ../cldf/values.csv
VALUES_COLSPEC <- c(
  ID = col_character(),
  Language_ID = col_character(),
  Parameter_ID = col_character(),
  Value = col_character(),
  Code_ID = col_character(),
  Comment = col_character(),
  Source = col_character()
)

LANGUAGES_COLSPEC = c(
  ID = col_character(),
  Name = col_character(),
  Macroarea = col_character(),
  Latitude = col_double(),
  Longitude = col_double(),
  Glottocode = col_character(),
  ISO639P3code = col_logical(),
  Coders = col_character(),
  provenance = col_character(),
  Family_name = col_character(),
  Family_level_ID = col_character(),
  Language_level_ID = col_character(),
  level = col_character(),
  lineage = col_character()
)

PARAMETERS_COLSPEC = c(
  ID = col_character(),
  Name = col_character(),
  Description = col_character(),
  patron = col_character(),
  name_in_french = col_character(),
  Grambank_ID_desc = col_character(),
  bound_morphology = col_character()
)

CODES_COLSPEC = c(
  ID = col_character(),
  Parameter_ID = col_character(),
  Name = col_character(),
  Description = col_character()
)

WIDE_COLSPEC = c(
  .default = col_integer(),
  Language_ID = col_character(),
  na_prop = col_double()
)

if (!dir.exists("output")) {
  dir.create("output")
}

if (!dir.exists("output/spatiophylogenetic_modelling/")) {
  dir.create("output/spatiophylogenetic_modelling/")
}
