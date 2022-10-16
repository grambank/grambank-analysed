source("functional_richness/requirements_fr.R")

gb_mfd <- read.delim('output/GB_wide/GB_wide_imputed_binarized.tsv', na.strings='NA') %>% 
  column_to_rownames("Language_ID")

# figure out trait categories

tr_cat <- data.frame(trait_name=colnames(gb_mfd), trait_type="N", trait_weight=1) # N for Nominal
# god this is annoying -- these need to be factors, and mFD complains that they're not
# ordinal given the binary data.
gb_mfd[gb_mfd == "0"] <- 'a'
gb_mfd[gb_mfd == "1"] <- 'p'

gb_mfd[] <- lapply(gb_mfd, function(x) as.factor(x))

gb.dist <- mFD::funct.dist(
    sp_tr         = gb_mfd,
    tr_cat        = tr_cat,
    metric        = "gower",
    scale_euclid  = "scale_center",
    ordinal_var   = "classic",
    weight_type   = "equal",
    stop_if_NA    = TRUE)

saveRDS(gb.dist, file = "output/functional_richness/gb.gower.mfd.RDS")
