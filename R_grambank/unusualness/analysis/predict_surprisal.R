# Load pkgs
source("requirements.R")

# Set working directory for output
#setup outpur dirs
OUTPUTDIR <- file.path("output/unusualness/")
if (!dir.exists(OUTPUTDIR)) { dir.create(OUTPUTDIR) }		

OUTPUTDIR_tables <- file.path("output/unusualness/tables")
if (!dir.exists(OUTPUTDIR_tables)) { dir.create(OUTPUTDIR_tables) }		

OUTPUTDIR_plots <- file.path("output/unusualness/plots")
if (!dir.exists(OUTPUTDIR_plots)) { dir.create(OUTPUTDIR_plots) }		

surprisal_fn <- paste0(OUTPUTDIR_tables, "/surprisal.tsv")
if(!file.exists(surprisal_fn)){
  source("unusualness/analysis/get_unusualness_bayesLCA.R")
}

#read in data
gb <- read.delim(file = surprisal_fn, sep = "\t") %>% 
  dplyr::select(Language_ID, aes, Surprisal, Estimator) %>% 
  filter(Estimator == "Kernel 30") %>% 
  inner_join(lgs_in_analysis, by = "Language_ID") %>% #subset to the lgs where we have phylo prec matrices
  mutate(Endangerement=ifelse(aes %in% c("threatened","moribund","nearly_extinct"),"endangered",aes)) # Recode endangerment

### NEXT PART REQUIRES MATRICES ETC

#########################################
## (5) Model unusualness in terms of genealogical, areal covariates, and endangerement status
#########################################


phylogenetic_tree <- phylogenetic_tree %>% keep.tip(inner_joined_df$Glottocode)

#brms
#making a covariance matrix of the tree
vcv_tree <- vcv.phylo(phylogenetic_tree)


###BRMS

formula_for_brms <- unusualness_score ~ L1_log10 + L2_log10 + Is_Written + Official +
  (1 | gr(Glottocode, cov = vcv_tree)) +
  (L1_log10 + L2_log10 + Is_Written + Official | Family_ID)

full_model <- brms::brm(formula = formula_for_brms,
                        data = filter(inner_joined_df, !is.na(L1), !is.na(L2)),
                        data2 = list(vcv_tree= vcv_tree),
                        iter = 7500,
                        iter = 10000,
                        cores = 4,
                        control = list(adapt_delta =0.99, max_treedepth=15)
) %>% add_criterion("waic")
full_model %>% broom.mixed::tidy() %>% write_csv("unusualness/analysis/full_model.csv")

simplified_model <- brms::brm(unusualness_score ~ 1 + (1 | gr(Glottocode, cov = vcv_tree)),
                              data = filter(inner_joined_df, !is.na(L1), !is.na(L2)),
                              data2 = list(vcv_tree= vcv_tree),
                              iter = 7500,
                              iter = 25000,
                              control = list(adapt_delta =0.99, max_treedepth=15)
) %>% add_criterion("waic")
simplified_model %>% broom.mixed::tidy() %>% write_csv("unusualness/analysis/simplified_model.csv")

loo_compare(full_model, simplified_model, criterion="waic") %>% as.tibble() %>% write_csv("unusualness/analysis/model_comparison.csv")
