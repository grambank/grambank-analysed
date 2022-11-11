source("requirements.R")
#script written by Luke Maurits

d <- read_tsv("output/spatiophylogenetic_modelling/featurewise/posteriors_df.tsv", show_col_types = F) %>%
	filter(str_detect(fn, "_kappa_2_sigma_1.15_pcprior0.1")) %>%
	group_by(Feature_ID) %>%
	summarise(spatial=mean(`Precision for spatial_id_in_dual`),
	          phylo=mean(`Precision for phylo_id_in_dual`))

domains <- read_tsv("output/GB_wide/table_theo_scores_supp.tsv") %>%
	dplyr::select(Feature_ID, Main_domain)

d <- left_join(d, domains)

#d %>%
#	group_by(Main_domain) %>%
#	summarise(N=length(Feature_ID))

ggplot(d) +
	geom_density(aes(x=phylo, fill=Main_domain, colour=Main_domain), alpha=0.25) +
  theme_classic()

null_prior <- prior(normal(0, 1.5), class="Intercept")
full_prior <- prior(normal(0, 1.5), class="b")

m_phylo_null <- brm(phylo ~ 1, prior=null_prior, data=d, family="beta")
m_phylo_null <- add_criterion(m_phylo_null, "waic")
m_phylo_full <- brm(phylo ~ 0 + Main_domain, prior=full_prior, data=d, family="beta")
m_phylo_full <- add_criterion(m_phylo_full, "waic")
loo_compare(m_phylo_null, m_phylo_full, criterion="waic")

m_spatial_null <- brm(spatial ~ 1, prior=null_prior, data=d, family="beta")
m_spatial_null <- add_criterion(m_spatial_null, "waic")
m_spatial_full <- brm(spatial ~ 0 + Main_domain, prior=full_prior, data=d, family="beta")
m_spatial_full <- add_criterion(m_spatial_full, "waic")
loo_compare(m_spatial_null, m_spatial_full, criterion="waic")

`Model` <- c("null model (spatial)",
                 "domain model (spatial)",
                 "null model (phylogenetic)",
                 "domain model (phylogenetic)")

WAIC <- c(waic(m_spatial_null)$estimates[3,1],
           waic(m_spatial_full)$estimates[3,1],  
           waic(m_phylo_null)$estimates[3,1],  
           waic(m_phylo_full)$estimates[3,1])

`SE (WAIC)` <- c(waic(m_spatial_null)$estimates[3,2],
              waic(m_spatial_full)$estimates[3,2],  
              waic(m_phylo_null)$estimates[3,2],  
              waic(m_phylo_full)$estimates[3,2])

table_for_sm <- cbind(Model, WAIC, `SE (WAIC)`)

table_for_sm %>% 
  as.data.frame() %>% 
  mutate(WAIC = as.numeric(WAIC) %>% round(2)) %>% 
  mutate(`SE (WAIC)` = as.numeric(`SE (WAIC)`) %>% round(2)) %>% 
write_tsv("output/spatiophylogenetic_modelling/predict_sp_table_sm.tsv", na = "")
                 