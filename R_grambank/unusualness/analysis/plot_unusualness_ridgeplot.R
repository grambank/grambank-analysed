OUTPUTDIR <- file.path("output/unusualness/")
if (!dir.exists(OUTPUTDIR)) { dir.create(OUTPUTDIR) }		

OUTPUTDIR_tables <- file.path("output/unusualness/tables")
if (!dir.exists(OUTPUTDIR_tables)) { dir.create(OUTPUTDIR_tables) }		

OUTPUTDIR_plots <- file.path("output/unusualness/plots")
if (!dir.exists(OUTPUTDIR_plots)) { dir.create(OUTPUTDIR_plots) }		

# Load pkgs
source("requirements.R")

#reading in data
surprisal_fn <- paste0(OUTPUTDIR_tables, "/surprisal.tsv")
if(!file.exists(surprisal_fn)){
  source("unusualness/analysis/get_unusualness_bayesLCA.R")
}
gb <- read_tsv(file = surprisal_fn, show_col_types = F) %>% 
  dplyr::select(Language_ID, Surprisal,Estimator,Language_ID, Probability, Name, aes, AUTOTYP_area, Macroarea, Family_ID) %>% 
  filter(Estimator == "Kernel 15") %>% 
  dplyr::select(Language_ID, Surprisal, Probability, AUTOTYP_area)

median_df <- gb %>% 
  group_by(AUTOTYP_area) %>% 
  mutate(median_surprisal = median(Surprisal))

gb <- gb %>% 
  left_join(median_df)

gb$AUTOTYP_area <- fct_reorder(gb$AUTOTYP_area, gb$median_surprisal)

gb %>% 
   ggplot() +
  geom_density_ridges(aes(x = Surprisal, y = AUTOTYP_area, fill = median_surprisal), bandwidth = 0.131,
                      quantile_lines = T, quantile_fun = median, jittered_points = TRUE, point_size = 2, point_shape = 21  ,  position = position_points_jitter(height = 0))  +
  geom_label(data = median_df , aes(x = median_surprisal, y = AUTOTYP_area,
                                     label = round(median_surprisal, 2)), size = 3, nudge_x = 0, nudge_y = 0.4, alpha = 0.7, label.padding = unit(0.1, "lines")) +
  theme_classic() +
  theme(axis.title = element_blank(), 
        legend.position = "None", 
        axis.text = element_text(size = 12)) +
  scale_fill_viridis(option="magma", discrete = F, begin = 0.3) 

ggsave(filename = paste0(OUTPUTDIR_plots, "/unusualness_ridgeplot.png"), height = 8, width = 6)
  
