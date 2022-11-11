source("requirements.R")

OUTPUTDIR <- "output/dists/"
if (!dir.exists(OUTPUTDIR)) {dir.create(OUTPUTDIR)}

cfx_macroarea_list <-   read_tsv("output/dists/cfx_Macroarea_cut_off_0_list.tsv", show_col_types = F)

cfx_macroarea_list$Vars <- fct_reorder(cfx_macroarea_list$Vars, cfx_macroarea_list$Value_cfx)

cfx_macroarea_list %>% 
  ggplot(aes(x = Vars, y = Value_cfx)) +
  geom_bar(aes(fill = 1 - Value_cfx), stat = "identity") +
  theme_classic() +
  theme(text = element_text(angle = 70, hjust = 1, size = 20), 
        legend.position = "None", 
        axis.title.x = element_blank()) +
  scale_fill_viridis()+
  ylab("Cultural Fixation Score")

mean(cfx_macroarea_list$Value_cfx)

ggsave(filename = file.path(OUTPUTDIR, "cfx_barplot_macroarea.png"), height =  7.89, width =  8.61)

fns <- list.files(path = OUTPUTDIR, pattern = ".*cfx_Family_ID_list_cut_off.*")

for(fn in fns){

cfx_Family_ID_list <- read_tsv(fn)

plot_title <- fn

cfx_Family_ID_list %>% 
  ggplot(aes(x = Vars, y = Value_cfx)) +
  geom_bar(aes(fill = 1 - Value_cfx), stat = "identity") +
  theme_classic() +
  theme(legend.position = "None", 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank()) +
  ggtitle(plot_title) 

mean(cfx_Family_ID_list$Value_cfx)

ggsave(filename = paste0(OUTPUTDIR, "cfx_barplot", fn, ".png"), height =  7.89, width =  8.61)

}
