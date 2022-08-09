source("fun_def_h_load.R")
h_load(c("tidyverse"))

#script written by Hedvig Skirgård

cat("Making bar plots of the coverage of Grambank per macroarea and largest families.\n")

#getting table from Glottolog which contains information on the med level of the languoids
glottolog_cldf_df <- read_tsv("output/non_GB_datasets/glottolog-cldf_wide_df.tsv",col_types = cols()) %>% 
  mutate(Language_level_ID = ifelse(is.na(Language_level_ID), Language_ID, Language_level_ID)) %>%  #making language-level entities their own parent, so that we can use this column for aggregation easier.
  mutate(Family_ID = if_else(is.na(Family_ID) & level != "family", "Isolates", Family_ID)) %>%  #Lump isolates for easier viz %>% 
  mutate(med_summarised = if_else(str_detect(med, "grammar"), "grammar exists", "no grammar exists/unknown")) %>% 
  mutate(med_summarised = if_else(is.na(med_summarised), "no grammar exists/unknown", med_summarised))

family_names_df <- glottolog_cldf_df %>% 
  filter(is.na(classification) & level == "family") %>% 
  dplyr::select(Family_ID = Language_ID, Family_name = Name) %>% 
  full_join(tibble( Family_ID = c("Isolates"), Family_name = c("Isolates")), by = c("Family_ID", "Family_name"))

##grambank import and aggregation to language-level (i.e. merge dialects)
GB_wide <-read_tsv(file.path("output", "GB_wide", "GB_wide_strict.tsv"), col_types = cols()) %>% 
  dplyr::select(Language_ID, na_prop) 

#per family
df_for_family_plot <-  glottolog_cldf_df %>% 
  filter(level == "language") %>% 
  filter(Family_ID != 'book1242') %>% #removing bookkeeping languages
  filter(Family_ID != 'unat1236') %>% #removing unattested
  filter(Family_ID != 'arti1236') %>% #removing artificial
  filter(Family_ID != 'spee1234')  %>% #removing speech register 
  filter(Family_ID != 'uncl1493')  %>% #removing unclassifiable
  filter(Family_ID != 'pidg1258')  %>% #removing pidgins
  full_join(GB_wide, by = "Language_ID") %>% 
  filter(!is.na(Family_ID)) %>% 
  mutate(plot_value = if_else(!is.na(na_prop), "in Grambank", med_summarised)) %>%
  distinct(Language_level_ID, Family_ID, plot_value) 

families_sum <- df_for_family_plot %>% 
  group_by(Family_ID, plot_value) %>% 
  dplyr::summarise(n = n(), .groups = "drop_last") %>% 
  mutate(sum = sum(n, na.rm = T)) %>% 
  ungroup() %>% 
  complete(Family_ID, plot_value, fill= list(n = 0)) %>% 
  distinct(Family_ID, plot_value, n, sum) %>% 
  arrange(desc(sum)) %>% 
  .[1:45,] %>% 
  left_join(family_names_df, by = "Family_ID") %>% 
  mutate(position_text = ifelse(plot_value == "in Grambank" & Family_name == "Sign Language", n-3, n))

families_sum$Family_name <- fct_reorder(families_sum $Family_name, families_sum $sum)
#families_sum$plot_value <- factor(families_sum$plot_value, levels = c("no grammar exists/unknown", "grammar exists","in Grambank"))

col_vector <- c("orange", "purple4", "turquoise3")

family_bar_plot <- families_sum %>% 
ggplot(aes(x = Family_name, fill = plot_value, y = n)) +
  geom_bar(position= position_fill(), stat = "identity") + 
  theme_classic() +
  theme(legend.title=element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y =element_blank(),
        text = element_text(size=25),
        axis.text.x = element_text(angle = 50, hjust=1)) +
  scale_fill_manual(values = col_vector) +
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(x = `Family_name`, y = n, label = `n`, hjust = 0.5), size=6.5, colour = "white", position= position_fill(0.5)) 

tiff("output/coverage_plots/coverage_top_fifteen_families.tiff", width= 15, height =  8.46, units = "in",res = 300)
plot(family_bar_plot )
x <- dev.off()

png("output/coverage_plots/coverage_top_fifteen_families.png", width= 15, height =  8.46, units = "in", res = 300)
plot(family_bar_plot )
x <- dev.off()

#making table
df_for_family_plot %>% 
  left_join(family_names_df, by = "Family_ID") %>% 
  group_by(Family_name, plot_value) %>% 
  dplyr::summarise(n = n(), .groups = "drop_last") %>%
  ungroup() %>%   
  complete(Family_name, plot_value, fill = list(n = 0)) %>% 
  reshape2::dcast(Family_name~plot_value, value.var = "n" ) %>% 
  arrange(desc(`in Grambank`)) %>%
  write_tsv("output/coverage_plots/coverage_table.tsv")
  
#per continent

df_for_macroarea_plot <-  glottolog_cldf_df %>% 
  filter(level == "language") %>% 
  filter(Family_ID != 'book1242') %>% #removing bookkeeping languages
  filter(Family_ID != 'unat1236') %>% #removing unattested
  filter(Family_ID != 'arti1236') %>% #removing artificial
  filter(Family_ID != 'spee1234')  %>% #removing speech register 
  left_join(GB_wide, by = "Language_ID") %>% 
  mutate(plot_value = if_else(!is.na(na_prop), "in Grambank", med_summarised)) %>%
  distinct(Language_level_ID, Macroarea, plot_value) %>% 
  group_by(Macroarea, plot_value) %>% 
  dplyr::summarise(n = n(), .groups = "drop_last") %>% 
  mutate(sum = sum(n, na.rm = T)) %>% 
  ungroup() %>% 
  complete(Macroarea, plot_value, fill= list(n = 0)) %>% 
  distinct(Macroarea, plot_value, n, sum) %>% 
  arrange(desc(sum))

macroarea_plot <- df_for_macroarea_plot%>% 
  ggplot(aes(x = Macroarea, fill = plot_value, y = n)) +
  geom_bar(position= position_fill(), stat = "identity") + 
  theme_classic() +
  theme(legend.title=element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y =element_blank(),
        text = element_text(size=25),
        axis.text.x = element_text(angle = 50, hjust=1)) +
  scale_fill_manual(values = col_vector) +
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(x = `Macroarea`, y = n, label = `n`, hjust = 0.5), size=8, colour = "white", position= position_fill(0.5)) 


tiff("output/coverage_plots/coverage_macroarea.tiff", width= 10, height =  8.46, units = "in", res = 300)
plot(macroarea_plot)
x <- dev.off()

png("output/coverage_plots/coverage_macroarea.png", width= 10, height =  8.46, units = "in",res = 300)
plot(macroarea_plot)
x <- dev.off()

cat("Coverage bar plots made.\n")
