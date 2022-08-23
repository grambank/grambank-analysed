
source("dists/fixation_scores_muthukrishna_GB_fun.R")

Language_meta_data <- Language_meta_data %>% 
  mutate(americas = ifelse(str_detect(Macroarea, "merica"), "America", "Not americas"))

fun_cfx(df = Language_meta_data, group = "Macroarea")
