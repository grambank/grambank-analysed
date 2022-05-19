source("requirements.R")

library(NCmisc, devtools)

fns <- list.files(path = ".", pattern = "*.R$", recursive = T)

df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(df) <- c("package", "functions", "filename")
df$package <- as.character(df$package)
df$functions <- as.list(df$functions)
df$filename <- as.character(df$filename)

index <- 0

for(fn in fns){
  #fn <- fns[2]

  index <- index +1

    cat(paste0("Wrangling  ", fn, ". It's number ", index, " out of ",   length(fns), ".\n"))
  df <-  NCmisc::list.functions.in.file(fn) %>% as.data.frame() %>%
    rownames_to_column("package") %>% 
    rename(functions = ".") %>%
    mutate(filename = fn) %>% 
    full_join(df, by = c("package", "functions", "filename"))
  }

df_used_functions <- df %>% 
  unnest("functions") %>% 
  mutate(package = str_replace_all(package, "\"", ""))  %>% 
  mutate(package = str_replace_all(package, "c\\(", "")) %>% 
  mutate(package = str_replace_all(package, "\\)", "")) %>% 
  mutate(package = str_split(package, ",")) %>%   
  unnest(package) %>% 
  mutate(package = str_replace_all(package, "package\\:", "")) %>%
  mutate(package = trimws(package)) %>% 
  mutate(functions = trimws(functions)) %>% 
  distinct()  %>% 
  mutate(used = "yes")

loaded_packages <- devtools::loaded_packages() %>% 
  as.data.frame() %>%
  dplyr::select(package) %>% 
  mutate(loaded = "yes")

df_used_functions %>% 
  group_by(filename) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) %>% 
  .[1:10,]

full_join(df_used_functions, loaded_packages) %>% 
  filter(loaded == "yes") %>% 
  filter(is.na(used)) %>% 
  write_tsv("unsused_pkgs.tsv")


