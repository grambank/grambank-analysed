# Simulate data 
set.seed(7573673)

source("fun_def_h_load.R")
h_load(pkg = c("ape", "dplyr"))

#This script takes a while to run. If you set the variable beep to 1 it will make a little pling sound when it's finished, which can be handy. By default it's set to 0. Do not set it to 1 if you are running the script on a cluster, as clusters don't have audio capabilities.
beep <- 1
if(beep == 1){
  h_load("beepr")
  }

# Simulation parameters
proportions = c(0.1, 0.25, 0.4)
# We want to get the proportions in 'proportions'
# But we are happy with proportions that are within the allowable variation
allowable_variation = 0.05 
lambda_transformations = seq(0, 1, by = 0.3)
lambda_transformations[1] = 0.01
iterations = 15
try_times = 300

# save simulations here
simulated_location = "output/spatiophylogenetic_modelling/simulated_data/"
if(!dir.exists(simulated_location)) {dir.create(simulated_location, recursive = TRUE)}

# Phylogeny used
tree_fn <- "output/spatiophylogenetic_modelling/processed_data/EDGE_pruned_tree.tree"
if(!(file.exists(tree_fn))){
  source("spatiophylogenetic_modelling/processing/pruning_EDGE_tree.R")}
tree = read.tree(tree_fn)

tree$edge.length = tree$edge.length / 1000

tree_tips_df <- tree$tip.label %>% 
  as.data.frame() %>% 
  rename("Language_ID"= ".")

# locations
locations_df = read.delim(
  file = "output/non_GB_datasets/glottolog-cldf_wide_df.tsv", sep = "\t") %>% 
  dplyr::select(Language_ID, Latitude, Longitude) %>% 
  inner_join(tree_tips_df, by = "Language_ID") #subset to matches in tree

#df to bind to in the for-loop
df <- locations_df %>% 
  dplyr::select(Language_ID)

for(iter in 1:iterations){
  cat("Iteration", iter, "of", max(iterations), "...\n")
  for(i in seq_along(proportions)){
    desired_proportion = proportions[i]
    for(j in seq_along(lambda_transformations)){
      searching = TRUE
      desired_lambda = lambda_transformations[j]
      try = 0
      cat("Searching for a proportion of", desired_proportion, "...\n")
      cat("\twith a lambda of", desired_lambda, "...\n")
      while(searching){
        if(try == try_times){ searching = FALSE }
        q = matrix(c(round(runif(n = 2, min = 0, max = 1), 1)),
                   nrow = 2,
                   ncol = 2) * c(-1, 1, 1, -1)
        y = 
          geiger::sim.char(
            geiger::rescale(tree,
                            desired_lambda,
                            model = "lambda"), 
            q, 
            model="discrete")[,1,]  
        
        observed_proportion = min(prop.table(table(y)))
        if(observed_proportion > desired_proportion - allowable_variation &
           observed_proportion < desired_proportion + allowable_variation)
          searching = FALSE
        
        cat(paste0("I'm on try ", try, " out of max ", try_times,".\n"))
        try = try + 1
        
      }
      
      if(length(table(y)) == 1){
        # In some cases the simulated data is all 1s or all 0s and 
        # we don't want that - so just skip this
        next 
      } else {
        cat(paste0("... stopped search.\n"))
        out_df = data.frame(y = y - 1, # make data 0 & 1 
                            Language_ID = tree$tip.label)
        
        value_col_name <- paste0("Prop", desired_proportion, "_Lambda",desired_lambda, "Iter", iter)
        colnames(out_df) <- c(value_col_name, "Language_ID")
        
        df = left_join(out_df, df, by = "Language_ID") %>% 
          dplyr::select(Language_ID, everything())
        
         
          write_delim(x = df, file = paste(simulated_location, "simulated_data_df.tsv"),delim = "\t")  
      }
        
     
    }
  }
}

if("y" %in% colnames(df)){
df %>% 
  dplyr::select(-y) %>% 
  write_delim(file = paste(simulated_location, "simulated_data_df.tsv"),delim = "\t")  
}

if(beep == 1){
beepr::beep()
}