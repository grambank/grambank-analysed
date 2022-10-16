#Site-species matrix with sites as rows and species as columns if not provided,
# This can be either a matrix, a data.frame, or a Matrix::Matrix() object.
make_group_matrix <- function(df, grouping_var, remove=c(), threshold=3) {
  df <- as.data.frame(df)[c('Language_ID', grouping_var)]
  groups <- names(table(df[grouping_var])[table(df[grouping_var]) >= threshold])
  
  mtx <- matrix(0, nrow=length(groups), ncol=length(df$Language_ID), dimnames=list(groups, df$Language_ID))
  
  for (g in groups) {
    mtx[g, df[df[grouping_var] == g, 'Language_ID']] <- 1
  }
  
  if (length(remove) > 0) mtx[, remove] <- 0   # zero the languages in remove
  
  mtx
}
