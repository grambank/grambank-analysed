#PCA analysis over the imputed data. Here we use the imputed data from random forests to run PCA.
library(readr)
library(tidyverse)
library(reshape2)
library(nFactors)

#script written by Hedvig Skirgård

OUTPUTDIR <- file.path("output/PCA")
if (!dir.exists(OUTPUTDIR)) { dir.create(OUTPUTDIR) }

cat("Running Principal Component Analysis on the cropped, dialect-merged, binarised and imputed dataset.\n")

GB_imputed <- read_tsv(file.path("output", "GB_wide", "GB_wide_imputed_binarized.tsv"), col_types= cols()) %>%
    column_to_rownames("Language_ID") %>%
    as.matrix()

#Shorter names for the Grambank parameters that can be used in graphs
Grambank_id_abbrev <- read_tsv(
    file.path("output", "GB_wide", "parameters_binary.tsv"),
    show_col_types=FALSE) %>% 
    dplyr::select(Parameter_ID = ID, Grambank_ID_desc)

# run PCA
GB_imputed_PCA <- prcomp(GB_imputed, scale. = TRUE)

#extracting the components dataframe from the prcomp object
GB_PCA_df <- as.data.frame(GB_imputed_PCA$x)

GB_PCA_df %>%
    rownames_to_column("Language_ID") %>%
    write_tsv(file.path(OUTPUTDIR, 'PCA_language_values.tsv'))


# tided dataframe based on the prcomp object, every row is a component
# with the propotion of variance information tided dataframe based on
# the prcomp object, every row is a component with the propotion of
# variance information
PCA_summary_importance <- t(summary(GB_imputed_PCA)$importance) %>%
    as.data.frame() %>%
    rownames_to_column("PC") %>%
    mutate(`Proportion of Variance` = round(`Proportion of Variance`, digits = 2)) %>%
    dplyr::select(`Proportion of Variance`, everything())


#Tidying the component contributions into a long format.
tidied_pca <- GB_imputed_PCA$rotation %>%
    as.data.frame() %>% 
    rownames_to_column("Parameter_ID") %>%
    dplyr::select(Parameter_ID, PC1, PC2, PC3, PC4) %>%
    reshape2::melt(id.vars = "Parameter_ID") %>% 
    dplyr::rename(PC = variable, Contribution = value) %>%
    left_join(PCA_summary_importance,  by = "PC") %>%
    left_join(Grambank_id_abbrev, by="Parameter_ID") %>%
    tidyr::separate(col = Grambank_ID_desc, into = c("Parameter_ID", "Grambank_ID_desc"), sep = " ") 

write_tsv(tidied_pca, file.path(OUTPUTDIR, 'PCA_rotations.tsv'))
cat("PCA tables written to dir PCA.\n")

#testing to evaluate the optimal number of components
ev <- eigen(cor(GB_imputed)) # get eigenEstimates
ap <- nFactors::parallel(
    subject=nrow(GB_imputed),
    var=ncol(GB_imputed),
    rep=100,
    cent=0.05)
nS <- nFactors::nScree(x=ev$values, aparallel=ap$eigen$qevpea)

optimal_components <- nS$Components$nparallel

plotnScree(nS)

png(filename = file.path(OUTPUTDIR, "PCA_nscree_plot.png"))
plotnScree(nS)
x <- dev.off()

prop_explained <- sum(summary(GB_imputed_PCA)$importance[2, 1:optimal_components])

nScree_summary_string <- paste(
    "The optimal number of components is", optimal_components, "and they explain", prop_explained * 100, "% of the variation.\n"
)

writeLines(
    nScree_summary_string,
    con = file.path(OUTPUTDIR, "PCA_nScree_summary.txt")
)

cat(nScree_summary_string)