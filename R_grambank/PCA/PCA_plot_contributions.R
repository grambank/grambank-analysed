#PCA analysis over the imputed data. Here we use the imputed data from random forests to run PCA.
source("requirements.R")

#script written by Hedvig Skirgård and Simon Greenhill

OUTPUTDIR <- file.path('.', "output", 'PCA')

# Part of this script is taken from Julia Silge's PCA visualisations of
# stack overflow data especially the tided contributions graphs:
#    https://juliasilge.com/blog/stack-overflow-pca/
tidied_PCA_descs <- read.delim(file.path("output", 'PCA', 'PCA_rotations.tsv'), sep = "\t")

PCA_prop_variance_df<- tidied_PCA_descs %>%
  dplyr::distinct(PC, `Proportion.of.Variance`)

# set up our plotting theme - starts with theme_classic and then modify some parts
theme_grambank_pca <- function(base_size = 12, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.ticks.x = element_blank()
    )
}

for (pcX in c("PC1", "PC2", "PC3", "PC4")) {
  prop_variance <- PCA_prop_variance_df[PCA_prop_variance_df$PC==pcX,2] %>%
      as.matrix() %>% as.vector()
  
  p <- tidied_PCA_descs %>%
    filter(PC == pcX) %>%
    top_n(40, abs(Contribution)) %>%
    mutate(Grambank_ID_desc = reorder(Grambank_ID_desc, Contribution)) %>%
    ggplot(aes(Grambank_ID_desc, Contribution, fill = Contribution)) +
    geom_col(show.legend = FALSE, alpha = 0.8) +
    theme_grambank_pca() +
    scale_fill_viridis(option="magma") +
    labs(title= paste(pcX),
      subtitle = paste0(" Proportion of Variance: ", prop_variance * 100, "%"),
      x = "Grambank features",
      y = sprintf("Relative importance in %s", pcX)
    ) +
    coord_flip()
  
  filename <- file.path(OUTPUTDIR, sprintf("PCA_contributions_%s.tiff", pcX))
  filename_png <- file.path(OUTPUTDIR, sprintf("PCA_contributions_%s.png", pcX))
  
    cat(paste("writing", filename, "\n"))
  ggsave(filename, p, width = 6, height = 8)
  ggsave(filename_png, p, width = 6, height = 8)
  
  }

PC1_plot <- tidied_PCA_descs %>%
  filter(PC == "PC1") %>%
  top_n(40, abs(Contribution)) %>%
  mutate(Grambank_ID_desc = reorder(Grambank_ID_desc, Contribution)) %>%
  ggplot(aes(Grambank_ID_desc, Contribution, fill = Contribution)) +
  geom_col(show.legend = FALSE, alpha = 0.8) +
  theme_grambank_pca() +
  scale_fill_viridis() +
  ylim(-0.25, 0.29) +
  ggtitle("a") +
  theme(title = element_text(face = "bold"), axis.title = element_blank()) +
  coord_flip()

PC2_plot <- tidied_PCA_descs %>%
  filter(PC == "PC2") %>%
  top_n(40, abs(Contribution)) %>%
  mutate(Grambank_ID_desc = reorder(Grambank_ID_desc, Contribution)) %>%
  ggplot(aes(Grambank_ID_desc, Contribution, fill = Contribution)) +
  geom_col(show.legend = FALSE, alpha = 0.8) +
  theme_grambank_pca() +
  scale_fill_viridis(option="magma") +
  ylim(-0.25, 0.29) +
  ggtitle("b") +
  theme(title  = element_text(face = "bold"), axis.title = element_blank()) +
  coord_flip()

GB_PCA_contributions_multiGrid <- ggarrange(
  plot(PC1_plot),
  plot(PC2_plot),
  ncol = 2, nrow = 1, align = "h")

GB_PCA_contributions_multiGrid <- annotate_figure(
  GB_PCA_contributions_multiGrid,
  bottom = text_grob("Relative importance", size = 14),
  left = text_grob("Grambank Features", size = 14, rot = 90, hjust = -0.05)
)

filename <- file.path(OUTPUTDIR, "PCA_contributions_PCA1_PCA2.tiff")
filename_png <- file.path(OUTPUTDIR, "PCA_contributions_PCA1_PCA2.png")
cat(paste("writing", filename, "\n"))
ggsave(filename, GB_PCA_contributions_multiGrid, width = 10, height = 8)
ggsave(filename_png, GB_PCA_contributions_multiGrid, width = 10, height = 8)