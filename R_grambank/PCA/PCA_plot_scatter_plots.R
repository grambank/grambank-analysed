source("requirements.R")

#script written by Hedvig Skirgård and Simon Greenhill

OUTPUTDIR <- file.path('.',"output", 'PCA')

Language_meta_data <-  read_csv(GRAMBANK_LANGUAGES, col_types=LANGUAGES_COLSPEC) %>%
  dplyr::select(Language_ID = Language_level_ID, Family_name, Name, Macroarea) %>%
  distinct(Language_ID, .keep_all = T) %>%
  mutate(Family_name = ifelse(is.na(Family_name), "Isolate", Family_name))

#reading in the dataframe with PCA values for each language
GB_PCA_df <- suppressMessages(read_tsv(file.path(OUTPUTDIR, 'PCA_language_values.tsv'))) %>%
  dplyr::select(Language_ID, PC1, PC2, PC3) %>%
  left_join(Language_meta_data, by = "Language_ID" ) %>%
  dplyr::select(Language_ID, everything())

top_families <- GB_PCA_df %>%
  group_by(Family_name) %>%
  dplyr::summarise(n = n()) %>%
  arrange(-n) %>%
  .[1:15,1] %>%
  arrange(Family_name)

tided_PCA_descs<- suppressMessages(read_tsv(
    file.path(OUTPUTDIR, 'PCA_rotations.tsv')
))


# function to plot the PCA locations of a given language `family``
# highlighting the selected family using `color`.
plotPCAFamily <- function(df, family, color="orange", key="Family_name", label=FALSE) {
    df$InGroup <- df[[key]] == family
    p <- ggplot(df, aes(x = PC1, y = PC2), alpha = 1)
    p <- p + geom_jitter(color="snow3", size = 0.5, alpha = 0.4, shape = 16)
    p <- p + geom_jitter(
      data = filter(df, df$InGroup),
      aes(x = PC1, y = PC2), color=color,
      size = 0.5, alpha = 1
    )
    p <- p + geom_density2d(
      data = filter(df, df$InGroup),
      aes(x = PC1, y = PC2), color = color, alpha=0.4
    )

    if (label == TRUE) {
        p <- p + geom_text_repel(
            data = filter(df, df$InGroup),
            aes(x = PC1, y = PC2, label = Name),
            size=0.9,
            max.overlaps=20)
    }

    p <- p + coord_fixed()
    p <- p + theme_classic()
    p <- p + theme(
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank()
    )
    p <- p + ggtitle(family)
    p
}


#Mayan, Turkic, Ta-Ne-Omotic Algic, Chibchan
GB_PCA_Families <- list(
  plotPCAFamily(GB_PCA_df, "Afro-Asiatic", "#FF5733", label=TRUE),
  plotPCAFamily(GB_PCA_df, "Atlantic-Congo", "steelblue2", label=TRUE),
  plotPCAFamily(GB_PCA_df, "Austroasiatic", "midnightblue", label=TRUE),
  plotPCAFamily(GB_PCA_df, "Austronesian", "turquoise3", label=TRUE),
  plotPCAFamily(GB_PCA_df, "Central Sudanic", "purple4", label=TRUE),
  plotPCAFamily(GB_PCA_df, "Chibchan", "#E08B00", label=TRUE),
  plotPCAFamily(GB_PCA_df, "Dravidian", "blue", label=TRUE),
  plotPCAFamily(GB_PCA_df, "Indo-European", "seagreen", label=TRUE),
  plotPCAFamily(GB_PCA_df, "Mayan", "#00ABFD", label=TRUE),
  plotPCAFamily(GB_PCA_df, "Nakh-Daghestanian", "orange", label=TRUE),
  plotPCAFamily(GB_PCA_df, "Nuclear Trans New Guinea", "springgreen", label=TRUE),
  plotPCAFamily(GB_PCA_df, "Pama-Nyungan", "darkgreen", label=TRUE),
  plotPCAFamily(GB_PCA_df, "Sino-Tibetan", "violetred3", label=TRUE),
  plotPCAFamily(GB_PCA_df, "Turkic", "#FF689F", label=TRUE),
  plotPCAFamily(GB_PCA_df, "Uralic", "#DC71FA", label=TRUE),
  plotPCAFamily(GB_PCA_df, "Uto-Aztecan", "red1", label=TRUE)
)

# save each one
for (p in GB_PCA_Families) {
    filename <- file.path(OUTPUTDIR, sprintf("PCA_family_%s.pdf",  p$labels$title))
    cat(paste("writing", filename, "\n"))
    ggsave(filename, p, width = 6, height = 6, dpi = 600)

    # warnings here are ggrepel not labelling:
    # ggrepel: xx unlabeled data points (too many overlaps). Consider increasing max.overlaps

}


GB_PCA_multiGrid_family <- ggarrange(
  plotPCAFamily(GB_PCA_df, "Afro-Asiatic", "#FF5733"),
  plotPCAFamily(GB_PCA_df, "Atlantic-Congo", "steelblue2"),
  plotPCAFamily(GB_PCA_df, "Austroasiatic", "midnightblue"),
  plotPCAFamily(GB_PCA_df, "Austronesian", "turquoise3"),
  plotPCAFamily(GB_PCA_df, "Central Sudanic", "purple4"),
  plotPCAFamily(GB_PCA_df, "Chibchan", "#E08B00"),
  plotPCAFamily(GB_PCA_df, "Dravidian", "blue"),
  plotPCAFamily(GB_PCA_df, "Indo-European", "seagreen"),
  plotPCAFamily(GB_PCA_df, "Mayan", "#00ABFD"),
  plotPCAFamily(GB_PCA_df, "Nakh-Daghestanian", "orange"),
  plotPCAFamily(GB_PCA_df, "Nuclear Trans New Guinea", "springgreen"),
  plotPCAFamily(GB_PCA_df, "Pama-Nyungan", "darkgreen"),
  plotPCAFamily(GB_PCA_df, "Sino-Tibetan", "violetred3"),
  plotPCAFamily(GB_PCA_df, "Turkic", "#FF689F"),
  plotPCAFamily(GB_PCA_df, "Uralic", "#DC71FA"),
  ncol = 3, nrow = 5)


GB_PCA_multiGrid_family <- GB_PCA_multiGrid_family + rremove("x.text")
GB_PCA_multiGrid_family <- annotate_figure(
    GB_PCA_multiGrid_family,
    bottom = text_grob("PC1", size = 18, vjust = 0.4),
    left = text_grob("PC2", size = 18, rot = 90, hjust = -1, vjust = 1)
)

filename <- file.path(OUTPUTDIR, "PCA_multigrid_major_families.tiff")
filename_600 <- file.path(OUTPUTDIR, "PCA_multigrid_major_families_600_dpi.tiff")
filename_png <- file.path(OUTPUTDIR, "PCA_multigrid_major_families.png")
cat(paste("writing", filename, "\n"))
ggsave(filename, GB_PCA_multiGrid_family, width = 8, height = 10, dpi = 400)
ggsave(filename_600, GB_PCA_multiGrid_family, width = 8, height = 10, dpi = 600)
ggsave(filename_png, GB_PCA_multiGrid_family, width = 8, height = 10, dpi = 400)


GB_PCA_all_families <- ggplot(GB_PCA_df, aes(x = PC1, y = PC2, color = Family_name), size = 2, alpha = 0.6, na.rm = FALSE)
GB_PCA_all_families <- GB_PCA_all_families + geom_jitter() + coord_fixed()
GB_PCA_all_families <- GB_PCA_all_families + theme_classic() + theme(
    legend.position = "none",
    axis.title = element_text(size=10),
    plot.title = element_text(size=10),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
)
GB_PCA_all_families <- GB_PCA_all_families + ggtitle("All families")

filename <- file.path(OUTPUTDIR, "PCA_all_families.tiff")
cat(paste("writing", filename, "\n"))
ggsave(filename, GB_PCA_all_families, width = 7, height = 7, dpi = 600)



###Macroarea plots
GB_PCA_multiGrid_macroareas <- ggarrange(
    plotPCAFamily(GB_PCA_df, "Eurasia", "#00BA38", key="Macroarea"),
    plotPCAFamily(GB_PCA_df, "Papunesia", "#619CFF", key="Macroarea"),
    plotPCAFamily(GB_PCA_df, "Australia", "#B79F00", key="Macroarea"),
    plotPCAFamily(GB_PCA_df, "Africa", "#F8766D", key="Macroarea"),
    plotPCAFamily(GB_PCA_df, "South America", "#F564E3", key="Macroarea"),
    plotPCAFamily(GB_PCA_df, "North America", "#00BFC4", key="Macroarea"),
    ncol = 3, nrow = 2) + rremove("x.text")

GB_PCA_multiGrid_macroareas <- annotate_figure(GB_PCA_multiGrid_macroareas,
    bottom = text_grob("PC1", size = 18, vjust = 0),
    left = text_grob("PC2", size = 18, rot = 90, hjust = -2)
)


filename <- file.path(OUTPUTDIR, "PCA_multigrid_macroareas.tiff")
cat(paste("writing", filename, "\n"))
ggsave(filename, GB_PCA_multiGrid_macroareas, width = 7, height = 7, dpi = 600)


GB_PCA_macroarea <- ggplot(GB_PCA_df,
    aes(x = PC1, y = PC2, color = Macroarea), size = 2, alpha = 0.5) +
    geom_jitter() +
    coord_fixed() +
    theme_classic() +
    ggtitle("Macroareas")

filename <- file.path(OUTPUTDIR, "PCA_macroarea.tiff")
cat(paste("writing", filename, "\n"))
ggsave(filename, GB_PCA_macroarea, width = 7, height = 7, dpi = 600)


##plots with specific languages highlighted

#redefine function to large plot sizes

plotPCAFamily_highlights <- function(df, family, color="orange", key="Family_name") {
  df$InGroup <- df[[key]] == family
  p <- ggplot(df, aes(x = PC1, y = PC2), alpha = 1)
  p <- p + geom_jitter(color="snow3", size = 1.1, alpha = 1, shape = 16)
  p <- p + geom_jitter(
    data = filter(df, df$InGroup),
    aes(x = PC1, y = PC2), color=color,
    size = 1.2, alpha = 1
  )
  p <- p + geom_density2d(
    data = filter(df, df$InGroup),
    aes(x = PC1, y = PC2), color = color, alpha=0.4
  )
  p <- p + coord_fixed()
  p <- p + theme_classic()
  p <- p + theme(
    legend.position = "none",
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank()
  )
  p <- p + ggtitle(family)
  p
}


#Indo-European
lgs_to_highlight<- tibble(Language_ID = c("stan1293",
                                           "tokp1240",
                                           "hait1244",
                                           "swed1254",
                                           "lati1261",
                                           "ukra1253",
                                           "mara1378",
                                          "dara1250",
                                           "kumz1235",
                                          "west2386",
                                          "kash1277",
                                          "russ1263",
                                          "stan1290",
                                          "iris1253",
                                          "sana1297"))
df_IE_highlights <- GB_PCA_df %>%
  inner_join(lgs_to_highlight, by = "Language_ID")

plotPCAFamily_highlights(GB_PCA_df, "Indo-European", "seagreen") +
    geom_point(data = df_IE_highlights, aes(x = PC1, y = PC2), color ="seagreen") +
  geom_label_repel(data = df_IE_highlights, mapping = aes(x = PC1, y = PC2, label = Name)) +
  theme(plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14,angle = 90))

ggsave(file.path(OUTPUTDIR, "PCA_family_Indo-European_highlights.tiff"), width = 7, height = 7, dpi = 600)
ggsave(file.path(OUTPUTDIR, "PCA_family_Indo-European_highlights.png"), width = 7, height = 7, dpi = 600)


#Austroasiatic
lgs_to_highlight <- tibble(Language_ID = c("oldk1249", "lave1248", "seme1247", "mund1320", "pnar1238", "sant1410", "bond1245"))

df_austroasiatic_highlights <- GB_PCA_df %>%
  inner_join(lgs_to_highlight, by = "Language_ID")

plotPCAFamily_highlights(GB_PCA_df, "Austroasiatic", "midnightblue")+
  geom_point(data = df_austroasiatic_highlights, aes(x = PC1, y = PC2), color = "midnightblue") +
  geom_label_repel(data = df_austroasiatic_highlights, mapping = aes(x = PC1, y = PC2, label = Name)) +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14,angle = 90))

ggsave(file.path(OUTPUTDIR, "PCA_family_Austroasiatic_highlights.tiff"), width = 7, height = 7, dpi = 600)
ggsave(file.path(OUTPUTDIR, "PCA_family_Austroasiatic_highlights.png"), width = 7, height = 7, dpi = 600)
