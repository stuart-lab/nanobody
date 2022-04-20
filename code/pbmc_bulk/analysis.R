library(Signac)
library(ggplot2)

bigWigs <- list(
  "H3K27me3 (multi)" = "data/pbmc_bulk/mapped/plex/me3.bw",
  "H3K27ac (multi)" = "data/pbmc_bulk/mapped/plex/ac.bw",
  "H3K27me3 (mono)" = "data/pbmc_bulk/mapped/mono/me3.bw",
  "H3K27ac (mono)" = "data/pbmc_bulk/mapped/mono/ac.bw"
)

bw <- BigwigTrack(
  region = "chr1-156630171-156920171",
  bigwig = bigWigs,
  y_label = "Normalized coverage",
  ymax = "q90",
  bigwig.scale = "common",
  smooth = 1000
) + theme(legend.position = "none")

bw <- bw + scale_fill_manual(values = c("#D3145A", "#F98401",
                                        "#676767", "#676767"))

ggsave("plots/pbmc_bulk/bulk_covplot_pbmc.png", plot = bw, height = 5, width = 10, units = "in")
