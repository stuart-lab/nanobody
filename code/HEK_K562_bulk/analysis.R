library(Signac)
library(ggplot2)

bigWigs <- list(
  "RNAPII (multi)" = "data/HEK_K562_bulk/mapped/K562-plex-Pol2.bw",
  "H3K27me3 (multi)" = "data/HEK_K562_bulk/mapped/K562-plex-K27me.bw",
  "H3K27ac (multi)" = "data/HEK_K562_bulk/mapped/K562-plex-K27ac.bw",
  "RNAPII (mono)" = "data/HEK_K562_bulk/mapped/K562-mono-Pol2.bw",
  "H3K27me3 (mono)" = "data/HEK_K562_bulk/mapped/K562-mono-K27me.bw",
  "H3K27ac (mono)" = "data/HEK_K562_bulk/mapped/K562-mono-K27ac.bw"
)

bw <- BigwigTrack(
  region = "chr11-69900000-70910000",
  bigwig = bigWigs,
  y_label = "Normalized coverage",
  ymax = "q90",
  bigwig.scale = "common",
  smooth = 1000
) + theme(legend.position = "none")

bw <- bw + scale_fill_manual(values = c("#036C9A", "#D3145A", "#F98401",
                                        "#676767", "#676767", "#676767"))

ggsave(
  filename = "plots/hek_k562_bulk/bulk_covplot_k562.png",
  plot = bw,
  height = 5,
  width = 10,
  units = "in"
)
