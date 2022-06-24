library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)

bmmc <- readRDS("objects/bmmc_dual.rds")

colormap <- list("H3K27me3" = "#D3145A", "H3K27ac" = "#F98401", "RNAPII" = "#036C9A")

p1 <- DimPlot(bmmc, reduction = "umap.me3", group.by = "celltype", pt.size = 0.1, label = TRUE, label.size = 3, repel = TRUE) +
  theme_void() +
  theme(text = element_text(size = 10), legend.position = "none") +
  ggtitle("H3K27me3")
p2 <- DimPlot(bmmc, reduction = "umap.ac", group.by = "celltype", pt.size = 0.1, label = TRUE, label.size = 3, repel = TRUE) +
  theme_void() +
  theme(text = element_text(size = 10), legend.position = "none") +
  ggtitle("H3K27ac") +
  NoLegend()
p3 <- DimPlot(bmmc, reduction = "umap.wnn", group.by = "celltype", label = TRUE, label.size = 4, repel = TRUE) + 
  ggtitle("BMMCs: H3K27me3 + H3K27ac") +
  theme_classic() +
  NoLegend() +
  ylab('UMAP 2') + xlab("UMAP 1")

pp <- ((p1 / p2) | p3) + plot_layout(widths = c(1, 2))
ggsave(filename = "plots/bmmc/dimplots_bmmc.png", plot = pp, height = 6, width = 8)

# fragment counts
fm <- CountFragments(
  fragments = Fragments(bmmc[["me3"]])[[1]]@path,
  cells = colnames(bmmc)
)

fa <- CountFragments(
  fragments = Fragments(bmmc[["ac"]])[[1]]@path,
  cells = colnames(bmmc)
)

mean(fm$frequency_count) # 1217
mean(fa$frequency_count) # 326
sd(fm$frequency_count)  # 1274
sd(fa$frequency_count)  # 334

# create plot
fm$mark <- "H3K27me3"
fa$mark <- "H3K27ac"
all_counts <- rbind(fm, fa)

p0 <- ggplot(all_counts, aes(mark, frequency_count, fill = mark)) +
  geom_violin(size = 0.2) +
  scale_y_log10() +
  theme_bw() + 
  ylab("Fragments") +
  xlab("Mark") +
  ggtitle("Total fragments") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c(colormap$H3K27ac, colormap$H3K27me3))

ggsave(filename = "plots/bmmc/fragments_bmmc.png", plot = p0, height = 3, width = 3)
