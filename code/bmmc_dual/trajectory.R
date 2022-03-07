library(Signac)
library(Seurat)
library(ggplot2)
library(Matrix)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
library(RColorBrewer)
library(RcppRoll)

if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
  remotes::install_github('satijalab/seurat-wrappers')
}

library(SeuratWrappers)

if (!requireNamespace("monocle3", quietly = TRUE)) {
  setRepositories(ind=1:2)
  remotes::install_github('cole-trapnell-lab/leidenbase')
  remotes::install_github('cole-trapnell-lab/monocle3')
}

library(monocle3)
library(future)
plan("multicore", workers = 8)
options(future.globals.maxSize = Inf)


bmmc <- readRDS("objects/bmmc_dual.rds")
Idents(bmmc) <- "celltype"
cells.use <- c("HSPC", "GMP/CLP", "Pre-B", "B", "Plasma")

bcell <- bmmc[, bmmc$celltype %in% cells.use]

DefaultAssay(bcell) <- "me3"
bcell <- FindTopFeatures(bcell)
bcell <- RunTFIDF(bcell)
bcell <- RunSVD(bcell, reduction.name = "lsi.me3")
bcell <- RunUMAP(bcell, reduction = "lsi.me3", reduction.name = "umap.me3", dims = 2:20)

## run monocle3
bcell@reductions$umap <- bcell@reductions$umap.me3 # needs to be named UMAP for monocle
bcell.cds <- as.cell_data_set(bcell)
bcell.cds <- cluster_cells(cds = bcell.cds, reduction_method = "UMAP")
plot_cells(bcell.cds)
bcell.cds <- learn_graph(bcell.cds, use_partition = FALSE)
root_cells <- WhichCells(bcell, idents = "HSPC")
bcell.cds <- order_cells(bcell.cds, reduction_method = "UMAP", root_cells = tail(root_cells, 10))

bcell <- AddMetaData(
  object = bcell,
  metadata = bcell.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "B.cell.pseudotime"
)

fp1 <- FeaturePlot(bcell, "B.cell.pseudotime", reduction = "umap.me3") +
  scale_color_viridis_c(option = "C") +
  theme_void() +
  ggtitle("") +
  theme(legend.position = "none")

ggsave("plots/bmmc/trajectory.png", plot = fp1, height = 4, width = 4)

## find features correlated with pseudotime
chr.use <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[1:22]

bins_me3 <- GenomeBinMatrix(
  fragments = Fragments(bcell[['me3']]),
  genome = chr.use,
  binsize = 10000,
  cells = colnames(bcell)
)

bins_ac <- GenomeBinMatrix(
  fragments = Fragments(bcell[['ac']]),
  genome = chr.use,
  binsize = 10000,
  cells = colnames(bcell)
)

DefaultAssay(bcell) <- "me3"
bcell <- FindNeighbors(bcell, reduction = "lsi.me3", dims = 2:20, k = 50)
knn <- bcell[['me3_nn']]

knn_smooth <- function(counts, knn) {
  knn_norm <- knn / rowSums(knn)
  smoothed <- Matrix::tcrossprod(counts, knn_norm)
  return(smoothed)
}

smoothed_me3 <- knn_smooth(counts = bins_me3, knn = knn)
smoothed_ac <- knn_smooth(counts = bins_ac, knn = knn)

# filter regions based on total counts
keep_region <- (rowSums(smoothed_ac) > 5) & (rowSums(smoothed_me3) > 5)

smoothed_me3 <- smoothed_me3[keep_region, ]
smoothed_ac <- smoothed_ac[keep_region, ]

# depth normalize
smoothed_me3 <- t(t(smoothed_me3) / colSums(smoothed_me3))
smoothed_ac <- t(t(smoothed_ac) / colSums(smoothed_ac))

## correlation
pstime <- bcell$B.cell.pseudotime[colnames(bins_me3)]
peak_cor_me3 <- cor(t(as.matrix(smoothed_me3)), pstime)[, 1]
peak_cor_ac <- cor(t(as.matrix(smoothed_ac)), pstime)[, 1]

cor_df <- data.frame(me3 = peak_cor_me3, ac = peak_cor_ac)
cor_df$diff <- cor_df$me3 - cor_df$ac

top_cor <- cor_df[abs(cor_df$diff) > 0.5, ]
tc <- top_cor

# order cells by pseudotime
porder <- order(pstime)
ac_peaks <- smoothed_ac[rownames(tc), porder]
me3_peaks <- smoothed_me3[rownames(tc), porder]

# order based on which point in trajectory has max value using broader smoothing
rowsmoothed <- t(roll_sum(as.matrix(t(me3_peaks)), n = 200))
rownames(rowsmoothed) <- rownames(me3_peaks)
rowmax <- apply(rowsmoothed, 1, which.max)
me3_peaks <- me3_peaks[order(rowmax), ]
ac_peaks <- ac_peaks[order(rowmax), ]

rowsmoothed_ac <- t(roll_sum(as.matrix(t(ac_peaks)), n = 100))
mm_ac <- t(apply(as.matrix(rowsmoothed_ac), 2, rev))

rowsmoothed_me3 <- t(roll_sum(as.matrix(t(me3_peaks)), n = 100))
mm_me3 <- t(apply(as.matrix(rowsmoothed_me3), 2, rev))

ac_max <- quantile(mm_ac, probs = seq(0, 1, 0.05))
me3_max <- quantile(mm_me3, probs = seq(0, 1, 0.05))

me3_plot <- mm_me3
me3_plot[me3_plot > me3_max[['95%']]] <- me3_max[['95%']]
me3_plot[me3_plot < me3_max[['5%']]] <- me3_max[['5%']]

ac_plot <- mm_ac
ac_plot[ac_plot > ac_max[['95%']]] <- ac_max[['95%']]
ac_plot[ac_plot < ac_max[['5%']]] <- ac_max[['5%']]

cols.use <- rev(colorRampPalette(brewer.pal(9, "YlGnBu"))(100))
png("./plots/hmap_ac.png", width = 5000, height = 5000, units = "px")
image(ac_plot, col = cols.use, main = "",
      xaxt = 'n', yaxt = 'n', ylab = "", xlab = "")
dev.off()

png("./plots/hmap_me3.png", width = 5000, height = 5000, units = "px")
image(me3_plot, col = cols.use, main = "",
      xaxt = 'n', yaxt = 'n', ylab = "", xlab = "")
dev.off()

# group peaks by whether they gain me3 or gain ac
tc <- tc[order(tc$me3, decreasing = TRUE), ]
peaks_me3 <- rownames(tc[tc$me3 > 0, ])
peaks_ac <- rownames(tc[tc$ac > 0, ])

cf_me3 <- ClosestFeature(bmmc, region = peaks_me3)
cf_ac <- ClosestFeature(bmmc, region = peaks_ac)

cf_me3 <- cf_me3[cf_me3$gene_biotype == "protein_coding", ]
cf_ac <- cf_ac[cf_ac$gene_biotype == "protein_coding", ]

maxdist <- 50000
closest_me3 <- unique(cf_me3[cf_me3$distance < maxdist, "gene_name"])
closest_ac <- unique(cf_ac[cf_ac$distance < maxdist, "gene_name"])

# save colorbar
png("plots/bmmc/colorbar_heatmap.png", width = 10, height = 4, units = "in", res = 400)
image(matrix(1:100), col=cols.use)
dev.off()

# RNA
bmmc_rna <- readRDS("objects/fullref.Rds")
ct.keep.rna <- c("HSC", "LMPP", "CLP", "pro B", "pre B", "transitional B", "Naive B", "Memory B")
bmmc_rna_bcell <- bmmc_rna[, bmmc_rna$celltype.l2 %in% ct.keep.rna]

# module score for each set of genes
bmmc_rna_bcell <- AddModuleScore(
  object = bmmc_rna_bcell,
  features = list(closest_me3),
  name = "Repressed"
)
bmmc_rna_bcell <- AddModuleScore(
  object = bmmc_rna_bcell,
  features = list(closest_ac),
  name = "Activated"
)

Idents(bmmc_rna_bcell) <- "celltype.l2"
levels(bmmc_rna_bcell) <- c("HSC", "LMPP", "CLP", "pro B", "pre B", "transitional B", "Naive B", "Memory B")

v1 <- VlnPlot(bmmc_rna_bcell, "Repressed1", pt.size = 0) + NoLegend() +
  geom_hline(yintercept = 0) +
  ylab("Module expression") +
  xlab("") +
  ylim(c(-0.15, 0.15)) +
  ggtitle("Repressed genes") +
  theme_classic() +
  theme(legend.position = "none", text = element_text(size = 20)) +
  scale_x_discrete(guide = guide_axis(angle = 45))

v2 <- VlnPlot(bmmc_rna_bcell, "Activated1", pt.size = 0) + NoLegend() +
  geom_hline(yintercept = 0) +
  ylab("Module expression") +
  xlab("") +
  ylim(c(-0.15, 0.15)) +
  ggtitle("Activated genes") +
  theme_classic() +
  theme(legend.position = "none", text = element_text(size = 20)) +
  scale_x_discrete(guide = guide_axis(angle = 45))

ggsave(filename = "plots/bmmc/rna_repressed.pdf", plot = v1, height = 4, width = 5)
ggsave(filename = "plots/bmmc/rna_activated.pdf", plot = v2, height = 4, width = 5)
