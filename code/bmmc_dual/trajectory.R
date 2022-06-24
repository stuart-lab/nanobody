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
plan("multicore", workers = 12)
options(future.globals.maxSize = Inf)


bmmc <- readRDS("objects/bmmc_dual.rds")
Idents(bmmc) <- "celltype"
cells.use <- c("HSPC", "GMP/CLP", "Pre-B", "B", "Plasma")

bcell <- bmmc[, bmmc$celltype %in% cells.use]

DefaultAssay(bcell) <- "me3"
feat.keep <- names(which(rowSums(bcell, slot = 'counts') > 1)) # remove low-count features
bcell[['me3']] <- subset(bcell[['me3']], features = feat.keep)
bcell <- FindTopFeatures(bcell, min.cutoff = 2)
bcell <- RunTFIDF(bcell, scale.factor = median(bcell$nCount_me3))
bcell <- RunSVD(bcell, reduction.name = "lsi.me3", scale.embeddings = FALSE)
bcell <- RunUMAP(bcell, reduction = "lsi.me3", reduction.name = 'umap.me3', dims = 2:20)

DefaultAssay(bcell) <- "ac"
feat.keep <- names(which(rowSums(bcell, slot = 'counts') > 1)) # remove low-count features
bcell[['ac']] <- subset(bcell[['ac']], features = feat.keep)
bcell <- FindTopFeatures(bcell, min.cutoff = 2)
bcell <- RunTFIDF(bcell, scale.factor = median(bcell$nCount_ac))
bcell <- RunSVD(bcell, reduction.name = "lsi.ac", scale.embeddings = FALSE)
bcell <- RunUMAP(bcell, reduction = "lsi.ac", reduction.name = 'umap.ac', dims = 2:20)

# WNN
bcell <- FindMultiModalNeighbors(
  object = bcell,
  reduction.list = list('lsi.me3', 'lsi.ac'),
  dims.list = list(2:20, 2:10)
)

bcell <- RunUMAP(bcell, nn.name = "weighted.nn", reduction.name = "umap.wnn")

## run monocle3
bcell@reductions$umap <- bcell@reductions$umap.wnn # needs to be named UMAP for monocle
bcell.cds <- as.cell_data_set(bcell)
bcell.cds <- cluster_cells(cds = bcell.cds, reduction_method = "UMAP")
bcell.cds <- learn_graph(bcell.cds, use_partition = FALSE)
root_cells <- WhichCells(bcell, idents = "HSPC")
bcell.cds <- order_cells(bcell.cds, reduction_method = "UMAP", root_cells = tail(root_cells, 10))

bcell <- AddMetaData(
  object = bcell,
  metadata = bcell.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "B.cell.pseudotime"
)

fp1 <- FeaturePlot(bcell, "B.cell.pseudotime", reduction = "umap.wnn") +
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
saveRDS(bins_me3, "bins_me3_bcell.rds")
bins_me3 <- readRDS("bins_me3_bcell.rds")

bins_ac <- GenomeBinMatrix(
  fragments = Fragments(bcell[['ac']]),
  genome = chr.use,
  binsize = 10000,
  cells = colnames(bcell)
)
saveRDS(bins_ac, "bins_ac_bcell.rds")
bins_ac <- readRDS("bins_ac_bcell.rds")

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

pv_me3 <- apply(smoothed_me3, 1, function(x) {
  cor.test(x, pstime)$p.value
})
pv_ac <- apply(smoothed_ac, 1, function(x) {
  cor.test(x, pstime)$p.value
})

cor_df <- data.frame(me3 = peak_cor_me3, ac = peak_cor_ac,
                     p_me3 = pv_me3, p_adj_me3 = p.adjust(pv_me3),
                     p_ac = pv_ac, p_adj_ac = p.adjust(pv_ac))
cor_df$diff <- cor_df$me3 - cor_df$ac

# cor 0.2 and difference 0.5
tc <- cor_df[abs(cor_df$diff) > 0.5 & abs(cor_df$me3) > 0.2 & abs(cor_df$ac) > 0.2, ]

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
ct.keep.rna <- c("HSC", "LMPP", "CLP", "pro B", "pre B", "transitional B", "Naive B", "Memory B", "Plasma")
bmmc_rna_bcell <- bmmc_rna[, bmmc_rna$celltype.l2 %in% ct.keep.rna]

# module score for each set of genes
bmmc_rna_bcell <- AddModuleScore(
  object = bmmc_rna_bcell,
  features = list(setdiff(closest_me3, closest_ac)),
  name = "Repressed"
)
bmmc_rna_bcell <- AddModuleScore(
  object = bmmc_rna_bcell,
  features = list(setdiff(closest_ac, closest_me3)),
  name = "Activated"
)

Idents(bmmc_rna_bcell) <- "celltype.l2"
levels(bmmc_rna_bcell) <- c("HSC", "LMPP", "CLP", "pro B", "pre B", "transitional B", "Naive B", "Memory B", "Plasma")

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

# add DE test between module scores at start vs end of trajectory
module_activated <- bmmc_rna_bcell$Activated1
module_repressed <- bmmc_rna_bcell$Repressed1

progenitor <- colnames(bmmc_rna_bcell)[bmmc_rna_bcell$celltype.l2 %in% c("HSC", "LMPP")]
mature <- colnames(bmmc_rna_bcell)[bmmc_rna_bcell$celltype.l2 %in% c("Memory B", "Plasma")]

t.test(module_activated[progenitor], module_activated[mature]) # <2.2e-16
t.test(module_repressed[progenitor], module_repressed[mature]) # <2.2e-16

########################
# BMMC atac comparison
########################

bmmc_atac <- readRDS("objects/bmmc_atac.rds")

# subset b cell lineage
ct.keep <- c("01_HSC", "06_CLP.1", "15_CLP.2", "05_CMP.LMPP",
             "17_B", "18_Plasma")
b_atac <- bmmc_atac[, bmmc_atac$BioClassification %in% ct.keep]

DefaultAssay(b_atac) <- "ATAC"
b_atac <- FindTopFeatures(b_atac, min.cutoff = 5)
b_atac <- RunTFIDF(b_atac)
b_atac <- RunSVD(b_atac)
b_atac <- RunUMAP(b_atac, reduction = "lsi", dims = 2:20)
b_atac <- FindNeighbors(b_atac, reduction = "lsi", dims = 2:20)
knn_atac <- b_atac[['ATAC_nn']]

# add numeric score for b trajectory
b_atac$trajectory <- as.numeric(factor(b_atac$BioClassification,
                                       levels = c("01_HSC",
                                                  "05_CMP.LMPP",
                                                  "06_CLP.1",
                                                  "15_CLP.2",
                                                  "17_B",
                                                  "18_Plasma")))

b_atac_traj <- b_atac$trajectory

b_atac$b <- factor(b_atac$BioClassification,
                                       levels = c("01_HSC",
                                                  "05_CMP.LMPP",
                                                  "06_CLP.1",
                                                  "15_CLP.2",
                                                  "17_B",
                                                  "18_Plasma"))


Idents(b_atac) <- "b"
levels(bcell) <- c('HSPC', 'GMP/CLP', 'Pre-B', 'B', 'Plasma')

# check dynamic regions in atac dataset
tc_regions_atac <- FeatureMatrix(
  fragments = Fragments(b_atac),
  features = rownames(tc),
  cells = colnames(b_atac)
)
smoothed_atac <- knn_smooth(counts = tc_regions_atac[, colnames(knn_atac)], knn = knn_atac)
smoothed_atac <- t(t(smoothed_atac) / colSums(smoothed_atac))
peak_cor_atac <- cor(t(as.matrix(smoothed_atac)), b_atac_traj)[, 1]

pv_atac <- apply(smoothed_atac, 1, function(x) {
  cor.test(x, b_atac_traj)$p.value
})

nc <- peak_cor_atac[abs(peak_cor_atac) < 0.05]
ns <- p.adjust(pv_atac[names(nc)]) > 0.01

sum(ns)
