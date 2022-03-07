library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)


args <- commandArgs(trailingOnly = TRUE)
ntt <- readRDS(args[[1]])
atac <- readRDS(args[[2]])
ct_ac <- readRDS(args[[3]])
ct_me <- readRDS(args[[4]])

colormap <- list("H3K27me3" = "#D3145A", "H3K27ac" = "#F98401", "RNAPII" = "#036C9A")
mixed <- list("mp"= "#6b407a", "ma" = "#e64c2e", "ap" = "#7e784e", "all" = "#9a5752")

# QC plots
me3_frag <- Fragments(ntt[['me3']])[[1]]@path
ac_frag <- Fragments(ntt[['ac']])[[1]]@path

tot_me3 <- CountFragments(
  fragments = me3_frag,
  cells = colnames(ntt)
)
tot_me3$mark <- "H3K27me3"
tot_ac <- CountFragments(
  fragments = ac_frag,
  cells = colnames(ntt)
)
tot_ac$mark <- "H3K27ac"

all_counts <- rbind(tot_me3, tot_ac)

mean(tot_me3$frequency_count)
mean(tot_ac$frequency_count)

p0 <- ggplot(all_counts, aes(mark, frequency_count, fill = mark)) +
  geom_violin() +
  scale_y_log10() +
  theme_bw() + 
  ylab("Fragments") +
  xlab("Mark") +
  ggtitle("Total fragments") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c(colormap$H3K27ac, colormap$H3K27me3))
  
ggsave(filename = "plots/pbmc/fragments.png", plot = p0, height = 4, width = 4)

t1 <- TSSPlot(ntt, assay = "me3", group.by = "orig.ident") +
  ggtitle("H3K27me3") + scale_color_manual(values = colormap$H3K27me3)
t2 <- TSSPlot(ntt, assay = "ac", group.by = "orig.ident") +
  ggtitle("H3K27ac") + scale_color_manual(values = colormap$H3K27ac)
tt <- (t1 | t2) & ylim(c(0, 7))
ggsave(filename = "plots/pbmc/tss.png", plot = tt, height = 8, width = 5)

p1 <- DimPlot(ntt, reduction = "umap.me3", label = TRUE, label.size = 3, repel = TRUE, pt.size = 0.1) + theme_void() + ggtitle("H3K27me3") + NoLegend()
p2 <- DimPlot(ntt, reduction = "umap.ac", label = TRUE, label.size = 3, repel = TRUE, pt.size = 0.1) + theme_void() + ggtitle("H3K27ac") + NoLegend()
p3 <- DimPlot(ntt, reduction = "umap.adt", label = TRUE, label.size = 3, repel = TRUE, pt.size = 0.1) + theme_void() + ggtitle("ADT") + NoLegend()
p4 <- DimPlot(ntt, reduction = "umap.wnn", label = TRUE, label.size = 4) + ggtitle("PBMCs: H3K27me3 + H3K27ac + protein") + theme_classic() + NoLegend() + ylab("UMAP 1") + xlab("UMAP 2") 

pp <- ((p1 / p2 / p3) | p4) + plot_layout(widths = c(1, 2))
ggsave(filename = "plots/pbmc/dimplots.png", plot = pp, height = 5, width = 7)

# wnn weight
v1 <- VlnPlot(ntt, c("ac.weight", "me3.weight", "ADT.weight"), ncol = 1, pt.size = 0)
ggsave(filename = "plots/pbmc/wnn_weight.png", plot = v1, height = 8, width = 12)

# adt sensitivity
df <- ntt[[]]
df$dataset <- "scNTT-seq"
df <- df[, c("nCount_ADT", "dataset")]

df2 <- ct_me[[]]
df2$dataset <- "CUT&Tag-pro (H3K27me3)"
df2 <- df2[, c("nCount_ADT", "dataset")]

df3 <- ct_ac[[]]
df3$dataset <- "CUT&Tag-pro (H3K27ac)"
df3 <- df3[, c("nCount_ADT", "dataset")]

df <- rbind(df, df2, df3)

v2 <- ggplot(data = df, aes(dataset, nCount_ADT, fill = dataset)) +
  geom_violin() +
  scale_y_log10() +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(filename = "plots/pbmc/ncount_adt.png", plot = v2, height = 4, width = 6)

# featureplots adt
DefaultAssay(ntt) <- "ADT"
fp <- FeaturePlot(
  ntt,
  c("CD14", "CD19", "CD4", "CD8A", "CD3D", "IL2RB"),
  pt.size = 0.1,
  cols = c("lightgrey", "darkgreen"),
  ncol = 3,
  reduction = 'umap.wnn'
) & xlab("UMAP 1") & ylab("UMAP 2") & NoLegend()

ggsave(filename = "plots/pbmc/featureplots.png", plot = fp, height = 6, width = 8)

DefaultAssay(ntt) <- "me3"
cp <- CoveragePlot(
  object = ntt,
  region = "PAX5",
  assay = list("me3", "ac"),
  split.assays = FALSE,
  assay.scale = "separate",
  features = c("CD19"),
  expression.assay = "ADT",
  idents = c("CD14+ Mono", "B cell"),
  extend.upstream = -30000,
  extend.downstream = 10000,
  window = 500,
  peaks = FALSE
)
cp[[1]][[1]] <- cp[[1]][[1]] + scale_fill_manual(values = c(colormap$H3K27me3, colormap$H3K27ac))
ggsave(filename = "plots/pbmc/pax5_covplot.png", plot = cp, height = 3, width = 8)

# colorscale
cols <- colorRampPalette(c("lightgrey", "darkgreen"))(100)

png("plots/pbmc/colorbar.png", width = 10, height = 4, units = "in", res = 400)
image(matrix(1:100), col=cols)
dev.off()

cp <- CoveragePlot(
  object = ntt,
  region = "CD33",
  assay = list("me3", "ac"),
  split.assays = FALSE,
  assay.scale = "separate",
  features = c("CD33"),
  expression.assay = "ADT",
  idents = c("CD14+ Mono", "B cell"),
  extend.upstream = 2000,
  extend.downstream = 2000,
  window = 500,
  peaks = FALSE
)
cp[[1]][[1]] <- cp[[1]][[1]] + scale_fill_manual(values = c(colormap$H3K27me3, colormap$H3K27ac))
ggsave(filename = "plots/pbmc/cd33_covplot.png", plot = cp, height = 3, width = 8)

# run WNN on pair of chromatin assays, quantify how well structure matches protein neighbor graph

ntt <- FindMultiModalNeighbors(
  object = ntt,
  reduction.list = list("lsi.me3", "lsi.ac"),
  dims.list = list(2:30, 2:30),
  knn.graph.name = "wknn.chrom",
  weighted.nn.name = "wnn.chrom",
  modality.weight.name = "wnn.chrom"
)

# compute KNN purity using ADT-defined cell types, compare different graphs
knn_purity_graph <- function(graph, clusters) {
  # pre-existing graph, can't change k
  nn_purity <- vector(mode = "numeric", length = length(x = clusters))
  for (i in seq_len(length.out = nrow(x = graph))) {
    nbr <- names(which(graph[i, ] == 1))
    nn_purity[i] <- sum(clusters[nbr] == clusters[i]) / sum(graph[i, ])
  }
  return(nn_purity)
}

clusters <- ntt$celltype

ntt <- FindNeighbors(ntt, reduction = "lsi.me3", dims = 2:30, k.param = 100)
ntt <- FindNeighbors(ntt, reduction = "lsi.ac", dims = 2:30, k.param = 100)
ntt <- FindNeighbors(ntt, reduction = "pca", dims = 1:30, k.param = 100)

graph.me3 <- ntt[['me3_nn']]
graph.ac <- ntt[['ac_nn']]
graph.chrom <- ntt[['wknn.chrom']]
graph.adt <- ntt[['ADT_nn']]
graph.wnn <- ntt[['wknn']]

kme3 <- data.frame(
  purity = knn_purity_graph(graph = graph.me3, clusters = clusters),
  assay = "H3K27me3",
  celltype = clusters
)
kac <- data.frame(
  purity = knn_purity_graph(graph = graph.ac, clusters = clusters),
  assay = "H3K27ac",
  celltype = clusters
)
kchrom <- data.frame(
  purity = knn_purity_graph(graph = graph.chrom, clusters = clusters),
  assay = "H3K27me3 + H3K27ac",
  celltype = clusters
)
kadt <- data.frame(
  purity = knn_purity_graph(graph = graph.adt, clusters = clusters),
  assay = "ADT",
  celltype = clusters
)
kwnn <- data.frame(
  purity = knn_purity_graph(graph = graph.wnn, clusters = clusters),
  assay = "ADT + H3K27me3 + H3K27ac",
  celltype = clusters
)

kdf <- rbind(kme3, kac, kchrom, kadt)
kdf$assay <- factor(kdf$assay, levels = c("H3K27ac", "H3K27me3", "H3K27me3 + H3K27ac", "ADT"))

df <- kdf %>% 
  group_by(assay) %>% 
  mutate(total_025 = sum(purity < 0.25) / n()) %>% 
  select(assay, total_025) %>% 
  unique()

knn_p <- ggplot(df, aes(assay, total_025, fill = assay)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("Assay") + ylab("Fraction of cells with \n<25% neighbors same celltype") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_fill_manual(
    values = c(colormap$H3K27ac, colormap$H3K27me3, mixed$ma, "darkgreen")
  )
  
ggsave(filename = "plots/pbmc/knn_purity.pdf", plot = knn_p, height = 4, width = 3)
