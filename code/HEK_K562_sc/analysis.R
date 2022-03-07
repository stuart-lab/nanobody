library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(Rsamtools)

obj <- readRDS("objects/hek_k562.rds")

# dimplot for figure 1
p_ac <- DimPlot(obj, reduction = "umap.k27ac", pt.size = 0.1) + ggtitle("H3K27ac") + theme_void() + theme(legend.position = "none")
p_me <- DimPlot(obj, reduction = "umap.k27me", pt.size = 0.1) + ggtitle("H3K2me3") + theme_void() + theme(legend.position = "none")
p_pol <- DimPlot(obj, reduction = "umap.pol2", pt.size = 0.1) + ggtitle("RNAPII") + theme_void() + theme(legend.position = "none")
p_wnn <- DimPlot(obj, reduction = "wnn.umap") + ggtitle("H3K27me3 + H3K27ac + RNAPII") + theme_classic() + xlab("UMAP 1") + ylab("UMAP 2") + theme(legend.position = "none")

p1 <- ((p_ac / p_me / p_pol) | p_wnn) + plot_layout(widths = c(1, 2))

ggsave("plots/hek_k562/dimplot_hek_k562.pdf", plot = p1, height = 5, width = 8)

# bulk data for comparison
hek_mono_k27ac <- "data/HEK_K562_bulk/mapped/HEK-mono-K27ac.bw"
hek_mono_k27me3 <- "data/HEK_K562_bulk/mapped/HEK-mono-K27me.bw"
hek_mono_s2s5p <- "data/HEK_K562_bulk/mapped/HEK-mono-Pol2.bw"
k562_mono_k27ac <- "data/HEK_K562_bulk/mapped/K562-mono-K27ac.bw"
k562_mono_k27me3 <- "data/HEK_K562_bulk/mapped/K562-mono-K27me.bw"
k562_mono_s2s5p <- "data/HEK_K562_bulk/mapped/K562-mono-Pol2.bw"

DefaultAssay(obj) <- "K27ac"

p <- CoveragePlot(
  object = obj,
  region = "SCMH1",
  idents = "K562",
  assay = list("K27ac", "K27me", "Pol2"),
  extend.upstream = 20000,
  extend.downstream = 100000,
  window = 50,
  split.assays = TRUE,
  peaks = FALSE,
  bigwig = list("H3K27ac" = k562_mono_k27ac,
                "H3K27me3" = k562_mono_k27me3,
                "RNAPII" = k562_mono_s2s5p)
) & NoLegend()

colors.use <- list("#F98401", "#D3145A", "#036C9A")
mixed <- list("mp"= "#6b407a", "ma" = "#e64c2e", "ap" = "#7e784e", "all" = "#9a5752")

p <- p[[1]]
p[[1]] <- p[[1]] + scale_fill_manual(values = colors.use)
p[[2]] <- p[[2]] & scale_fill_manual(values =  c("darkgrey", "darkgrey", "darkgrey"))

ggsave(
    filename = "plots/SCMH1_k562.png",
    plot = p,
    dpi = 600,
    units = "in",
    height = 5,
    width = 15
)

# additional regions for supplement
p <- CoveragePlot(
  object = obj,
  region = "NASP",
  idents = "K562",
  assay = list("K27ac", "K27me", "Pol2"),
  extend.upstream = 20000,
  extend.downstream = 100000,
  window = 50,
  split.assays = TRUE,
  peaks = FALSE,
  bigwig = list("H3K27ac" = k562_mono_k27ac,
                "H3K27me3" = k562_mono_k27me3,
                "RNAPII" = k562_mono_s2s5p)
) & NoLegend()

colors.use <- list("#F98401", "#D3145A", "#036C9A")
mixed <- list("mp"= "#6b407a", "ma" = "#e64c2e", "ap" = "#7e784e", "all" = "#9a5752")

p <- p[[1]]
p[[1]] <- p[[1]] + scale_fill_manual(values = colors.use)
p[[2]] <- p[[2]] & scale_fill_manual(values =  c("darkgrey", "darkgrey", "darkgrey"))

ggsave(
  filename = "plots/NASP_k562.png",
  plot = p,
  dpi = 600,
  units = "in",
  height = 5,
  width = 15
)

p <- CoveragePlot(
  object = obj,
  region = "FAF1",
  idents = "K562",
  assay = list("K27ac", "K27me", "Pol2"),
  extend.upstream = 50000,
  extend.downstream = 100000,
  window = 50,
  split.assays = TRUE,
  peaks = FALSE,
  bigwig = list("H3K27ac" = k562_mono_k27ac,
                "H3K27me3" = k562_mono_k27me3,
                "RNAPII" = k562_mono_s2s5p)
) & NoLegend()

colors.use <- list("#F98401", "#D3145A", "#036C9A")
mixed <- list("mp"= "#6b407a", "ma" = "#e64c2e", "ap" = "#7e784e", "all" = "#9a5752")

p <- p[[1]]
p[[1]] <- p[[1]] + scale_fill_manual(values = colors.use)
p[[2]] <- p[[2]] & scale_fill_manual(values =  c("darkgrey", "darkgrey", "darkgrey"))

ggsave(
  filename = "plots/FAF1_k562.png",
  plot = p,
  dpi = 600,
  units = "in",
  height = 5,
  width = 15
)

p <- CoveragePlot(
  object = obj,
  region = "SDAD1",
  idents = "K562",
  assay = list("K27ac", "K27me", "Pol2"),
  extend.upstream = 100000,
  extend.downstream = 100000,
  window = 50,
  split.assays = TRUE,
  peaks = FALSE,
  bigwig = list("H3K27ac" = k562_mono_k27ac,
                "H3K27me3" = k562_mono_k27me3,
                "RNAPII" = k562_mono_s2s5p)
) & NoLegend()

colors.use <- list("#F98401", "#D3145A", "#036C9A")
mixed <- list("mp"= "#6b407a", "ma" = "#e64c2e", "ap" = "#7e784e", "all" = "#9a5752")

p <- p[[1]]
p[[1]] <- p[[1]] + scale_fill_manual(values = colors.use)
p[[2]] <- p[[2]] & scale_fill_manual(values =  c("darkgrey", "darkgrey", "darkgrey"))

ggsave(
  filename = "plots/SDAD1_k562.png",
  plot = p,
  dpi = 600,
  units = "in",
  height = 5,
  width = 15
)

# multiCUT&Tag comparison
# same
me3_me3 <- "data/multict/fragments/H3K27me3-H3K27me3.tsv.gz"
ac_ac <- "data/multict/fragments/H3K27ac-H3K27ac.tsv.gz"

# mixed
ac_me3 <- "data/multict/fragments/H3K27ac-H3K27me3.tsv.gz"
me3_ac <- "data/multict/fragments/H3K27me3-H3K27ac.tsv.gz"

totals_mct <- CountFragments(
  fragments = list(me3_me3, ac_ac, ac_me3, me3_ac)
)

cells_keep_mct <- totals_mct[totals_mct$frequency_count > 200, "CB"]

# count fragments for each mark
total_me3_mct <- CountFragments(fragments = me3_me3, cells = cells_keep_mct)
total_ac_mct <- CountFragments(fragments = ac_ac, cells = cells_keep_mct)

# count mixed 1/2 for each
total_me3_mixed_mct <- CountFragments(fragments = me3_ac, cells = cells_keep_mct)
total_ac_mixed_mct <- CountFragments(fragments = ac_me3, cells = cells_keep_mct)
total_me3_mixed_mct$frequency_count <- total_me3_mixed_mct$frequency_count / 2
total_ac_mixed_mct$frequency_count <- total_ac_mixed_mct$frequency_count / 2

mct_counts_me3 <- total_me3_mct$frequency_count
names(mct_counts_me3) <- total_me3_mct$CB
mct_counts_me3[total_me3_mixed_mct$CB] <- mct_counts_me3[total_me3_mixed_mct$CB] + total_me3_mixed_mct$frequency_count
mct_counts_me3[total_ac_mixed_mct$CB] <- mct_counts_me3[total_ac_mixed_mct$CB] + total_ac_mixed_mct$frequency_count

mct_counts_ac <- total_ac_mct$frequency_count
names(mct_counts_ac) <- total_ac_mct$CB
mct_counts_ac[total_ac_mixed_mct$CB] <- mct_counts_ac[total_ac_mixed_mct$CB] + total_ac_mixed_mct$frequency_count
mct_counts_ac[total_me3_mixed_mct$CB] <- mct_counts_ac[total_me3_mixed_mct$CB] + total_me3_mixed_mct$frequency_count

count_df <- data.frame(count = c(mct_counts_me3, mct_counts_ac),
                       mark = c(rep("H3K27me3", length(mct_counts_me3)),
                                rep("H3K27ac", length(mct_counts_ac))),
                       dataset = "multiCUT&Tag"
                       )

# ntt
total_me3 <- CountFragments(
    fragments = "data/HEK_K562_sc/sinto/K27me.bed.gz",
    cells = colnames(obj)
)
total_ac <- CountFragments(
    fragments = "data/HEK_K562_sc/sinto/K27ac.bed.gz",
    cells = colnames(obj)
)
total_pol <- CountFragments(
    fragments = "data/HEK_K562_sc/sinto/Pol2.bed.gz",
    cells = colnames(obj)
)

# construct total vectors
me3 <- total_me3$frequency_count
names(me3) <- total_me3$CB
df1 <- data.frame(
    count = me3,
    mark = "H3K27me3",
    dataset = "scNTT-seq"
)

ac <- total_ac$frequency_count
names(ac) <- total_ac$CB
df2 <- data.frame(
    count = ac,
    mark = "H3K27ac",
    dataset = "scNTT-seq"
)

pol <- total_pol$frequency_count
names(pol) <- total_pol$CB
df3 <- data.frame(
    count = pol,
    mark = "RNAPII",
    dataset = "scNTT-seq"
)

count_df <- rbind(count_df, df1, df2, df3)
count_df <- rbind(count_df, data.frame(count=c(1,1), mark="RNAPII", dataset="multiCUT&Tag"))

# total fragments, split by mark
p3 <- ggplot(data = count_df, mapping = aes(mark, count, fill=dataset)) +
  geom_boxplot(outlier.size = 0.1) +
  scale_y_log10() +
  theme_bw() +
  ylab("Total fragments per cell") +
  xlab("Antibody") +
  ggtitle("Sensitivity") 

ggsave(
    filename = "plots/hek_k562/sensitivity.png",
    plot = p3,
    height = 4,
    width = 5,
    dpi = 500
)

# Region heatmaps
atac_peaks <- rtracklayer::import("data/K562_ATAC/ENCFF006OFA.bigBed")
atac_peaks <- reduce(atac_peaks)
chr.keep <- c(paste0("chr", 1:22), "chrX")
atac_peaks <- keepSeqlevels(atac_peaks, value = chr.keep, pruning.mode = "coarse")
atac_peaks <- subsetByOverlaps(atac_peaks, ranges = blacklist_hg38_unified, invert = TRUE)

regions_me3 <- read.table("data/k562_peaks/ENCFF031FSF.bed.gz", sep = "\t",
                          col.names = c("chromosome", "start", "end", "name", "score",
                                        "strand", "a", "b" , "c", "d"))
regions_me3 <- makeGRangesFromDataFrame(regions_me3, keep.extra.columns = TRUE)
regions_me3$peak <- "H3K27me3"

regions_ac <- read.table("data/k562_peaks/ENCFF038DDS.bed.gz", sep = "\t",
                         col.names = c("chromosome", "start", "end", "name", "score",
                                       "strand", "a", "b" , "c", "d"))
regions_ac <- makeGRangesFromDataFrame(regions_ac, keep.extra.columns = TRUE)
regions_ac$peak <- "H3K27ac"

regions_s5 <- read.table("data/k562_peaks/ENCFF053XYZ.bed.gz", sep = "\t",
                         col.names = c("chromosome", "start", "end", "name", "score",
                                       "strand", "a", "b" , "c", "d"))
regions_s5 <- makeGRangesFromDataFrame(regions_s5, keep.extra.columns = TRUE)

regions_s2 <- read.table("data/k562_peaks/ENCFF266OPF.bed.gz", sep = "\t",
                         col.names = c("chromosome", "start", "end", "name", "score",
                                       "strand", "a", "b" , "c", "d"))
regions_s2 <- makeGRangesFromDataFrame(regions_s2, keep.extra.columns = TRUE)

regions_rna <- reduce(c(regions_s2, regions_s5))
regions_rna$peak <- "RNAPII"

regions <- c(regions_me3, regions_ac, regions_rna)
regions <- keepStandardChromosomes(regions, pruning.mode = "coarse")
regions <- dropSeqlevels(regions, value = "chrM", pruning.mode = "coarse")
regions <- dropSeqlevels(regions, value = "chrY", pruning.mode = "coarse")

atac_chr1 <- atac_peaks[seqnames(atac_peaks) == 'chr1']
me3_peaks <- regions_me3[seqnames(regions_me3) == 'chr1']
ac_peaks <- regions_ac[seqnames(regions_ac) == 'chr1']
pol2_peaks <- regions_rna[seqnames(regions_rna) == 'chr1']

# ATAC
obj <- RegionMatrix(
    object = obj,
    regions = atac_chr1,
    key = "ATAC",
    assay = "K27ac"
)
obj <- RegionMatrix(
    object = obj,
    regions = atac_chr1,
    key = "ATAC",
    assay = "K27me"
)
obj <- RegionMatrix(
    object = obj,
    regions = atac_chr1,
    key = "ATAC",
    assay = "Pol2"
)

# me3
obj <- RegionMatrix(
  object = obj,
  regions = me3_peaks,
  key = "me3",
  assay = "K27ac"
)
obj <- RegionMatrix(
  object = obj,
  regions = me3_peaks,
  key = "me3",
  assay = "K27me"
)
obj <- RegionMatrix(
  object = obj,
  regions = me3_peaks,
  key = "me3",
  assay = "Pol2"
)

# ac
obj <- RegionMatrix(
  object = obj,
  regions = ac_peaks,
  key = "ac",
  assay = "K27ac"
)
obj <- RegionMatrix(
  object = obj,
  regions = ac_peaks,
  key = "ac",
  assay = "K27me"
)
obj <- RegionMatrix(
  object = obj,
  regions = ac_peaks,
  key = "ac",
  assay = "Pol2"
)

# pol2
obj <- RegionMatrix(
  object = obj,
  regions = pol2_peaks,
  key = "pol2",
  assay = "K27ac"
)
obj <- RegionMatrix(
  object = obj,
  regions = pol2_peaks,
  key = "pol2",
  assay = "K27me"
)
obj <- RegionMatrix(
  object = obj,
  regions = pol2_peaks,
  key = "pol2",
  assay = "Pol2"
)

# plots
cm <- list("K27ac" = "#F98401", "K27me" = "#D3145A", "Pol2" = "#036C9A")
p4 <- RegionHeatmap(
  object = obj,
  assay = c("K27ac", "K27me", "Pol2"),
  cols = cm,
  key = "ATAC",
  min.counts = 20,
  normalize = FALSE,
  window = 50,
  nrow = 1,
  idents = "K562"
) & ylab("ATAC peaks")

p5 <- RegionHeatmap(
  object = obj,
  assay = c("K27me", "K27ac", "Pol2"),
  cols = cm,
  key = "me3",
  min.counts = 20,
  normalize = FALSE,
  window = 50,
  nrow = 1,
  idents = "K562"
) & ylab("H3K27me3 peaks")

p6 <- RegionHeatmap(
  object = obj,
  assay = c("K27ac", "K27me", "Pol2"),
  cols = cm,
  key = "ac",
  min.counts = 20,
  normalize = FALSE,
  window = 50,
  nrow = 1,
  idents = "K562"
) & ylab("H3K27ac peaks")

p7 <- RegionHeatmap(
  object = obj,
  assay = c("Pol2", "K27ac", "K27me"),
  cols = cm,
  key = "pol2",
  min.counts = 20,
  normalize = FALSE,
  window = 50,
  nrow = 1,
  idents = "K562"
) & ylab("RNAPII peaks")

update_plot <- function(p, y) {
  p[[1]] <- p[[1]] + ylab(y)
  p[[2]] <- p[[2]] + theme(axis.title.y=element_blank(),
                           axis.text.y=element_blank())
  p[[3]] <- p[[3]] + theme(axis.title.y=element_blank(),
                           axis.text.y=element_blank())
  return(p)
}

p4 <- update_plot(p4, y = "ATAC peaks")
p5 <- update_plot(p5, y = "H3K27me3 peaks")
p6 <- update_plot(p6, y = "H3K27ac peaks")
p7 <- update_plot(p7, y = "RNAPII peaks")

pp <- (p4 | p5 | p6 | p7) & theme(legend.position = "none")
ggsave(filename = "plots/hek_k562/region_heatmaps.png", plot = pp, height = 6, width = 12)

# genome-wide correlation plot

# bulk fragfiles
k562_mono_k27ac_frag <- "data/HEK_K562_bulk/mapped/K562-mono-K27ac.bed.gz"
k562_mono_k27me3_frag <- "data/HEK_K562_bulk/mapped/K562-mono-K27me.bed.gz"
k562_mono_s2s5p_frag <- "data/HEK_K562_bulk/mapped/K562-mono-Pol2.bed.gz"

# split single-cell
sc_k27_ac <- "data/HEK_K562_sc/split/K27ac_k562_filter.bed.gz"
sc_k27_me <- "data/HEK_K562_sc/split/K27me_k562_filter.bed.gz"
sc_pol <- "data/HEK_K562_sc/split/Pol2_k562_filter.bed.gz"

QuantifyRegions <- function(fragfile, regions) {
  results <- vector(mode = "numeric", length = length(regions))
  
  # get total fragments for normalization
  # just read whole thing
  f <- read.table(fragfile, sep = "\t", header = FALSE)
  total <- nrow(f)
  rm(f)
  gc()
  
  # open tabix connection
  tabix.file <- TabixFile(file = fragfile)
  open(con = tabix.file)
  
  # get fragments in each region
  for (j in seq_along(regions)) {
    reads <- scanTabix(file = tabix.file, param = regions[j])
    results[j] <- sum(sapply(X = reads, FUN = length))
  }
  
  # close connection
  close(con = tabix.file)
  cmp <- (results / total) * 10^6
  return(cmp)
}


# need bulk vs bulk
# bulk vs single cell
# single cell vs single cell

k27ac_cov <- QuantifyRegions(k562_mono_k27ac_frag, regions)
k27me3_cov <- QuantifyRegions(k562_mono_k27me3_frag, regions)
pol2_cov <- QuantifyRegions(k562_mono_s2s5p_frag, regions)

pol2_sc_cov <- QuantifyRegions(sc_pol, regions)
ac_sc_cov <- QuantifyRegions(sc_k27_ac, regions)
me_sc_cov <- QuantifyRegions(sc_k27_me, regions)

# create matrix
mat <- cbind(k27ac_cov, ac_sc_cov, k27me3_cov, me_sc_cov, pol2_cov, pol2_sc_cov)

correlation_matrix <- cor(mat)
colnames(correlation_matrix) <- c("Bulk H3K27ac", "sc H3K27ac",
                                  "Bulk H3K27me3", "sc H3K27me3",
                                  "Bulk RNAPII", "sc RNAPII")
rownames(correlation_matrix) <- colnames(correlation_matrix)
df <- as.data.frame(as.table(correlation_matrix))

p8 <- ggplot(df, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu") +
  xlab("") + ylab("") +
  guides(fill = guide_legend(title = "Pearson\ncorrelation")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(filename = "plots/hek_k562/correlation.png", plot = p8, height = 5, width = 6)

# scatterplots

# 1. use encode peak regions for each mark
# 2. plot one mark vs another
# 3. fit linear model

get_cod_expression <- function(m) {
  eq <- substitute(~~italic(R)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 2)))
  cod <- as.expression(eq)
  return(cod)
}

mat <- as.data.frame(mat)
mat$region <- regions$peak
colormap <- list("H3K27me3" = "#D3145A", "H3K27ac" = "#F98401", "RNAPII" = "#036C9A")

# bulk-bulk
m.use <- mat[mat$region %in% c("H3K27me3", "H3K27ac"), ]
bulk_ac_me3 <- ggplot(m.use, aes(k27ac_cov, k27me3_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 400)) + ylim(c(0, 400)) +
  xlab("H3K27ac (bulk, single antibody)") +
  ylab("H3K27me3 (bulk, single antibody)") +
  ggtitle(get_cod_expression(lm(k27ac_cov ~ k27me3_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

m.use <- mat[mat$region %in% c("RNAPII", "H3K27ac"), ]
bulk_ac_pol <- ggplot(m.use,
                      aes(k27ac_cov, pol2_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 300)) + ylim(c(0, 300)) +
  xlab("H3K27ac (bulk, single antibody)") +
  ylab("RNAPII (bulk, single antibody)") +
  ggtitle(get_cod_expression(lm(k27ac_cov ~ pol2_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

m.use <- mat[mat$region %in% c("RNAPII", "H3K27me3"), ]
bulk_pol_me <- ggplot(m.use,
                      aes(pol2_cov, k27me3_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 400)) + ylim(c(0, 400)) +
  xlab("RNAPII (bulk, single antibody)") +
  ylab("H3K27me3 (bulk, single antibody)") +
  ggtitle(get_cod_expression(lm(pol2_cov ~ k27me3_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

# bulk-sc
m.use <- mat[mat$region %in% c("H3K27ac", "H3K27me3"), ]
bulk_ac_sc_me3 <- ggplot(m.use,
                         aes(k27ac_cov, me_sc_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 800)) + ylim(c(0, 800)) +
  xlab("H3K27ac (bulk, single antibody)") +
  ylab("H3K27me3 (single-cell, multiplexed)") +
  ggtitle(get_cod_expression(lm(k27ac_cov ~ me_sc_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

m.use <- mat[mat$region %in% c("H3K27ac", "RNAPII"), ]
bulk_ac_sc_pol <- ggplot(m.use,
                         aes(k27ac_cov, pol2_sc_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 400)) + ylim(c(0, 400)) +
  xlab("H3K27ac (bulk, single antibody)") +
  ylab("RNAPII (single-cell, multiplexed)") +
  ggtitle(get_cod_expression(lm(k27ac_cov ~ pol2_sc_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

m.use <- mat[mat$region %in% c("H3K27me3", "RNAPII"), ]
bulk_pol_sc_me <- ggplot(m.use,
                         aes(pol2_cov, me_sc_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 400)) + ylim(c(0, 400)) +
  xlab("RNAPII (bulk, single antibody)") +
  ylab("H3K27me3 (single-cell, multiplexed)") +
  ggtitle(get_cod_expression(lm(pol2_cov ~ me_sc_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

m.use <- mat[mat$region %in% c("H3K27ac"), ]
bulk_ac_sc_ac <- ggplot(m.use,
                        aes(k27ac_cov, ac_sc_cov, color = region)) +
  geom_point(size = 0.1) +
  theme_classic() +
  scale_color_manual(values = colormap) +
  xlim(c(0, 600)) + ylim(c(0, 600)) +
  xlab("H3K27ac (bulk, single antibody)") +
  ylab("H3K27ac (single-cell, multiplexed)") +
  ggtitle(get_cod_expression(lm(k27ac_cov ~ ac_sc_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

m.use <- mat[mat$region %in% c("H3K27me3"), ]
bulk_me_sc_me <- ggplot(m.use,
                        aes(k27me3_cov, me_sc_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 300)) + ylim(c(0, 300)) +
  xlab("H3K27me3 (bulk, single antibody)") +
  ylab("H3K27me3 (single-cell, multiplexed)") +
  ggtitle(get_cod_expression(lm(k27me3_cov ~ me_sc_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

m.use <- mat[mat$region %in% c("RNAPII"), ]
bulk_pol_sc_pol <- ggplot(m.use,
                          aes(pol2_cov, pol2_sc_cov, color = region)) +
  geom_point(size = 0.1) +
  theme_classic() +
  scale_color_manual(values = colormap) +
  xlim(c(0, 200)) + ylim(c(0, 200)) +
  xlab("RNAPII (bulk, single antibody)") +
  ylab("RNAPII (single-cell, multiplexed)") +
  ggtitle(get_cod_expression(lm(pol2_cov ~ pol2_sc_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

# sc-sc
m.use <- mat[mat$region %in% c("H3K27ac", "H3K27me3"), ]
sc_ac_sc_me3 <- ggplot(m.use,
                       aes(ac_sc_cov, me_sc_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 400)) + ylim(c(0, 400)) +
  xlab("H3K27ac (single-cell, multiplexed)") +
  ylab("H3K27me3 (single-cell, multiplexed)") +
  ggtitle(get_cod_expression(lm(ac_sc_cov ~ me_sc_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

m.use <- mat[mat$region %in% c("H3K27ac", "RNAPII"), ]
sc_ac_sc_pol <- ggplot(m.use,
                       aes(ac_sc_cov, pol2_sc_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 400)) + ylim(c(0, 400)) +
  xlab("H3K27ac (single-cell, multiplexed)") +
  ylab("RNAPII (single-cell, multiplexed)") +
  ggtitle(get_cod_expression(lm(ac_sc_cov ~ pol2_sc_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

m.use <- mat[mat$region %in% c("H3K27me3", "RNAPII"), ]
sc_pol_sc_me <- ggplot(m.use,
                       aes(pol2_sc_cov, me_sc_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 200)) + ylim(c(0, 200)) +
  xlab("RNAPII (single-cell, multiplexed)") +
  ylab("H3K27me3 (single-cell, multiplexed)") +
  ggtitle(get_cod_expression(lm(pol2_sc_cov ~ me_sc_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

pp <- (bulk_ac_me3 / bulk_ac_pol / bulk_pol_me) | (bulk_ac_sc_me3 / bulk_ac_sc_pol / bulk_pol_sc_me) | (sc_ac_sc_me3 / sc_ac_sc_pol / sc_pol_sc_me) | (bulk_me_sc_me / bulk_ac_sc_ac / bulk_pol_sc_pol)
pp <- pp & NoLegend()
ggsave(filename = "plots/hek_k562/scatterplots.png", plot = pp, height = 10, width = 12)

# ternary plot for single-cell data
if (!requireNamespace(package = "ggtern", quietly = TRUE)) {
  install.packages("ggtern")
}
library(ggtern)

region_counts <- mat[, c("region", "ac_sc_cov", "me_sc_cov", "pol2_sc_cov")]
colnames(region_counts) <- c("region", "H3K27ac", "H3K27me3", "RNAPII")

# randomize plotting order
region_counts <- region_counts[sample(1:nrow(region_counts)), ]

gp <- ggtern(region_counts, aes(x = H3K27ac, y = H3K27me3, z = RNAPII, color = region)) +
  geom_point(size = 0.1, alpha = 0.2) +
  scale_color_manual(values = colors.use) +
  theme_classic()

ggsave(filename = "plots/hek_k562/ternary.png", plot = gp, height = 8, width = 8)

## KNN purity ## 
if (!requireNamespace(package = "RANN", quietly = TRUE)) {
  install.packages("RANN")
}
library(RANN)

knn_purity_emb <- function(embeddings, dims, clusters, k) {
  nn <- nn2(data = embeddings[, dims], k = k + 1)$nn.idx[, 2:k]
  nn_purity <- vector(mode = "numeric", length = length(x = clusters))
  for (i in seq_len(length.out = nrow(x = nn))) {
    nn_purity[i] <- sum(clusters[nn[i, ]] == clusters[i]) / k
  }
  return(nn_purity)
}

knn_purity_graph <- function(graph, clusters) {
  # pre-existing graph, can't change k
  nn_purity <- vector(mode = "numeric", length = length(x = clusters))
  for (i in seq_len(length.out = nrow(x = graph))) {
    nbr <- names(which(graph[i, ] == 1))
    nn_purity[i] <- sum(clusters[nbr] == clusters[i]) / sum(graph[i, ])
  }
  return(nn_purity)
}

# compute wnn for each combination of modalities
obj <- FindMultiModalNeighbors(
  object = obj,
  knn.graph.name = "wknn_all",
  reduction.list = list("lsi.k27ac", "lsi.k27me", "lsi.pol2"), 
  dims.list = list(2:10, 2:10, 2:10)
)

obj <- FindMultiModalNeighbors(
  object = obj,
  knn.graph.name = "wknn_ac_me3",
  reduction.list = list("lsi.k27ac", "lsi.k27me"), 
  dims.list = list(2:10, 2:10)
)

obj <- FindMultiModalNeighbors(
  object = obj,
  knn.graph.name = "wknn_ac_pol2",
  reduction.list = list("lsi.k27ac", "lsi.pol2"), 
  dims.list = list(2:10, 2:10)
)

obj <- FindMultiModalNeighbors(
  object = obj,
  knn.graph.name = "wknn_me3_pol2",
  reduction.list = list("lsi.k27me", "lsi.pol2"), 
  dims.list = list(2:10, 2:10)
)

dims <- 2:10
k <- 50
emb_me3 <- Embeddings(obj[['lsi.k27me']])
emb_ac <- Embeddings(obj[['lsi.k27ac']])
emb_pol <- Embeddings(obj[['lsi.pol2']])

ct <- obj$celltype

kme3 <- data.frame(
  purity = knn_purity_emb(embeddings = emb_me3, dims = dims, clusters = ct, k = k),
  assay = "H3K27me3",
  celltype = ct
)
kac <- data.frame(
  purity = knn_purity_emb(embeddings = emb_ac, dims = dims, clusters = ct, k = k),
  assay = "H3K27ac",
  celltype = ct
)
kpol <- data.frame(
  purity = knn_purity_emb(embeddings = emb_pol, dims = dims, clusters = ct, k = k),
  assay = "RNAPII",
  celltype = ct
)

kme3_ac <- data.frame(
  purity = knn_purity_graph(graph = obj[['wknn_ac_me3']], clusters = ct),
  assay = "H3K27me3 + H3K27ac",
  celltype = ct
)
kac_pol <- data.frame(
  purity = knn_purity_graph(graph = obj[['wknn_ac_pol2']], clusters = ct),
  assay = "H3K27ac + RNAPII",
  celltype = ct
)
kpol_me <- data.frame(
  purity = knn_purity_graph(graph = obj[['wknn_me3_pol2']], clusters = ct),
  assay = "H3K27me3 + RNAPII",
  celltype = ct
)
kpol_me_ac <- data.frame(
  purity = knn_purity_graph(graph = obj[['wknn_all']], clusters = ct),
  assay = "H3K27me3 + H3K27ac + RNAPII",
  celltype = ct
)

kdf <- rbind(kme3, kac, kpol, kme3_ac, kac_pol, kpol_me, kpol_me_ac)
kdf$assay <- factor(kdf$assay, levels = c("RNAPII", "H3K27ac", "H3K27me3",
                                             "H3K27ac + RNAPII", "H3K27me3 + RNAPII", "H3K27me3 + H3K27ac",
                                             "H3K27me3 + H3K27ac + RNAPII"))

knn_p <- ggplot(kdf, aes(y = assay, x = purity, fill = assay)) +
  ggridges::geom_density_ridges(bandwidth = 0.008) +
  xlab("Fraction of neighbors of same cell type") + ylab("Assay") +
  theme_classic() +
  scale_fill_manual(values = c("#036C9A", "#F98401", "#D3145A", "#7e784e", "#6b407a", "#e64c2e", "#9a5752")) +
  xlim(c(0.8, 1.01)) +
  theme(legend.position = "none", text = element_text(size = 12)) + 
  geom_vline(xintercept = 1)

ggsave(filename = "plots/hek_k562/knn_purity_wnn.png", plot = knn_p, height = 4, width = 7)