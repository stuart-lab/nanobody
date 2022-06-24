library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)

obj <- readRDS("objects/hek_k562.rds")
obj <- obj[, obj$celltype == "K562"]

get_total <- function(frags, cells) {
  tot <- CountFragments(fragments = frags, cells = cells)
  totals <- tot$frequency_count
  names(totals) <- tot$CB
  return(totals)
}

get_in_peaks <- function(obj, regions) {
  counts <- FeatureMatrix(
    fragments = Fragments(obj),
    features = regions,
    cells = colnames(obj)
  )
  return(colSums(counts))
}

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

# frip
totals_ac <- get_total(Fragments(obj[['K27ac']])[[1]]@path, cells = colnames(obj))
totals_me <- get_total(Fragments(obj[['K27me']])[[1]]@path, cells = colnames(obj))
totals_rna <- get_total(Fragments(obj[['Pol2']])[[1]]@path, cells = colnames(obj))

in_peak_ac <- get_in_peaks(obj = obj[['K27ac']], regions = regions_ac)
in_peak_me <- get_in_peaks(obj = obj[['K27me']], regions = regions_me3)
in_peak_rna <- get_in_peaks(obj = obj[['Pol2']], regions = regions_rna)

make_df <- function(totals, inpeak, assay, pk, dataset) {
  d <- data.frame(
    cell = names(totals),
    total = totals,
    in_peak = inpeak[names(totals)],
    assay = assay,
    peak = pk,
    dataset = dataset
  )
  d$FRiP <- d$in_peak / d$total
  return(d)
}

df <- rbind(
  make_df(totals_ac, in_peak_ac, "H3K27ac", "H3K27ac", "NTT"),
  make_df(totals_me, in_peak_me, "H3K27me3", "H3K27me3", "NTT"),
  make_df(totals_rna, in_peak_rna, "RNAPII", 'RNAPII', "NTT")
)

p <- ggplot(df, aes(x = assay, y = FRiP, fill = assay)) +
  geom_boxplot() +
  theme_bw()
ggsave("plots/hek_k562/frip_k562.png", plot = p, height = 4, width = 6)

# load mct
mct_frip <- readRDS("objects/frip_mct.rds")
df$dataset <- NULL
df$dataset <- "scNTT-seq"
mct_frip$dataset <- "multiCUT&Tag"

df_all <- rbind(df, mct_frip)

p <- ggplot(df_all[df_all$assay != "RNAPII", ], aes(x = assay, y = FRiP, fill = dataset)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  xlab("Assay") +
  ggtitle("FRiP") +
  theme(legend.position = 'none')
ggsave("plots/hek_k562/mct_vs_ntt_frip.png", plot = p, height = 4, width = 3)

# total fragments
p <- ggplot(df, aes(x = assay, y = total, fill = assay)) +
  geom_violin() +
  ylab("Total fragments") +
  xlab("Assay") +
  scale_y_log10() +
  theme_bw()

ggsave("plots/hek_k562/total_frags_k562.png", plot = p, height = 4, width = 6)
saveRDS(object = df, file = "data/HEK_K562_sc/frip.rds")

mean(df[df$peak == "H3K27ac", "FRIP"])  # 0.59
sd(df[df$peak == "H3K27ac", "FRIP"])  # 0.09
mean(df[df$peak == "H3K27me3", "FRIP"])  # 0.40
sd(df[df$peak == "H3K27me3", "FRIP"])  # 0.05
mean(df[df$peak == "RNAPII", "FRIP"])  # 0.20
sd(df[df$peak == "RNAPII", "FRIP"])  # 0.04
