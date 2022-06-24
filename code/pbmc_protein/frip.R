library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)

colormap <- list("H3K27me3" = "#D3145A", "H3K27ac" = "#F98401")

pbmc <- readRDS("objects/pbmc_protein.rds")
pbmc_sc <- readRDS("objects/pbmc_sc.rds")
henikoff_me3_filt <- readRDS("objects/henikoff/me3_filt.rds")
henikoff_ac_filt <- readRDS("objects/henikoff/ac_filt.rds")

# ENCODE peaks
pk_ac <- read.table("data/encode/ENCFF832RWT.bed.gz")
pk_me3 <- read.table("data/encode/ENCFF291LVP.bed.gz")

pk_ac <- makeGRangesFromDataFrame(df = pk_ac, seqnames.field = "V1", start.field = "V2", end.field = "V3")
pk_me3 <- makeGRangesFromDataFrame(df = pk_me3, seqnames.field = "V1", start.field = "V2", end.field = "V3")

# total fragments for each dataset

get_total <- function(frags, cells) {
  tot <- CountFragments(fragments = frags, cells = cells)
  totals <- tot$frequency_count
  names(totals) <- tot$CB
  return(totals)
}

# pbmc
totals_pbmc_ac <- get_total(Fragments(pbmc[['ac']])[[1]]@path, colnames(pbmc))
totals_pbmc_me3 <- get_total(Fragments(pbmc[['me3']])[[1]]@path, cells = colnames(pbmc))

# new replicate
totals_pbmc_ac_n <- get_total(Fragments(pbmc_sc[['ac']])[[1]]@path, cells = colnames(pbmc_sc))
totals_pbmc_me3_n <- get_total(Fragments(pbmc_sc[['me3']])[[1]]@path, cells = colnames(pbmc_sc))

# henikoff
r1_cells_filt <- henikoff_me3_filt$orig.ident == '1'
r2_cells_filt <- henikoff_me3_filt$orig.ident == '2'
totals_h_me3_1_f <- get_total(Fragments(henikoff_me3_filt)[[1]]@path, cells = henikoff_me3_filt$barcode[r1_cells_filt])
totals_h_me3_2_f <- get_total(Fragments(henikoff_me3_filt)[[2]]@path, cells = henikoff_me3_filt$barcode[r2_cells_filt])
names(totals_h_me3_1_f) <- paste0(names(totals_h_me3_1_f), ":1")
names(totals_h_me3_2_f) <- paste0(names(totals_h_me3_2_f), ":2")
totals_h_me3_f <- c(totals_h_me3_1_f, totals_h_me3_2_f)

totals_h_ac_f <- get_total(Fragments(henikoff_ac_filt)[[1]]@path, cells = colnames(henikoff_ac_filt))

# total counts in peak regions
get_in_peaks <- function(obj, regions) {
  counts <- FeatureMatrix(
    fragments = Fragments(obj),
    features = regions,
    cells = colnames(obj)
  )
  return(colSums(counts))
}

# ntt
in_peak_pbmc_ac_ac <- get_in_peaks(obj = pbmc[['ac']], regions = pk_ac)
in_peak_pbmc_ac_me <- get_in_peaks(obj = pbmc[['ac']], regions = pk_me3)
in_peak_pbmc_me_me <- get_in_peaks(obj = pbmc[['me3']], regions = pk_me3)
in_peak_pbmc_me_ac <- get_in_peaks(obj = pbmc[['me3']], regions = pk_ac)

# rep2
in_peak_pbmc_ac_ac_n <- get_in_peaks(obj = pbmc_sc[['ac']], regions = pk_ac)
in_peak_pbmc_ac_me_n <- get_in_peaks(obj = pbmc_sc[['ac']], regions = pk_me3)
in_peak_pbmc_me_me_n <- get_in_peaks(obj = pbmc_sc[['me3']], regions = pk_me3)
in_peak_pbmc_me_ac_n <- get_in_peaks(obj = pbmc_sc[['me3']], regions = pk_ac)

# henikoff filt
in_peak_h_ac_ac_f <- get_in_peaks(obj = henikoff_ac_filt, regions = pk_ac)
in_peak_h_ac_me_f <- get_in_peaks(obj = henikoff_ac_filt, regions = pk_me3)
in_peak_h_me_me_f <- get_in_peaks(obj = henikoff_me3_filt, regions = pk_me3)
in_peak_h_me_ac_f <- get_in_peaks(obj = henikoff_me3_filt, regions = pk_ac)

# frip
make_df <- function(totals, inpeak, assay, pk, dataset) {
  d <- data.frame(
    cell = names(totals),
    total = totals,
    in_peak = inpeak[names(totals)],
    assay = assay,
    peak = pk,
    dataset = dataset
  )
  d$FRIP <- d$in_peak / d$total
  return(d)
}

df <- rbind(
  make_df(totals_pbmc_ac, in_peak_pbmc_ac_ac, "H3K27ac", "H3K27ac", "scNTT-seq"),
  make_df(totals_pbmc_ac, in_peak_pbmc_ac_me, "H3K27ac", "H3K27me3", "scNTT-seq"),
  make_df(totals_pbmc_me3, in_peak_pbmc_me_me, "H3K27me3", 'H3K27me3', "scNTT-seq"),
  make_df(totals_pbmc_me3, in_peak_pbmc_me_ac, "H3K27me3", "H3K27ac", "scNTT-seq"),
  
  make_df(totals_h_ac_f, in_peak_h_ac_ac_f, "H3K27ac", "H3K27ac", "scCUT&Tag"),
  make_df(totals_h_ac_f, in_peak_h_ac_me_f, "H3K27ac", "H3K27me3", "scCUT&Tag"),
  make_df(totals_h_me3_f, in_peak_h_me_me_f, "H3K27me3", "H3K27me3", "scCUT&Tag"),
  make_df(totals_h_me3_f, in_peak_h_me_ac_f, "H3K27me3", "H3K27ac", "scCUT&Tag")
)

df2 <- rbind(
  make_df(totals_pbmc_ac_n, in_peak_pbmc_ac_ac_n, "H3K27ac", "H3K27ac", "scNTT-seq"),
  make_df(totals_pbmc_ac_n, in_peak_pbmc_ac_me_n, "H3K27ac", "H3K27me3", "scNTT-seq"),
  make_df(totals_pbmc_me3_n, in_peak_pbmc_me_me_n, "H3K27me3", 'H3K27me3', "scNTT-seq"),
  make_df(totals_pbmc_me3_n, in_peak_pbmc_me_ac_n, "H3K27me3", "H3K27ac", "scNTT-seq"),
)

df_filt <- df[df$assay == df$peak, ]
df_filt <- df_filt[!is.na(df_filt$FRIP), ]

p <- ggplot(df_filt, aes(x = assay, y = FRIP, fill = dataset)) +
  geom_boxplot(outlier.size = 0.5) +
  xlab("Mark") +
  theme_bw() +
  ggtitle("FRiP") +
  theme(legend.position = 'none')

ggsave("plots/pbmc/frip_pbmc.png", plot = p, height = 4, width = 2.5)


p2 <- ggplot(df2, aes(x = assay, y = FRIP, fill = assay)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = colormap) +
  xlab("Mark") +
  theme_bw() +
  ggtitle("FRiP") +
  theme(legend.position = 'none')
ggsave("plots/pbmc/frip_pbmc_rep2.png", plot = p2, height = 4, width = 2.5)

saveRDS(object = df_filt, file = "data/pbmc_protein/frip.rds")

df1 <- df_filt[df_filt$dataset == "scNTT-seq", ]

mean(df1[df1$peak == "H3K27ac", "FRIP"])  # 0.21
sd(df1[df1$peak == "H3K27ac", "FRIP"])  # 0.09
mean(df1[df1$peak == "H3K27me3", "FRIP"])  # 0.11
sd(df1[df1$peak == "H3K27me3", "FRIP"])  # 0.05

mean(df2[df2$peak == "H3K27ac", "FRIP"])  # 0.28
sd(df2[df2$peak == "H3K27ac", "FRIP"])  # 0.08
mean(df2[df2$peak == "H3K27me3", "FRIP"])  # 0.10
sd(df2[df2$peak == "H3K27me3", "FRIP"])  # 0.06


# revision
df3 <- rbind(
  make_df(totals_pbmc_ac, in_peak_pbmc_ac_ac, "H3K27ac", "H3K27ac", "scNTT-seq (PBMC 1)"),
  make_df(totals_pbmc_ac, in_peak_pbmc_ac_me, "H3K27ac", "H3K27me3", "scNTT-seq (PBMC 1)"),
  make_df(totals_pbmc_me3, in_peak_pbmc_me_me, "H3K27me3", 'H3K27me3', "scNTT-seq (PBMC 1)"),
  make_df(totals_pbmc_me3, in_peak_pbmc_me_ac, "H3K27me3", "H3K27ac", "scNTT-seq (PBMC 1)"),
  
  make_df(totals_h_ac_f, in_peak_h_ac_ac_f, "H3K27ac", "H3K27ac", "scCUT&Tag (PBMC)"),
  make_df(totals_h_ac_f, in_peak_h_ac_me_f, "H3K27ac", "H3K27me3", "scCUT&Tag (PBMC)"),
  make_df(totals_h_me3_f, in_peak_h_me_me_f, "H3K27me3", "H3K27me3", "scCUT&Tag (PBMC)"),
  make_df(totals_h_me3_f, in_peak_h_me_ac_f, "H3K27me3", "H3K27ac", "scCUT&Tag (PBMC)"),

  make_df(totals_pbmc_ac_n, in_peak_pbmc_ac_ac_n, "H3K27ac", "H3K27ac", "scNTT-seq (PBMC 2)"),
  make_df(totals_pbmc_ac_n, in_peak_pbmc_ac_me_n, "H3K27ac", "H3K27me3", "scNTT-seq (PBMC 2)"),
  make_df(totals_pbmc_me3_n, in_peak_pbmc_me_me_n, "H3K27me3", 'H3K27me3', "scNTT-seq (PBMC 2)"),
  make_df(totals_pbmc_me3_n, in_peak_pbmc_me_ac_n, "H3K27me3", "H3K27ac", "scNTT-seq (PBMC 2)")
)


bmmc_frip <- readRDS(file = "data/bmmc_dual/frip.rds")
bmmc_frip$dataset <- "scNTT-seq (BMMC)"

df3 <- rbind(bmmc_frip, df3)
df_filt <- df3[df3$assay == df3$peak, ]
df_filt <- df_filt[!is.na(df_filt$FRIP), ]

p <- ggplot(df_filt, aes(x = assay, y = FRIP, fill = dataset)) +
  geom_boxplot(outlier.size = 0.5) +
  xlab("Mark") +
  theme_bw() +
  ggtitle("Fragments in peaks")

ggsave("plots/pbmc/frip_pbmc_revision.png", plot = p, height = 4, width = 6)

