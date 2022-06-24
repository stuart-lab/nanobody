library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)


bmmc <- readRDS("objects/bmmc_dual.rds")

# ENCODE peaks
pk_ac <- read.table("data/encode/ENCFF832RWT.bed.gz")
pk_me3 <- read.table("data/encode/ENCFF291LVP.bed.gz")

pk_ac <- makeGRangesFromDataFrame(df = pk_ac, seqnames.field = "V1", start.field = "V2", end.field = "V3")
pk_me3 <- makeGRangesFromDataFrame(df = pk_me3, seqnames.field = "V1", start.field = "V2", end.field = "V3")

t_bmmc_ac <- CountFragments(fragments = Fragments(bmmc[['ac']])[[1]]@path, cells = colnames(bmmc))
totals_bmmc_ac <- t_bmmc_ac$frequency_count
names(totals_bmmc_ac) <- t_bmmc_ac$CB

t_bmmc_me3 <- CountFragments(fragments = Fragments(bmmc[['me3']])[[1]]@path, cells = colnames(bmmc))
totals_bmmc_me3 <- t_bmmc_me3$frequency_count
names(totals_bmmc_me3) <- t_bmmc_me3$CB

# total counts in peak regions
get_total <- function(obj, regions) {
  counts <- FeatureMatrix(
    fragments = Fragments(obj),
    features = regions,
    cells = colnames(obj)
  )
  return(colSums(counts))
}

# ntt
in_peak_bmmc_ac_ac <- get_total(obj = bmmc[['ac']], regions = pk_ac)
in_peak_bmmc_ac_me <- get_total(obj = bmmc[['ac']], regions = pk_me3)
in_peak_bmmc_me_me <- get_total(obj = bmmc[['me3']], regions = pk_me3)
in_peak_bmmc_me_ac <- get_total(obj = bmmc[['me3']], regions = pk_ac)

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
  make_df(totals_bmmc_ac, in_peak_bmmc_ac_ac, "H3K27ac", "H3K27ac", "NTT"),
  make_df(totals_bmmc_ac, in_peak_bmmc_ac_me, "H3K27ac", "H3K27me3", "NTT"),
  make_df(totals_bmmc_me3, in_peak_bmmc_me_me, "H3K27me3", 'H3K27me3', "NTT"),
  make_df(totals_bmmc_me3, in_peak_bmmc_me_ac, "H3K27me3", "H3K27ac", "NTT")
)

df_filt <- df[df$assay == df$peak, ]

p <- ggplot(df_filt, aes(x = peak, y = FRIP, fill = peak)) +
  geom_boxplot(outlier.size = 0.2) +
  scale_fill_manual(values = c("#F98401", "#D3145A")) +
  theme_bw() +
  ylab("FRiP") +
  xlab("Assay") +
  ggtitle("FRiP") +
  theme(legend.position = "none")

ggsave("plots/bmmc/frip_bmmc.png", plot = p, height = 4, width = 3)
saveRDS(object = df_filt, file = "data/bmmc_dual/frip.rds")

mean(df_filt[df_filt$peak == "H3K27me3", "FRIP" ]) # 0.18
sd(df_filt[df_filt$peak == "H3K27me3", "FRIP" ]) # 0.09
mean(df_filt[df_filt$peak == "H3K27ac", "FRIP" ]) # 0.26
sd(df_filt[df_filt$peak == "H3K27ac", "FRIP" ]) # 0.09