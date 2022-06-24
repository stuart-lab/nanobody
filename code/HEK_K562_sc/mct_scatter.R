library(Signac)
library(ggplot2)
library(patchwork)
library(Rsamtools)


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

get_cod_expression <- function(m) {
  eq <- substitute(~~italic(R)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 2)))
  cod <- as.expression(eq)
  return(cod)
}

colormap <- list("H3K27me3" = "#D3145A", "H3K27ac" = "#F98401")

# multiCUT&Tag
# same
me3_me3 <- "data/multict/fragments/H3K27me3-H3K27me3.tsv.gz"
ac_ac <- "data/multict/fragments/H3K27ac-H3K27ac.tsv.gz"

# mixed
ac_me3 <- "data/multict/fragments/H3K27ac-H3K27me3.tsv.gz"
me3_ac <- "data/multict/fragments/H3K27me3-H3K27ac.tsv.gz"

# peaks
h3k27me3 <- read.table("data/mesc/ENCFF008XKX.bed.gz", sep = "\t") # H3K27me3 ES-Bruce4 mm10
h3k27ac <- read.table("data/mesc/ENCFF360VIS.bed.gz", sep = "\t") # H3K27ac ES-Bruce4 mm10

h3k27me3 <- makeGRangesFromDataFrame(
  df = h3k27me3, seqnames.field = "V1", start.field = "V2", end.field = "V3"
)
h3k27ac <- makeGRangesFromDataFrame(
  df = h3k27ac, seqnames.field = "V1", start.field = "V2", end.field = "V3"
)

h3k27me3 <- keepStandardChromosomes(h3k27me3, pruning.mode = "coarse")
h3k27ac <- keepStandardChromosomes(h3k27ac, pruning.mode = "coarse")

totals_mct <- CountFragments(
  fragments = list(me3_me3, ac_ac, ac_me3, me3_ac)
)

cells_keep_mct <- totals_mct[totals_mct$frequency_count > 200, "CB"]

me3_me3_filt <- "data/multict/fragments/H3K27me3-H3K27me3_filtered.tsv.gz"
ac_ac_filt <- "data/multict/fragments/H3K27ac-H3K27ac_filtered.tsv.gz"
ac_me3_filt <- "data/multict/fragments/H3K27ac-H3K27me3_filtered.tsv.gz"
me3_ac_filt <- "data/multict/fragments/H3K27me3-H3K27ac_filtered.tsv.gz"

# filter cells
FilterCells(
  fragments = me3_me3,
  cells = cells_keep_mct,
  outfile = me3_me3_filt
)
FilterCells(
  fragments = ac_ac,
  cells = cells_keep_mct,
  outfile = ac_ac_filt
)
FilterCells(
  fragments = ac_me3,
  cells = cells_keep_mct,
  outfile = ac_me3_filt
)
FilterCells(
  fragments = me3_ac,
  cells = cells_keep_mct,
  outfile = me3_ac_filt
)

h3k27me3$peak <- "H3K27me3"
h3k27ac$peak <- "H3K27ac"
regions <- c(h3k27me3, h3k27ac)

###########
# scatter #
###########

# quantify regions
me3_me3 <- QuantifyRegions(me3_me3_filt, regions)
me3_ac <- QuantifyRegions(me3_ac_filt, regions)
ac_ac <- QuantifyRegions(ac_ac_filt, regions)
ac_me3 <- QuantifyRegions(ac_me3_filt, regions)

mat <- cbind(me3_me3, me3_ac, ac_ac, ac_me3)
mat <- as.data.frame(mat)
mat$region <- regions$peak

mat$me3 <- mat$me3_me3 + (mat$me3_ac/2) + (mat$ac_me3/2)
mat$ac <- mat$ac_ac + (mat$me3_ac/2) + (mat$ac_me3/2)

p_me3_ac <- ggplot(mat, aes(x = ac, y = me3, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  theme(legend.position = 'none') +
  xlim(c(0, 400)) + ylim(c(0, 400)) +
  ylab("H3K27me3") +
  xlab("H3K27ac") +
  ggtitle(get_cod_expression(lm(me3 ~ ac, data = mat))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

ggsave(filename = "plots/hek_k562/multi_cuttag_scatter.png", plot = p_me3_ac, height = 4, width = 4)

############
### FRiP ###
############

get_total <- function(frags, cells) {
  tot <- CountFragments(fragments = frags@path, cells = cells)
  totals <- tot$frequency_count
  names(totals) <- tot$CB
  return(totals)
}

get_in_peaks <- function(obj, regions) {
  counts <- FeatureMatrix(
    fragments = obj,
    features = regions,
    cells = obj@cells
  )
  return(colSums(counts))
}


me3_me3_obj <- CreateFragmentObject(me3_me3_filt, cells = cells_keep_mct)
ac_ac_obj <- CreateFragmentObject(ac_ac_filt, cells = cells_keep_mct)
me3_ac_obj <- CreateFragmentObject(me3_ac_filt, cells = cells_keep_mct)
ac_me3_obj <- CreateFragmentObject(ac_me3_filt, cells = cells_keep_mct)

totals_me3_me3 <- get_total(me3_me3_obj, cells = cells_keep_mct)
totals_ac_me3 <- get_total(ac_me3_obj, cells = cells_keep_mct)
totals_me3_ac <- get_total(me3_ac_obj, cells = cells_keep_mct)
totals_ac_ac <- get_total(ac_ac_obj, cells = cells_keep_mct)

in_me3_peak_me3_me3 <- get_in_peaks(obj = me3_me3_obj, regions = h3k27me3)
in_me3_peak_ac_ac <- get_in_peaks(obj = ac_ac_obj, regions = h3k27me3)
in_me3_peak_me3_ac <- get_in_peaks(obj = me3_ac_obj, regions = h3k27me3)
in_me3_peak_ac_me3 <- get_in_peaks(obj = ac_me3_obj, regions = h3k27me3)

in_ac_peak_me3_me3 <- get_in_peaks(obj = me3_me3_obj, regions = h3k27ac)
in_ac_peak_ac_ac <- get_in_peaks(obj = ac_ac_obj, regions = h3k27ac)
in_ac_peak_me3_ac <- get_in_peaks(obj = me3_ac_obj, regions = h3k27ac)
in_ac_peak_ac_me3 <- get_in_peaks(obj = ac_me3_obj, regions = h3k27ac)

me3_frip <- data.frame(
  cell = names(totals_me3_me3),
  total = totals_me3_me3 + (totals_ac_me3/2) + (totals_me3_ac/2),
  in_peak = in_me3_peak_me3_me3 + (in_me3_peak_me3_ac/2) + (in_me3_peak_ac_me3/2),
  assay = "H3K27me3",
  peak = "H3K27me3"
)
me3_frip$FRiP <- me3_frip$in_peak / me3_frip$total

ac_frip <- data.frame(
  cell = names(totals_ac_ac),
  total = totals_ac_ac + (totals_ac_me3/2) + (totals_me3_ac/2),
  in_peak = in_ac_peak_ac_ac + (in_ac_peak_me3_ac/2) + (in_ac_peak_ac_me3/2),
  assay = "H3K27ac",
  peak = "H3K27ac"
)
ac_frip$FRiP <- ac_frip$in_peak / ac_frip$total

df <- rbind(me3_frip, ac_frip)

p <- ggplot(df, aes(x = assay, y = FRiP, fill = assay)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = colormap) +
  xlab("Antibody") +
  ggtitle("multiCUT&Tag FRiP") +
  theme_bw() +
  theme(legend.position = "none")

ggsave(filename = "plots/hek_k562/multi_cuttag_frip.png", plot = p, height = 4, width = 3)

saveRDS(object = df, file = "objects/frip_mct.rds")
