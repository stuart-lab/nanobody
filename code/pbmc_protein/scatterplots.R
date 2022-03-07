library(Seurat)
library(Signac)
library(Rsamtools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(patchwork)


colormap <- list("H3K27me3" = "#D3145A", "H3K27ac" = "#F98401", "RNAPII" = "#036C9A")

regions_me3 <- read.table("data/pbmc_bulk/mapped/me3_peaks.broadPeak", sep = "\t",
                          col.names = c("chromosome", "start", "end", "name", "score",
                                        "strand", "a", "b" , "c"))
regions_me3 <- makeGRangesFromDataFrame(regions_me3, keep.extra.columns = TRUE)
regions_me3$peak <- "H3K27me3"
regions_me3 <- keepStandardChromosomes(regions_me3, pruning.mode = "coarse")
regions_me3 <- dropSeqlevels(regions_me3, value = "chrM", pruning.mode = "coarse")

regions_ac <- read.table("data/pbmc_bulk/mapped/ac_peaks.narrowPeak", sep = "\t",
                          col.names = c("chromosome", "start", "end", "name", "score",
                                        "strand", "a", "b" , "c", "d"))
regions_ac <- makeGRangesFromDataFrame(regions_ac, keep.extra.columns = TRUE)
regions_ac$peak <- "H3K27ac"
regions_ac <- keepStandardChromosomes(regions_ac, pruning.mode = "coarse")
regions_ac <- dropSeqlevels(regions_ac, value = "chrM", pruning.mode = "coarse")

regions <- c(regions_me3, regions_ac)

QuantifyRegions <- function(fragfile, regions) {
  results <- vector(mode = "numeric", length = length(regions))
  
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
  cmp <- (results / sum(results)) * 10^6
  return(cmp)
}

# get filtered fragment files
# single-cell
frags_me3 <- "data/pbmc_protein/bigwig/pbmc_me3.bed.gz"
frags_ac <- "data/pbmc_protein/bigwig/pbmc_ac.bed.gz"

# bulk_cell
bulk_me3 <- "data/pbmc_bulk/mapped/mono/me3.bed.gz"
bulk_ac <- "data/pbmc_bulk/mapped/mono/ac.bed.gz"

# quantify
ac_sc_cov <- QuantifyRegions(frags_ac, regions)
me_sc_cov <- QuantifyRegions(frags_me3, regions)

me_bulk_cov <- QuantifyRegions(bulk_me3, regions)
ac_bulk_cov <- QuantifyRegions(bulk_ac, regions)

# create matrix
mat <- cbind(ac_sc_cov, me_sc_cov, ac_bulk_cov, me_bulk_cov)
mat <- as.data.frame(mat)
mat$region <- regions$peak

# randomize order
mat <- mat[sample(nrow(mat)), ]

get_cod_expression <- function(m) {
    eq <- substitute(~~italic(R)^2~"="~r2, 
     list(r2 = format(summary(m)$r.squared, digits = 2)))
    cod <- as.expression(eq)
    return(cod)
}

ac_me3_sc <- ggplot(mat, aes(x = ac_sc_cov, y = me_sc_cov, color = region)) +
  geom_point(size = 0.1) +
  theme_classic() +
  ylim(c(0, 600)) + xlim(c(0, 600)) +
  xlab("") +
  ylab("") +
  ggtitle(get_cod_expression(lm(ac_sc_cov ~ me_sc_cov, data = mat))) +
  scale_color_manual(values = c(colormap$H3K27ac, colormap$H3K27me3)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  theme(legend.position = "none")

ac_sc_ac_bulk <- ggplot(mat[mat$region=="H3K27ac", ], aes(x = ac_sc_cov, y = ac_bulk_cov, color = region)) +
  geom_point(size = 0.1) +
  theme_classic() +
  xlab("H3K27ac (single-cell)") +
  ylab("H3K27ac (bulk)") +
  ggtitle(get_cod_expression(lm(ac_sc_cov ~ ac_bulk_cov, data = mat))) +
  scale_color_manual(values = c(colormap$H3K27ac, colormap$H3K27me3)) +
  theme(legend.position = "none")

ggsave(filename = "plots/pbmc/scatterplot_pbmc.png", plot = ac_me3_sc, height = 3.5, width = 3.5)

## BMMC
frags_me3_bmmc <- "data/bmmc_dual/bulk_fragments/me3.bed.gz"
frags_ac_bmmc <- "data/bmmc_dual/bulk_fragments/ac.bed.gz"

# quantify
me_sc_cov <- QuantifyRegions(frags_me3_bmmc, regions)
ac_sc_cov <- QuantifyRegions(frags_ac_bmmc, regions)

# create matrix
mat_bmmc <- cbind(ac_sc_cov, me_sc_cov)
mat_bmmc <- as.data.frame(mat_bmmc)
mat_bmmc$region <- regions$peak

mat_bmmc <- mat_bmmc[sample(nrow(mat_bmmc)), ]

ac_me3_sc_bmmc <- ggplot(mat_bmmc, aes(x = ac_sc_cov, y = me_sc_cov, color = region)) +
  geom_point(size = 0.1) +
  theme_classic() +
  ylim(c(0, 600)) + xlim(c(0, 600)) +
  xlab("H3K27ac (multiplexed single-cell)") +
  scale_color_manual(values = c(colormap$H3K27ac, colormap$H3K27me3)) +
  ylab("H3K27me3 (multiplexed single-cell)") +
  ggtitle(get_cod_expression(lm(ac_sc_cov ~ me_sc_cov, data = mat_bmmc))) +
  theme(legend.position = "none") +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

ggsave(filename = "plots/bmmc/scatterplot_bmmc.png", plot = ac_me3_sc_bmmc, height = 4, width = 4)
