library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(Rsamtools)

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

get_cod_expression <- function(m) {
  eq <- substitute(~~italic(R)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 2)))
  cod <- as.expression(eq)
  return(cod)
}

# bulk fragments
me3_mono <- "data/HEK_K562_bulk/mapped/K562-mono-K27me.bed.gz"
me3_plex <- "data/HEK_K562_bulk/mapped/K562-plex-K27me.bed.gz"
ac_mono <- "data/HEK_K562_bulk/mapped/K562-mono-K27ac.bed.gz"
ac_plex <- "data/HEK_K562_bulk/mapped/K562-plex-K27ac.bed.gz"
rna_mono <- "data/HEK_K562_bulk/mapped/K562-mono-Pol2.bed.gz"
rna_plex <- "data/HEK_K562_bulk/mapped/K562-plex-Pol2.bed.gz"

# load peaks
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

# quantify
ac_mono_cov <- QuantifyRegions(ac_mono, regions)
ac_plex_cov <- QuantifyRegions(ac_plex, regions)
me3_mono_cov <- QuantifyRegions(me3_mono, regions)
me3_plex_cov <- QuantifyRegions(me3_plex, regions)
rna_mono_cov <- QuantifyRegions(rna_mono, regions)
rna_plex_cov <- QuantifyRegions(rna_plex, regions)

mat <- cbind(ac_mono_cov, ac_plex_cov, me3_mono_cov, me3_plex_cov, rna_mono_cov, rna_plex_cov)

mat <- as.data.frame(mat)
mat$region <- regions$peak
colormap <- list("H3K27me3" = "#D3145A", "H3K27ac" = "#F98401", "RNAPII" = "#036C9A")

# p-m combinations

## AC
m.use <- mat[mat$region %in% c("H3K27ac", "H3K27me3"), ]
ac_plex_me3_mono <- ggplot(
  data = m.use,
  mapping = aes(ac_plex_cov, me3_mono_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 500)) + ylim(c(0, 500)) +
  xlab("H3K27ac (multiple antibody)") +
  ylab("H3K27me3 (single antibody)") +
  ggtitle(get_cod_expression(lm(ac_plex_cov ~ me3_mono_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

m.use <- mat[mat$region %in% c("H3K27ac", "RNAPII"), ]
ac_plex_rna_mono <- ggplot(
  data = m.use,
  mapping = aes(ac_plex_cov, rna_mono_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 500)) + ylim(c(0, 500)) +
  xlab("H3K27ac (multiple antibody)") +
  ylab("RNAPII (single antibody)") +
  ggtitle(get_cod_expression(lm(ac_plex_cov ~ rna_mono_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

m.use <- mat[mat$region == "H3K27ac", ]
ac_plex_ac_mono <- ggplot(
  data = m.use,
  mapping = aes(ac_plex_cov, ac_mono_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 800)) + ylim(c(0, 800)) +
  xlab("H3K27ac (multiple antibody)") +
  ylab("H3K27ac (single antibody)") +
  ggtitle(get_cod_expression(lm(ac_plex_cov ~ ac_mono_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

## ME3
m.use <- mat[mat$region %in% c("H3K27me3", "RNAPII"), ]
me3_plex_rna_mono <- ggplot(
  data = m.use,
  mapping = aes(me3_plex_cov, rna_mono_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 600)) + ylim(c(0, 600)) +
  xlab("H3K27me3 (multiple antibody)") +
  ylab("RNAPII (single antibody)") +
  ggtitle(get_cod_expression(lm(me3_plex_cov ~ rna_mono_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

m.use <- mat[mat$region %in% c("H3K27me3", "H3K27ac"), ]
me3_plex_ac_mono <- ggplot(
  data = m.use,
  mapping = aes(me3_plex_cov, ac_mono_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 800)) + ylim(c(0, 800)) +
  xlab("H3K27me3 (multiple antibody)") +
  ylab("H2K27ac (single antibody)") +
  ggtitle(get_cod_expression(lm(me3_plex_cov ~ ac_mono_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

m.use <- mat[mat$region == "H3K27me3", ]
me3_plex_me3_mono <- ggplot(
  data = m.use,
  mapping = aes(me3_plex_cov, me3_mono_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 1000)) + ylim(c(0, 1000)) +
  xlab("H3K27me3 (multiple antibody)") +
  ylab("H3K27me3 (single antibody)") +
  ggtitle(get_cod_expression(lm(me3_plex_cov ~ me3_mono_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

## RNA
m.use <- mat[mat$region %in% c("RNAPII", "H3K27ac"), ]
rna_plex_ac_mono <- ggplot(
  data = m.use,
  mapping = aes(rna_plex_cov, ac_mono_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 800)) + ylim(c(0, 800)) +
  xlab("RNAPII (multiple antibody)") +
  ylab("H2K27ac (single antibody)") +
  ggtitle(get_cod_expression(lm(rna_plex_cov ~ ac_mono_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

m.use <- mat[mat$region %in% c("RNAPII", "H3K27me3"), ]
rna_plex_me3_mono <- ggplot(
  data = m.use,
  mapping = aes(rna_plex_cov, me3_mono_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 500)) + ylim(c(0, 500)) +
  xlab("RNAPII (multiple antibody)") +
  ylab("H2K27me3 (single antibody)") +
  ggtitle(get_cod_expression(lm(rna_plex_cov ~ me3_mono_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

m.use <- mat[mat$region == "RNAPII", ]
rna_plex_rna_mono <- ggplot(
  data = m.use,
  mapping = aes(rna_plex_cov, rna_mono_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 600)) + ylim(c(0, 600)) +
  xlab("RNAPII (multiple antibody)") +
  ylab("RNAPII (single antibody)") +
  ggtitle(get_cod_expression(lm(rna_plex_cov ~ rna_mono_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

# m-m combinations
## AC
m.use <- mat[mat$region %in% c("H3K27ac", "H3K27me3"), ]
ac_mono_me3_mono <- ggplot(
  data = m.use,
  mapping = aes(ac_mono_cov, me3_mono_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 800)) + ylim(c(0, 800)) +
  xlab("H3K27ac (single antibody)") +
  ylab("H3K27me3 (single antibody)") +
  ggtitle(get_cod_expression(lm(ac_mono_cov ~ me3_mono_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

m.use <- mat[mat$region %in% c("H3K27ac", "RNAPII"), ]
ac_mono_rna_mono <- ggplot(
  data = m.use,
  mapping = aes(ac_mono_cov, rna_mono_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 800)) + ylim(c(0, 800)) +
  xlab("H3K27ac (single antibody)") +
  ylab("RNAPII (single antibody)") +
  ggtitle(get_cod_expression(lm(ac_mono_cov ~ rna_mono_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

## ME3
m.use <- mat[mat$region %in% c("H3K27me3", "RNAPII"), ]
me3_mono_rna_mono <- ggplot(
  data = m.use,
  mapping = aes(me3_mono_cov, rna_mono_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 800)) + ylim(c(0, 800)) +
  xlab("H3K27me3 (single antibody)") +
  ylab("RNAPII (single antibody)") +
  ggtitle(get_cod_expression(lm(me3_mono_cov ~ rna_mono_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

# p-p
m.use <- mat[mat$region %in% c("H3K27ac", "H3K27me3"), ]
ac_plex_me3_plex <-  ggplot(
  data = m.use,
  mapping = aes(ac_plex_cov, me3_plex_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 600)) + ylim(c(0, 600)) +
  xlab("H3K27ac (multiple antibody)") +
  ylab("H3K27me3 (multiple antibody)") +
  ggtitle(get_cod_expression(lm(ac_plex_cov ~ me3_plex_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

m.use <- mat[mat$region %in% c("H3K27ac", "RNAPII"), ]
ac_plex_rna_plex <-  ggplot(
  data = m.use,
  mapping = aes(ac_plex_cov, rna_plex_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 500)) + ylim(c(0, 500)) +
  xlab("H3K27ac (multiple antibody)") +
  ylab("RNAPII (multiple antibody)") +
  ggtitle(get_cod_expression(lm(ac_plex_cov ~ rna_plex_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

m.use <- mat[mat$region %in% c("H3K27me3", "RNAPII"), ]
me3_plex_rna_plex <-  ggplot(
  data = m.use,
  mapping = aes(me3_plex_cov, rna_plex_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 600)) + ylim(c(0, 600)) +
  xlab("H3K27me3 (multiple antibody)") +
  ylab("RNAPII (multiple antibody)") +
  ggtitle(get_cod_expression(lm(me3_plex_cov ~ rna_plex_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

pp <- (ac_mono_me3_mono / ac_mono_rna_mono / me3_mono_rna_mono) | (ac_plex_me3_mono / ac_plex_rna_mono / me3_plex_rna_mono) | (me3_plex_ac_mono / rna_plex_ac_mono / rna_plex_me3_mono) | (ac_plex_me3_plex / ac_plex_rna_plex / me3_plex_rna_plex) | (me3_plex_me3_mono / ac_plex_ac_mono / rna_plex_rna_mono)
pp <- pp & theme(legend.position = "none")
ggsave(filename = "plots/hek_k562/scatterplots_bulk.png", plot = pp, height = 12, width = 17)
