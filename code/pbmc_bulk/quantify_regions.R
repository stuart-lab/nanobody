library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(Rsamtools)


me3_mono <- "data/pbmc_bulk/mapped/mono/me3.bed.gz"
ac_mono <- "data/pbmc_bulk/mapped/mono/ac.bed.gz"
me3_plex <- "data/pbmc_bulk/mapped/plex/me3.bed.gz"
ac_plex <- "data/pbmc_bulk/mapped/plex/ac.bed.gz"

regions_ac <- read.table("data/encode/ENCFF832RWT.bed.gz", sep = "\t",
                         col.names = c("chromosome", "start", "end", "name", "score",
                                       "strand", "a", "b" , "c", "d"))
regions_ac <- makeGRangesFromDataFrame(regions_ac, keep.extra.columns = TRUE)
regions_ac$peak <- "H3K27ac"

regions_me3 <- read.table("data/encode/ENCFF291LVP.bed.gz", sep = "\t",
                         col.names = c("chromosome", "start", "end", "name", "score",
                                       "strand", "a", "b" , "c", "d"))
regions_me3 <- makeGRangesFromDataFrame(regions_me3, keep.extra.columns = TRUE)
regions_me3$peak <- "H3K27me3"

regions <- c(regions_me3, regions_ac)

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

ac_plex_cov <- QuantifyRegions(ac_plex, regions)
ac_mono_cov <- QuantifyRegions(ac_mono, regions)
me3_plex_cov <- QuantifyRegions(me3_plex, regions)
me3_mono_cov <- QuantifyRegions(me3_mono, regions)

mat <- cbind(ac_plex_cov, ac_mono_cov, me3_plex_cov, me3_mono_cov)

mat <- as.data.frame(mat)
mat$region <- regions$peak
colormap <- list("H3K27me3" = "#D3145A", "H3K27ac" = "#F98401", "RNAPII" = "#036C9A")

mat <- mat[sample(nrow(mat)), ]

# p-m
m.use <- mat[mat$region == "H3K27ac", ]
ac_plex_mono <- ggplot(
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

m.use <- mat[mat$region == "H3K27me3", ]
me3_plex_mono <- ggplot(
  data = m.use,
  mapping = aes(me3_plex_cov, me3_mono_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 100)) + ylim(c(0, 100)) +
  xlab("H3K27me3 (multiple antibody)") +
  ylab("H3K27me3 (single antibody)") +
  ggtitle(get_cod_expression(lm(me3_plex_cov ~ me3_mono_cov, data = m.use))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

ac_plex_me3_mono <- ggplot(
  data = mat,
  mapping = aes(ac_plex_cov, me3_mono_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 400)) + ylim(c(0, 400)) +
  xlab("H3K27ac (multiple antibody)") +
  ylab("H3K27me3 (single antibody)") +
  ggtitle(get_cod_expression(lm(ac_plex_cov ~ me3_mono_cov, data = mat))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

ac_mono_me3_plex <- ggplot(
  data = mat,
  mapping = aes(ac_mono_cov, me3_plex_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 400)) + ylim(c(0, 400)) +
  xlab("H3K27ac (single antibody)") +
  ylab("H3K27me3 (multiple antibody)") +
  ggtitle(get_cod_expression(lm(ac_mono_cov ~ me3_plex_cov, data = mat))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

# m-m
ac_mono_me3_mono <- ggplot(
  data = mat,
  mapping = aes(ac_mono_cov, me3_mono_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 400)) + ylim(c(0, 400)) +
  xlab("H3K27ac (single antibody)") +
  ylab("H3K27me3 (single antibody)") +
  ggtitle(get_cod_expression(lm(ac_mono_cov ~ me3_mono_cov, data = mat))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

# p-p
ac_plex_me3_plex <- ggplot(
  data = mat,
  mapping = aes(ac_plex_cov, me3_plex_cov, color = region)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = colormap) +
  theme_classic() +
  xlim(c(0, 400)) + ylim(c(0, 400)) +
  xlab("H3K27ac (multiple antibody)") +
  ylab("H3K27me3 (multiple antibody)") +
  ggtitle(get_cod_expression(lm(ac_plex_cov ~ me3_plex_cov, data = mat))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2)

pp <- (ac_mono_me3_mono | ac_mono_me3_plex | ac_plex_me3_mono | ac_plex_me3_plex | ac_plex_mono | me3_plex_mono) & theme(legend.position = "none")
ggsave(filename = "plots/pbmc/bulk_scatter.png", plot = pp, height = 3.5, width = 16, dpi = 600)
