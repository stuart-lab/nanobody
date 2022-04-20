library(ggplot2)
library(GenomicRanges)
library(Rsamtools)

getFrip <- function(fragfile, regions) {
  results <- vector(mode = "numeric", length = length(regions))
  
  # open tabix connection
  tabix.file <- TabixFile(file = fragfile)
  open(con = tabix.file)
  
  in_peak <- countTabix(tabix.file, param = regions)
  total_in_peak <- sum(unlist(in_peak))
  
  # close connection
  close(con = tabix.file)
  
  total <- system(
    command = paste0("gzip -dc ", fragfile, " | wc -l"),
    wait = TRUE,
    intern = TRUE
  )
  total <- as.numeric(total)
  
  frip <- total_in_peak / total
  return(frip)
}


# bulk fragments
me3_mono <- "data/pbmc_bulk/mapped/mono/me3.bed.gz"
ac_mono <- "data/pbmc_bulk/mapped/mono/ac.bed.gz"
me3_plex <- "data/pbmc_bulk/mapped/plex/me3.bed.gz"
ac_plex <- "data/pbmc_bulk/mapped/plex/ac.bed.gz"

regions_ac <- read.table("data/encode/ENCFF832RWT.bed.gz", sep = "\t",
                         col.names = c("chromosome", "start", "end", "name", "score",
                                       "strand", "a", "b" , "c", "d"))
regions_ac <- makeGRangesFromDataFrame(regions_ac, keep.extra.columns = TRUE)

regions_me3 <- read.table("data/encode/ENCFF291LVP.bed.gz", sep = "\t",
                          col.names = c("chromosome", "start", "end", "name", "score",
                                        "strand", "a", "b" , "c", "d"))
regions_me3 <- makeGRangesFromDataFrame(regions_me3, keep.extra.columns = TRUE)

frip_me3_mono_me3_peak <- getFrip(fragfile = me3_mono, regions = regions_me3)
frip_me3_plex_me3_peak <- getFrip(fragfile = me3_plex, regions = regions_me3)
frip_ac_mono_me3_peak <- getFrip(fragfile = ac_mono, regions = regions_me3)
frip_ac_plex_me3_peak <- getFrip(fragfile = ac_plex, regions = regions_me3)

frip_me3_mono_ac_peak <- getFrip(fragfile = me3_mono, regions = regions_ac)
frip_me3_plex_ac_peak <- getFrip(fragfile = me3_plex, regions = regions_ac)
frip_ac_mono_ac_peak <- getFrip(fragfile = ac_mono, regions = regions_ac)
frip_ac_plex_ac_peak <- getFrip(fragfile = ac_plex, regions = regions_ac)

me3_peak <- data.frame(
  frip = c(frip_me3_mono_me3_peak, frip_me3_plex_me3_peak, frip_ac_mono_me3_peak, frip_ac_plex_me3_peak),
  assay = c("H3K27me3_mono", "H3K27me3_plex", "H3K27ac_mono", "H3K27ac_plex")
)
me3_peak$assay <- factor(me3_peak$assay, levels = c("H3K27me3_plex", "H3K27ac_plex", "H3K27me3_mono", "H3K27ac_mono"))

p1 <- ggplot(me3_peak, aes(assay, frip, fill = assay)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#D3145A", "#F98401", "#676767", "#676767")) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("FRiP") +
  ggtitle("Fraction of reads in H3K27me3 peaks")

ac_peak <- data.frame(
  frip = c(frip_me3_mono_ac_peak, frip_me3_plex_ac_peak, frip_ac_mono_ac_peak, frip_ac_plex_ac_peak),
  assay = c("H3K27me3_mono", "H3K27me3_plex", "H3K27ac_mono", "H3K27ac_plex")
)
ac_peak$assay <- factor(ac_peak$assay, levels = c("H3K27me3_plex", "H3K27ac_plex", "H3K27me3_mono", "H3K27ac_mono"))

p2 <- ggplot(ac_peak, aes(assay, frip, fill = assay)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#D3145A", "#F98401", "#676767", "#676767")) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("FRiP") +
  ggtitle("Fraction of reads in H3K27ac peaks")

ggsave(filename = "plots/pbmc_bulk/frip_me3.png", plot = p1, height = 4, width = 4)
ggsave(filename = "plots/pbmc_bulk/frip_ac.png", plot = p2, height = 4, width = 4)

