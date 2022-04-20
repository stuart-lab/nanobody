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
  
  frip <- total_in_peak / total
  return(frip)
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

frip_me3_mono_me3_peak <- getFrip(fragfile = me3_mono, regions = regions_me3)
frip_me3_plex_me3_peak <- getFrip(fragfile = me3_plex, regions = regions_me3)
frip_ac_mono_me3_peak <- getFrip(fragfile = ac_mono, regions = regions_me3)
frip_ac_plex_me3_peak <- getFrip(fragfile = ac_plex, regions = regions_me3)
frip_rna_mono_me3_peak <- getFrip(fragfile = rna_mono, regions = regions_me3)
frip_rna_plex_me3_peak <- getFrip(fragfile = rna_plex, regions = regions_me3)
