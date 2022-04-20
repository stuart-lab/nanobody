library(GenomicRanges)

regions_s5 <- read.table("data/k562_peaks/ENCFF053XYZ.bed.gz", sep = "\t",
                         col.names = c("chromosome", "start", "end", "name", "score",
                                       "strand", "a", "b" , "c", "d"))
regions_s5 <- makeGRangesFromDataFrame(regions_s5, keep.extra.columns = TRUE)

regions_s2 <- read.table("data/k562_peaks/ENCFF266OPF.bed.gz", sep = "\t",
                         col.names = c("chromosome", "start", "end", "name", "score",
                                       "strand", "a", "b" , "c", "d"))
regions_s2 <- makeGRangesFromDataFrame(regions_s2, keep.extra.columns = TRUE)

regions_rna <- reduce(c(regions_s2, regions_s5))

regions_rna <- as.data.frame(regions_rna)
write.table(regions_rna, file = "data/k562_peaks/rna.bed", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
