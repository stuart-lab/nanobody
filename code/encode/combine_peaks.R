library(GenomicRanges)

me3_peaks <- c("ENCFF138FVQ", "ENCFF720AKK",
               "ENCFF407VKP","ENCFF588JEN",
               "ENCFF570DGJ", "ENCFF071JFL",
               "ENCFF276IMO")

ac_peaks <- c("ENCFF240LSH", "ENCFF614FJO",
              "ENCFF449ITG", "ENCFF321NYQ",
              "ENCFF544QYY", "ENCFF059WXH",
              "ENCFF412NPA")

bulk_me3 <- "data/encode/ENCFF291LVP.bed.gz"
bulk_ac <- "data/encode/ENCFF832RWT.bed.gz"

me3_bedfiles <- paste0("data/encode/", me3_peaks, ".bed.gz")
ac_bedfiles <- paste0("data/encode/", ac_peaks, ".bed.gz")

read_bed <- function(filepath) {
    df <- read.table(file = filepath,
                     col.names = c("chr", "start", "end",
                                   "name", "score", "strand",
                                   "fold_change", "neg_log10pvalue_summit",
                                   "neg_log10qvalue_summit", "relative_summit_position"))
    gr <- makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE, 
        starts.in.df.are.0based = TRUE)
    return(gr)
}

me3_ranges <- sapply(me3_bedfiles, read_bed)
ac_ranges <- sapply(ac_bedfiles, read_bed)

all_me3 <- Reduce(f = c, x = me3_ranges)
me3 <- reduce(all_me3)
me3 <- keepStandardChromosomes(me3, pruning.mode = "coarse")
me3 <- me3[width(me3) < 20000]

all_ac <- Reduce(f = c, x = ac_ranges)
ac <- reduce(all_ac)
ac <- keepStandardChromosomes(ac, pruning.mode = "coarse")
ac <- ac[width(ac) < 20000]

combined <- c(ac, me3)

ac <- as.data.frame(ac)
write.table(x = ac, file = "data/encode/h3k27ac.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
me3 <- as.data.frame(me3)
write.table(x = me3, file = "data/encode/h3k27me3.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
combined <- as.data.frame(combined)
write.table(x = combined, file = "data/encode/all.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# bulk peaks
ac_bulk <- read_bed(bulk_ac)
me3_bulk <- read_bed(bulk_me3)
all_bulk <- c(ac_bulk, me3_bulk)
all_bulk <- as.data.frame(all_bulk)
write.table(x = all_bulk, file = "data/encode/all_bulk.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
