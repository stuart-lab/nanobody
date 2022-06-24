library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(future)
library(GenomicRanges)


args <- commandArgs(trailingOnly = TRUE)
threads <- args[[1]]
outfile <- args[[2]]

plan("multicore", workers = as.integer(threads))
options(future.globals.maxSize = Inf)

# annotation
annot <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)
seqlevelsStyle(annot) <- "UCSC"
genome.use <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[1:23]

fpath <- "data/bmmc_dual/DNA"
fragfiles <- list.files(
  path = fpath,
  pattern = "*.bed.gz$",
  full.names = TRUE
)
fname <- gsub(pattern = paste0(fpath, "/"),
              replacement = "",
              x = fragfiles)
fname <- gsub(pattern = ".bed.gz",
              replacement = "",
              x = fname)
names(fragfiles) <- fname

total_counts <- CountFragments(fragments = as.list(fragfiles))
cells.keep <- total_counts[total_counts$frequency_count > 150, "CB"]

frags_me3 <- CreateFragmentObject(
  path = fragfiles[['H3K27me3']],
  cells = cells.keep
)
frags_ac <- CreateFragmentObject(
  path = fragfiles[['H3K27ac']],
  cells = cells.keep
)

counts_me3 <- AggregateTiles(
  object = frags_me3,
  cells = cells.keep,
  min_counts = 1,
  binsize = 5000,
  genome = genome.use
)

counts_ac <- AggregateTiles(
  object = frags_ac,
  cells = cells.keep,
  min_counts = 3,
  binsize = 1000,
  genome = genome.use
)

# remove features overlapping blacklist
remove_blacklist <- function(counts) {
  olap <- findOverlaps(query = StringToGRanges(rownames(counts)), subject = blacklist_hg38_unified)
  is_bl <- queryHits(olap)
  counts <- counts[setdiff(1:nrow(counts), is_bl), ]
  return(counts)
}

# create seurat object
me3 <- CreateChromatinAssay(
  counts = remove_blacklist(counts_me3),
  fragments = fragfiles[['H3K27me3']],
  annotation = annot
)

ac <- CreateChromatinAssay(
  counts = remove_blacklist(counts_ac),
  fragments = fragfiles[['H3K27ac']],
  annotation = annot
)

bmmc <- CreateSeuratObject(counts = me3, assay = "me3")
bmmc[['ac']] <- ac

bmmc <- subset(bmmc, subset = nCount_me3 < 10000 &
                 nCount_ac < 10000 &
                 nCount_me3 > 100 &
                 nCount_ac > 75)

bmmc <- TSSEnrichment(bmmc, assay = "me3")
bmmc$tss.me3 <- bmmc$TSS.enrichment
bmmc <- TSSEnrichment(bmmc, assay = "ac")
bmmc$tss.ac <- bmmc$TSS.enrichment

bmmc <- NucleosomeSignal(bmmc, assay = "me3")
bmmc$nucleosome.me3 <- bmmc$nucleosome_signal
bmmc <- NucleosomeSignal(bmmc, assay = "ac")
bmmc$nucleosome.ac <- bmmc$nucleosome_signal

bmmc <- subset(bmmc, subset = tss.ac > 1)

DefaultAssay(bmmc) <- "me3"
feat.keep <- names(which(rowSums(bmmc, slot = 'counts') > 1)) # remove low-count features
bmmc[['me3']] <- subset(bmmc[['me3']], features = feat.keep)
bmmc <- FindTopFeatures(bmmc)
bmmc <- RunTFIDF(bmmc, scale.factor = median(bmmc$nCount_me3))
bmmc <- RunSVD(bmmc, reduction.name = "lsi.me3", n = 100, scale.embeddings = FALSE)
bmmc <- RunUMAP(bmmc, reduction = "lsi.me3", reduction.name = "umap.me3", dims = 2:100)

DefaultAssay(bmmc) <- "ac"
feat.keep <- names(which(rowSums(bmmc, slot = 'counts') > 1))
bmmc[['ac']] <- subset(bmmc[['ac']], features = feat.keep)
bmmc <- FindTopFeatures(bmmc)
bmmc <- RunTFIDF(bmmc, scale.factor = median(bmmc$nCount_ac))
bmmc <- RunSVD(bmmc, reduction.name = "lsi.ac", n = 100, scale.embeddings = FALSE)
bmmc <- RunUMAP(bmmc, reduction = "lsi.ac", reduction.name = "umap.ac", dims = 2:60)

bmmc <- FindMultiModalNeighbors(
  object = bmmc,
  reduction.list = list("lsi.me3", "lsi.ac"),
  dims.list = list(2:100, 2:60)
)
bmmc <- RunUMAP(bmmc, nn.name = "weighted.nn", reduction.name = "umap.wnn")
bmmc <- FindClusters(bmmc, graph.name = "wsnn", algorithm = 3, resolution = 4)

renaming <- list(
  '0' = "NK",
  '1' = "NK",
  '2' = "GMP/CLP",
  '3' = 'CD8 Naive',
  '4' = 'CD4 Naive',
  '5' = 'CD14 Mono',
  '6' = 'pDC',
  '7' = 'B',
  '8' = 'CD4 Memory',
  '9' = 'CD8 Memory',
  '10' = 'CD14 Mono',
  '11' = 'Pre-B',
  '12' = 'CD8 Memory',
  '13' = 'Late erythroid',
  '14' = 'Early erythroid',
  '15' = 'CD14 Mono',
  '16' = 'B',
  '17' = 'HSPC',
  '18' = 'CD14 Mono',
  '19' = 'CD14 Mono',
  '20' = 'Plasma',
  '21' = 'Early basophil'
)
Idents(bmmc) <- "seurat_clusters"
bmmc <- RenameIdents(bmmc, renaming)
bmmc$celltype <- Idents(bmmc)

saveRDS(object = bmmc, file = "objects/bmmc_dual.rds")
