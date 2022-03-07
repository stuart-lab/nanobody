library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(future)


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
  min_counts = 1,
  binsize = 5000,
  genome = genome.use
)

# combine matrices for same barcode pairs
common.cells <- intersect(colnames(counts_me3), colnames(counts_ac))
all.counts <- rbind(counts_me3[, common.cells], counts_ac[, common.cells])

# create seurat object
bmmc <- CreateSeuratObject(
  counts = all.counts,
  min.cells = 2,
  assay = "all"
)
bmmc[['me3']] <- CreateChromatinAssay(
  counts = counts_me3[, colnames(bmmc)],
  fragments = fragfiles[['H3K27me3']],
  annotation = annot,
  min.features = -1
)
bmmc[['ac']] <- CreateChromatinAssay(
  counts = counts_ac[, colnames(bmmc)],
  fragments = fragfiles[['H3K27ac']],
  annotation = annot,
  min.features = -1
)

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

DefaultAssay(bmmc) <- "me3"
bmmc <- FindTopFeatures(bmmc)
bmmc <- RunTFIDF(bmmc, scale.factor = median(bmmc$nCount_me3))
bmmc <- RunSVD(bmmc, reduction.name = "lsi.me3", scale.embeddings = FALSE)
bmmc <- RunUMAP(bmmc, reduction = "lsi.me3", reduction.name = "umap.me3", dims = 2:50)

DefaultAssay(bmmc) <- "ac"
bmmc <- FindTopFeatures(bmmc)
bmmc <- RunTFIDF(bmmc, scale.factor = median(bmmc$nCount_ac))
bmmc <- RunSVD(bmmc, reduction.name = "lsi.ac", n = 100, scale.embeddings = FALSE)
bmmc <- RunUMAP(bmmc, reduction = "lsi.ac", reduction.name = "umap.ac", dims = 2:80)

bmmc <- FindMultiModalNeighbors(
  object = bmmc,
  reduction.list = list("lsi.me3", "lsi.ac"),
  dims.list = list(2:50, 2:80)
)
bmmc <- RunUMAP(bmmc, nn.name = "weighted.nn", reduction.name = "umap.wnn")
bmmc <- FindClusters(bmmc, graph.name = "wsnn", algorithm = 3, resolution = 3)

renaming <- list(
  '0' = "CD14 Mono",
  '1' = "GMP/CLP",
  '2' = "NK",
  '3' = 'NK',
  '4' = 'CD8 Memory',
  '5' = 'CD4 Naive',
  '6' = 'CD14 Mono',
  '7' = 'pDC',
  '8' = 'CD4 Memory',
  '9' = 'B',
  '10' = 'CD8 Naive',
  '11' = 'CD8 Memory',
  '12' = 'Late erythroid',
  '13' = 'Pre-B',
  '14' = 'Early erythroid',
  '15' = 'CD14 Mono',
  '16' = 'HSPC',
  '17' = 'B',
  '18' = 'Plasma',
  '19' = 'CD14 Mono',
  '20' = 'Early basophil'
)
Idents(bmmc) <- "seurat_clusters"
bmmc <- RenameIdents(bmmc, renaming)
bmmc$celltype <- Idents(bmmc)

saveRDS(object = bmmc, file = "objects/bmmc_dual.rds")
