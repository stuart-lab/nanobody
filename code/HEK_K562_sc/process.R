library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(future)
library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
threads <- args[[1]]
outfile <- args[[2]]
fpath <- args[[3]]

plan("multicore", workers = as.integer(threads))
options(future.globals.maxSize = Inf)

# annotation
annot <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)
seqlevelsStyle(annot) <- "UCSC"
genome.use <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[1:23]

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

create_assay <- function(frag, cells, genome.use) {
  frag_obj <- CreateFragmentObject(
    path = frag,
    cells = cells
  )
  counts <- AggregateTiles(
    object = frag_obj,
    cells = cells,
    min_counts = 1,
    binsize = 10000,
    genome = genome.use
  )
  return(counts)
}

all.counts <- lapply(
  X = fragfiles,
  FUN = create_assay,
  cells = cells.keep,
  genome.use = genome.use
)

# combine matrices for same barcode pairs
names(all.counts) <- fname

all.cells <- lapply(all.counts, colnames)
common.cells <- Reduce(intersect, all.cells)
all.counts <- lapply(all.counts, function(x) x[, common.cells])
combined.matrix <- do.call(what = rbind, args = all.counts)

# create seurat object
obj <- CreateSeuratObject(
    counts = combined.matrix,
    min.cells = 10,
    assay = "all"
)

obj <- obj[, obj$nCount_all < 10000 & obj$nCount_all > 100]

for (i in seq_along(all.counts)) {
  mat <- all.counts[[i]][, colnames(obj)]
  assay <- CreateChromatinAssay(
    counts = mat,
    min.features = -1,
    fragments = fragfiles[[i]],
    annotation = annot
  )
  obj[[names(all.counts)[i]]] <- assay
}

obj <- subset(obj, subset = nCount_K27ac > 75 &
                nCount_K27me > 150 &
                nCount_Pol2 > 100)
# process
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = 50)
obj <- RunSVD(obj)
obj <- RunUMAP(obj, reduction = 'lsi', dims = 2:10)

# create lsi, umap from each assay
DefaultAssay(obj) <- "K27ac"
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = 10)
obj <- RunSVD(obj, reduction.name = "lsi.k27ac")
obj <- RunUMAP(obj, reduction = 'lsi.k27ac', reduction.name = "umap.k27ac", dims = 2:10)

DefaultAssay(obj) <- "K27me"
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = 10)
obj <- RunSVD(obj, reduction.name = "lsi.k27me")
obj <- RunUMAP(obj, reduction = 'lsi.k27me', reduction.name = "umap.k27me", dims = 2:10)

DefaultAssay(obj) <- "Pol2"
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = 10)
obj <- RunSVD(obj, reduction.name = "lsi.pol2")
obj <- RunUMAP(obj, reduction = 'lsi.pol2', reduction.name = "umap.pol2", dims = 2:10)

# WNN
obj <- FindMultiModalNeighbors(
  object = obj,
  reduction.list = list("lsi.k27ac", "lsi.k27me", "lsi.pol2"), 
  dims.list = list(2:10, 2:10, 2:10)
)

obj <- RunUMAP(
  object = obj,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap"
)

obj <- FindClusters(
  object = obj,
  graph.name = "wsnn",
  algorithm = 3,
  resolution = 0.05
)

obj <- RenameIdents(obj, list("0" = "K562", "1" = "HEK"))
obj$celltype <- Idents(obj)

# save object
saveRDS(object = obj, file = outfile)

# write cell names for vireo
writeLines(colnames(obj), "data/HEK_K562_sc/cells.txt")
                     
# split fragment files by cell type
hek <- WhichCells(object = obj, idents = "HEK")
k562 <- WhichCells(object = obj, idents = "K562")

filter_chrom <- function(infile, outfile, chrom.keep) {
  # open input and output files
  con <- gzfile(infile, 'rb')
  outf <- file(description = outfile, "w")

  # read line by line, writing to output file
  while (TRUE) {
    line = readLines(con, n = 1)
    if (length(line) == 0) {
      break
    } else if (substring(text = line, first = 1, last = 1) == "#") {
      next # skip header
    } else {
      # check chromosome
      fields <- unlist(x = strsplit(x = line, "\t"))
      if (!(fields[[1]] %in% chrom.keep)) {
        next
      } else {
        writeLines(text = line, con = outf)
      }
    }
  }
  close(con)
  close(outf)
  system(
    command = paste0("bgzip ", outfile),
    wait = TRUE,
    ignore.stderr = FALSE,
    ignore.stdout = FALSE
  )
  system(
    command = paste0("tabix -p bed ", paste0(outfile, '.gz')),
    wait = TRUE,
    ignore.stderr = FALSE,
    ignore.stdout = FALSE
  )
}

chr.use <- names(genome.use)
for (i in seq_along(fragfiles)) {
  outf <- paste0("data/HEK_K562_sc/split/",
                 names(fragfiles)[[i]])
  FilterCells(
    fragments = fragfiles[[i]],
    cells = k562,
    outfile = paste0(outf, "_k562.bed.gz")
  )
  FilterCells(
    fragments = fragfiles[[i]],
    cells = hek,
    outfile = paste0(outf, "_hek.bed.gz")
  )

  # filter chromosomes
  filter_chrom(
    infile = paste0(outf, "_k562.bed.gz"),
    outfile = paste0(outf, "_k562_filter.bed"),
    chrom.keep = chr.use
  )
  filter_chrom(
    infile = paste0(outf, "_hek.bed.gz"),
    outfile = paste0(outf, "_hek_filter.bed"),
    chrom.keep = chr.use
  )
}
