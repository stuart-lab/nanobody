library(Signac)
library(Seurat)

obj <- readRDS("objects/bmmc_dual.rds")
outdir <- "data/bmmc_dual/bulk_fragments/"

# remove non-cell fragments for each assay for genome bin quant
DefaultAssay(obj) <- "ac"
frags <- Fragments(obj)[[1]]@path
FilterCells(
  fragments = frags,
  cells = colnames(obj),
  outfile = paste0(outdir, "ac.bed.gz"),
  verbose = TRUE
)

DefaultAssay(obj) <- "me3"
frags <- Fragments(obj)[[1]]@path
FilterCells(
  fragments = frags,
  cells = colnames(obj),
  outfile = paste0(outdir, "me3.bed.gz"),
  verbose = TRUE
)
