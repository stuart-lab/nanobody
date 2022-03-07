library(Signac)
library(Seurat)

obj <- readRDS("objects/pbmc_protein.rds")

outdir <- "data/pbmc_protein/bigwig/"
celltypes.split <- c("CD14+ Mono", "B cell", "NK", "CD8+ T cell", "CD4+ T cell", "HSC")

Idents(obj) <- "celltype"

# split fragments for each assay
DefaultAssay(obj) <- "ac"
frags <- Fragments(obj)[[1]]@path
for (i in celltypes.split) {
    ct_safe <- gsub(" ", "_", i, fixed = TRUE)
    ct_safe <- gsub("+", "", ct_safe, fixed = TRUE)
    message(ct_safe)
    FilterCells(
      fragments = frags,
      cells = WhichCells(object = obj, idents = i),
      outfile = paste0(outdir, ct_safe, "_ac.bed.gz"),
      verbose = TRUE
    )
}

DefaultAssay(obj) <- "me3"
frags <- Fragments(obj)[[1]]@path
for (i in celltypes.split) {
    ct_safe <- gsub(" ", "_", i, fixed = TRUE)
    ct_safe <- gsub("+", "", ct_safe, fixed = TRUE)
    message(ct_safe)
    FilterCells(
      fragments = frags,
      cells = WhichCells(object = obj, idents = i),
      outfile = paste0(outdir, ct_safe, "_me3.bed.gz"),
      verbose = TRUE
    )
}

# remove non-cell fragments for each assay for genome bin quant
DefaultAssay(obj) <- "ac"
frags <- Fragments(obj)[[1]]@path
FilterCells(
  fragments = frags,
  cells = colnames(obj),
  outfile = paste0(outdir, "pbmc_ac.bed.gz"),
  verbose = TRUE
)

DefaultAssay(obj) <- "me3"
frags <- Fragments(obj)[[1]]@path
FilterCells(
  fragments = frags,
  cells = colnames(obj),
  outfile = paste0(outdir, "pbmc_me3.bed.gz"),
  verbose = TRUE
)
