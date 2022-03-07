library(Signac)
library(Seurat)

atac <- readRDS("objects/bmmc_atac.rds")
bmmc_dual <- readRDS("objects/bmmc_dual.rds")

bmmc_atac <- atac[, atac$tissue == "BMMC"]

DefaultAssay(bmmc_atac) <- "ATAC"
bmmc_atac <- FindTopFeatures(bmmc_atac)
bmmc_atac <- RunTFIDF(bmmc_atac)
bmmc_atac <- RunSVD(bmmc_atac)
bmmc_atac <- RunUMAP(bmmc_atac, reduction = "lsi", dims = 2:50)

# quantify atac peaks in ntt
counts <- FeatureMatrix(
  fragments = Fragments(bmmc_dual[['ac']])[[1]],
  features = granges(bmmc_atac[['ATAC']]),
  cells = colnames(bmmc_dual)
)

bmmc_dual[['ATAC']] <- CreateChromatinAssay(
  counts = counts, fragments = Fragments(bmmc_dual[['ac']])[[1]]
)

DefaultAssay(bmmc_dual) <- "ATAC"

# normalize with reference IDF
get_idf <- function(x) {
  rsums <- rowSums(x = x)
  idf <- ncol(x = x) / rsums
  idf <- log(1 + idf)
  return(idf)
}

idf <- get_idf(x = GetAssayData(bmmc_atac, assay = "ATAC", slot = "counts"))

bmmc_dual <- RunTFIDF(bmmc_dual, idf = idf)

# find transfer anchors
anchors <- FindTransferAnchors(
  reference = bmmc_atac,
  query = bmmc_dual,
  reduction = "lsiproject",
  dims = 2:30,
  reference.reduction = "lsi",
  reference.assay = "ATAC",
  query.assay = "ATAC"
)

pred <- TransferData(
  anchorset = anchors,
  refdata = bmmc_atac$BioClassification,
  weight.reduction = bmmc_dual[['lsi.me3']],
  dims = 2:50
)

write.table(pred, "datasets/bmmc_dual_predictions.tsv", quote = FALSE, sep = "\t")