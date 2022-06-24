library(Seurat)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(tximport)
library(future)


args <- commandArgs(trailingOnly = TRUE)
threads <- as.numeric(args[[1]])
adt <- args[[2]]
ac_frag <- args[[3]]
me3_frag <- args[[4]]
outfile <- args[[5]]

plan("multicore", workers = as.numeric(threads))
options(future.globals.maxSize = Inf)

# process chromatin data
annot <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)
seqlevelsStyle(annot) <- "UCSC"

total_me3 <- CountFragments(fragments = me3_frag)
total_ac <- CountFragments(fragments = ac_frag)

total_ac <- total_ac[total_ac$frequency_count > 100, ]
total_me3 <- total_me3[total_me3$frequency_count > 100, ]
cells_keep <- intersect(total_ac$CB, total_me3$CB)

frag_obj_me3 <- CreateFragmentObject(path = me3_frag, cells = cells_keep)
frag_obj_ac <- CreateFragmentObject(path = ac_frag, cells = cells_keep)

chr_use <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[1:22]
counts_me3 <- AggregateTiles(
  object = frag_obj_me3,
  genome = chr_use,
  cells = cells_keep,
  min_counts = 1
)

counts_ac <- AggregateTiles(
  object = frag_obj_ac,
  genome = chr_use,
  cells = cells_keep,
  min_counts = 1
)

# remove features overlapping blacklist
remove_blacklist <- function(counts) {
  olap <- findOverlaps(query = StringToGRanges(rownames(counts)), subject = blacklist_hg38_unified)
  is_bl <- queryHits(olap)
  counts <- counts[setdiff(1:nrow(counts), is_bl), ]
  return(counts)
}

assay_me3 <- CreateChromatinAssay(
  counts = remove_blacklist(counts_me3),
  fragments = me3_frag,
  annotation = annot
)
assay_ac <- CreateChromatinAssay(
  counts = remove_blacklist(counts_ac),
  fragments = ac_frag,
  annotation = annot
)
obj <- CreateSeuratObject(counts = assay_me3, assay = "me3")
obj[['ac']] <- assay_ac

# add adt data
adt.txi <- tximport(files = adt, type = "alevin")
cells <- intersect(colnames(obj), colnames(adt.txi$counts))
obj <- obj[, cells]
obj[['ADT']] <- CreateAssayObject(counts = adt.txi$counts[, cells])

obj <- subset(obj, subset = nCount_me3 < 40000 &
                nCount_ac < 10000 &
                nCount_me3 > 300 &
                nCount_ac > 100 &
                nCount_ADT < 10000 &
                nCount_ADT > 120)

DefaultAssay(obj) <- "ADT"
obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2)
igg_features <- rownames(obj)[grepl("^Ig", x = rownames(obj))]
var_feat <- setdiff(rownames(obj), igg_features)
VariableFeatures(obj) <- var_feat
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:40)
obj <- FindClusters(obj)
obj <- RunUMAP(obj, reduction = 'pca', reduction.name = "umap.adt", dims = 1:40)

# remove artifact clusters
obj <- subset(obj, idents = c(4, 8), invert = TRUE)

obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:40)
obj <- FindClusters(obj)
obj <- RunUMAP(obj, reduction = 'pca', reduction.name = "umap.adt", dims = 1:40)

DefaultAssay(obj) <- "me3"
feat.keep <- names(which(rowSums(obj, slot = 'counts') > 1)) # remove low-count features
obj[['me3']] <- subset(obj[['me3']], features = feat.keep)
obj <- FindTopFeatures(obj)
obj <- RunTFIDF(obj)
obj <- RunSVD(obj, scale.embeddings = TRUE, reduction.name = 'lsi.me3')
obj <- RunUMAP(obj, reduction = 'lsi.me3', dims = 2:30, reduction.name = 'umap.me3')

DefaultAssay(obj) <- "ac"
feat.keep <- names(which(rowSums(obj, slot = 'counts') > 1)) # remove low-count features
obj[['ac']] <- subset(obj[['ac']], features = feat.keep)
obj <- FindTopFeatures(obj)
obj <- RunTFIDF(obj)
obj <- RunSVD(obj, scale.embeddings = TRUE, reduction.name = 'lsi.ac')
obj <- RunUMAP(obj, reduction = 'lsi.ac', dims = 2:30, reduction.name = 'umap.ac')

obj <- FindMultiModalNeighbors(
  object = obj,
  reduction.list = list("lsi.me3", "lsi.ac", "pca"),
  dims.list = list(2:30, 2:30, 1:30)
)
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "umap.wnn")
obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 3, resolution = 1)

Idents(obj) <- "seurat_clusters"

obj <- RenameIdents(obj, list("0" = "CD14+ Mono",
                              "1" = "CD14+ Mono",
                              "2" = "B cell",
                              "3" = "CD14+ Mono",
                              "4" = "CD4 T cell",
                              "5" = "CD16+ Mono",
                              "6" = "Late erythroid",
                              "7" = "CD8 T cell",
                              "8" = "B cell",
                              "9" = "NK",
                              "10" = "cDC"
                             )
)
obj$celltype <- Idents(obj)

obj <- TSSEnrichment(obj, fast = FALSE, assay = "me3")
obj <- TSSEnrichment(obj, fast = FALSE, assay = "ac")

saveRDS(object = obj, file = outfile)
