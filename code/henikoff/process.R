library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)

annot <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)
seqlevelsStyle(annot) <- "UCSC"

me3 <- list("data/henikoff/GSM5034342_K27me3_R1_PBMC.fragments.HG38.tsv.gz",
            "data/henikoff/GSM5034343_K27me3_R2_PBMC.fragments.HG38.tsv.gz")
ac <- "data/henikoff/GSM5034344_K27ac_PBMC.fragments.HG38.tsv.gz"

# count fragments per cell
frags_me3_1 <- CountFragments(me3[[1]])
frags_me3_2 <- CountFragments(me3[[2]])
frags_ac <- CountFragments(ac)

min_frags <- 500
cells_me3_1 <- frags_me3_1[frags_me3_1$frequency_count > min_frags, "CB"]
cells_me3_2 <- frags_me3_2[frags_me3_2$frequency_count > min_frags, "CB"]
cells_ac <- frags_ac[frags_ac$frequency_count > min_frags, "CB"]

f_me3_1 <- CreateFragmentObject(me3[[1]], cells = cells_me3_1)
f_me3_2 <- CreateFragmentObject(me3[[2]], cells = cells_me3_2)
f_ac <- CreateFragmentObject(ac, cells = cells_ac)

me3_counts_1 <- FeatureMatrix(
  fragments = f_me3_1,
  features = granges(pbmc[['me3']]),
  cells = cells_me3_1
)
me3_counts_2 <- FeatureMatrix(
  fragments = f_me3_2,
  features = granges(pbmc[['me3']]),
  cells = cells_me3_2
)
ac_counts <- FeatureMatrix(
  fragments = f_ac,
  features = granges(pbmc[['ac']]),
  cells = cells_ac
)

# change cell barcodes to add batch info
colnames(me3_counts_1) <- paste0(colnames(me3_counts_1), ":1")
colnames(me3_counts_2) <- paste0(colnames(me3_counts_2), ":2")
names(f_me3_1@cells) <- paste0(names(f_me3_1@cells), ":1")
names(f_me3_2@cells) <- paste0(names(f_me3_2@cells), ":2")

me3 <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = cbind(me3_counts_1, me3_counts_2),
    fragments = list(f_me3_1, f_me3_2),
    annotation = annot
  ),
  names.field = 2,
  names.delim = ":",
  assay = "me3",
  meta.data = md
)

me3 <- FindTopFeatures(me3)
me3 <- RunTFIDF(me3)
me3 <- RunSVD(me3, scale.embeddings = TRUE, reduction.name = 'lsi.me3')
me3 <- RunUMAP(me3, reduction = 'lsi.me3', dims = 2:30, reduction.name = 'umap.me3')

saveRDS(object = me3, file = "objects/henikoff/me3_filt.rds")

ac <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = ac_counts,
    fragments = ac,
    annotation = annot
  ),
  assay = "ac",
  meta.data = md_ac
)

ac <- FindTopFeatures(ac)
ac <- RunTFIDF(ac)
ac <- RunSVD(ac, scale.embeddings = TRUE, reduction.name = 'lsi.ac')
ac <- RunUMAP(ac, reduction = 'lsi.ac', dims = 2:30, reduction.name = 'umap.ac')

saveRDS(object = ac, file = "objects/henikoff/ac_filt.rds")
