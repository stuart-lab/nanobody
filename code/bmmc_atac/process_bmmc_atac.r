library(Signac)
library(Seurat)
library(SummarizedExperiment)
library(future)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)

plan("multicore", workers = 8)
options(future.globals.maxSize = 100 * 1024 ^ 3)

### Load processed object for the metadata
orig.obj <- readRDS("data/bmmc_atac/scATAC-Healthy-Hematopoiesis-191120.rds")
annot <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)
seqlevelsStyle(annot) <- "UCSC"

# get all fragment files
frag.files <- list.files("data/bmmc_atac/", pattern = "^scATAC_", full.names = TRUE)
frag.files <- paste0(frag.files, "/fragments.tsv.gz")
names(frag.files) <-  list.files("data/bmmc_atac/", pattern = "^scATAC_")

# call peaks on each and combine
pk.all <- list()
for (i in seq_along(frag.files)) {
  pks <- CallPeaks(
    object = frag.files[[i]]
  )
  pk.all[[i]] <- pks
}

gr.combined <- Reduce(f = c, x = pk.all)
gr.combined <- GenomicRanges::reduce(gr.combined)
gr.combined <- keepStandardChromosomes(gr.combined, pruning.mode = "coarse")

# quantify and create objects
md <- colData(orig.obj)

obj.list <- list()
for (i in seq_along(frag.files)) {
  obj_ident <- gsub(pattern = "scATAC_", x = names(frag.files)[[i]], replacement = "")
  cells.keep <- md[md$Group == obj_ident, "Barcode"]
  frags <- CreateFragmentObject(
    path = frag.files[[i]],
    cells = cells.keep
  )
  peakcounts <- FeatureMatrix(
    fragments = frags,
    features = gr.combined,
    cells = cells.keep
  )
  obj <- CreateChromatinAssay(counts = peakcounts, fragments = frags)
  obj <- CreateSeuratObject(counts = obj, assay = "ATAC")
  obj <- RenameCells(obj, add.cell.id = obj_ident)
  obj.list[[i]] <- obj
}

rownames(md) <- paste0(md$Group, "_", md$Barcode)
md <- as.data.frame(md)
bmmc <- merge(x = obj.list[[1]], y = obj.list[2:length(obj.list)])
bmmc <- AddMetaData(object = bmmc, metadata = md)
Annotation(bmmc) <- annot

bmmc$tissue <- sapply(strsplit(bmmc$Group, "_"), `[[`, 1)
bmmc$sample <- sapply(strsplit(bmmc$Group, "_"), `[[`, 2)

# remove PBMC
Idents(bmmc) <- "tissue"
bmmc <- subset(bmmc, idents = "PBMC", invert = TRUE)

bmmc <- NucleosomeSignal(bmmc)

# process
bmmc <- FindTopFeatures(bmmc)
bmmc <- RunTFIDF(bmmc)
bmmc <- RunSVD(bmmc)
bmmc <- RunUMAP(bmmc, reduction = "lsi", dims = 2:40)
ga <- GeneActivity(bmmc)
bmmc[["GA"]] <- CreateAssayObject(counts = ga)

DefaultAssay(bmmc) <- "GA"
bmmc <- NormalizeData(bmmc)

DefaultAssay(bmmc) <- "ATAC"

# integrate embeddings
atac <- SplitObject(bmmc, split.by = "Group")
atac <- lapply(atac, RunSVD, scale.embeddings = FALSE)

atac_anchors <- FindIntegrationAnchors(
  object.list = atac,
  reduction = "rlsi",
  dims = 2:30
)

atac_integrated <- IntegrateEmbeddings(
  anchorset = atac_anchors,
  new.reduction.name = "integrated_lsi",
  reductions = bmmc[["lsi"]],
  dims.to.integrate = 1:30
)
bmmc[["integrated_lsi"]] <- atac_integrated[["integrated_lsi"]]

saveRDS(object = bmmc, file = "objects/bmmc_atac.rds")