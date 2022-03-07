library(Signac)
library(Seurat)

ac <- readRDS("data/ct_pro/H3K27ac.rds")
Fragments(ac) <- NULL
Fragments(ac) <- CreateFragmentObject(
    path = "data/ct_pro/H3K27ac_fragments.tsv.gz",
    cells = Cells(ac)
)
saveRDS(object = ac, file = "data/ct_pro/H3K27ac_updated.rds")

me <- readRDS("data/ct_pro/H3K27me3.rds")
Fragments(me) <- NULL
Fragments(me) <- CreateFragmentObject(
    path = "data/ct_pro/H3K27me3_fragments.tsv.gz",
    cells = Cells(me)
)
saveRDS(object = me, file = "data/ct_pro/H3K27me3_updated.rds")