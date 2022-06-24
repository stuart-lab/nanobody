library(ggplot2)
library(patchwork)

mat_ac <- read.table("data/pbmc_protein/bigwig/raw_ac.tsv", sep = "\t", header = TRUE, comment.char = "")
mat_ac <- mat_ac[, 4:ncol(mat_ac)] # remove region
mat_ac <- as.matrix(mat_ac)
mat_ac[is.nan(mat_ac)] <- 0
mat_ac <- as.data.frame(mat_ac)

mat_me3 <- read.table("data/pbmc_protein/bigwig/raw_me3.tsv", sep = "\t", header = TRUE, comment.char = "")
mat_me3 <- mat_me3[, 4:ncol(mat_me3)] # remove region
mat_me3 <- as.matrix(mat_me3)
mat_me3[is.nan(mat_me3)] <- 0
mat_me3 <- as.data.frame(mat_me3)

encode <- list("ENCFF569IPX.bigWig" = "NK_H3K27me3_chip",
               "ENCFF293ETP.bigWig" = "NK_H3K27ac_chip",
               "ENCFF842JLZ.bigWig" = "CD8_memory_H3K27me3_chip",
               "ENCFF611GRL.bigWig" = "CD8_memory_H3K27ac_chip",
               "ENCFF181OXO.bigWig" = "CD4_H3K27ac_chip",
               "ENCFF811VAX.bigWig" = "CD4_H3K27me3_chip",
               "ENCFF526VJO.bigWig" = "CD8_H3K27ac_chip",
               "ENCFF499VWN.bigWig" = "CD8_H3K27me3_chip",
               "ENCFF777EGG.bigWig" = "CD14_mono_H3K27me3_chip",
               "ENCFF601NLG.bigWig" = "CD14_mono_H3K27ac_chip",
               "ENCFF951ZBV.bigWig" = "B_cell_H3K27ac_chip",
               "ENCFF046VLL.bigWig" = "B_cell_H3K27me3_chip",
               "ENCFF085YMB.bigWig" = "CMP_CD34_H3K27me3_chip",
               "ENCFF064JOI.bigWig" = "CMP_CD34_H3K27ac_chip")

process_matrix <- function(mat, encode, method = "pearson") {
    encode_mat <- colnames(mat) %in% names(encode)
    colnames(mat)[encode_mat] <- unname(encode[colnames(mat)[encode_mat]])
    colnames(mat) <- gsub(".bw", "", colnames(mat), fixed = TRUE)
    
    # work out me3 or ac
    if (any(grepl("me3", colnames(mat)))) {
        mat <- mat[, c("B_cell_me3", "B_cell_H3K27me3_chip",
               "NK_me3", "NK_H3K27me3_chip",
               "CD14_Mono_me3", "CD14_mono_H3K27me3_chip",
               "CD8_T_cell_me3", "CD8_H3K27me3_chip",
               "CD4_T_cell_me3", "CD4_H3K27me3_chip",
               "Late_erythroid_me3", "CMP_CD34_H3K27me3_chip")]
    } else {
        mat <- mat[, c("B_cell_ac", "B_cell_H3K27ac_chip",
               "NK_ac", "NK_H3K27ac_chip",
               "CD14_Mono_ac", "CD14_mono_H3K27ac_chip",
               "CD8_T_cell_ac", "CD8_H3K27ac_chip",
               "CD4_T_cell_ac", "CD4_H3K27ac_chip",
               "Late_erythroid_ac", "CMP_CD34_H3K27ac_chip")]
    }
    all_cor <- cor(mat, method = method)
    ac <- reshape2::melt(all_cor)
    ac <- ac[ac$Var1 %in% encode, ]
    ac <- ac[ac$Var2 %in% setdiff(colnames(mat), encode), ]
    ac$`Pearson correlation` <- ac$value
    return(ac)
}

me3 <- process_matrix(mat=mat_me3, encode=encode, method = "spearman")
ac <- process_matrix(mat=mat_ac, encode=encode, method = "spearman")

cells_keep <- c("CD14_Mono_me3", "Late_erythroid_me3", "B_cell_me3")
encode_keep <- c("B_cell_H3K27me3_chip", "CD14_mono_H3K27me3_chip", "CMP_CD34_H3K27me3_chip")

me3 <- me3[me3$Var2 %in% cells_keep, ]
me3 <- me3[me3$Var1 %in% encode_keep, ]

cells_keep <- c("CD14_Mono_ac", "Late_erythroid_ac", "B_cell_ac")
encode_keep <- c("B_cell_H3K27ac_chip", "CD14_mono_H3K27ac_chip", "CMP_CD34_H3K27ac_chip")

ac <- ac[ac$Var2 %in% cells_keep, ]
ac <- ac[ac$Var1 %in% encode_keep, ]

replacement_encode <- list("B_cell_H3K27ac_chip" = "B cell",
                           "CD14_mono_H3K27ac_chip" = "CD14 Monocyte",
                           "CMP_CD34_H3K27ac_chip" = "CD34+ CMP")

replacement_sc <- list("B_cell_ac" = "B cell",
                       "CD14_Mono_ac" = "CD14 Monocyte",
                       "Late_erythroid_ac" = "Late_erythroid")

ac$Var1 <- as.character(replacement_encode[as.character(ac$Var1)])
ac$Var2 <- as.character(replacement_sc[as.character(ac$Var2)])

replacement_encode <- list("B_cell_H3K27me3_chip" = "B cell",
                           "CD14_mono_H3K27me3_chip" = "CD14 Monocyte",
                           "CMP_CD34_H3K27me3_chip" = "CD34+ CMP")

replacement_sc <- list("B_cell_me3" = "B cell",
                       "CD14_Mono_me3" = "CD14 Monocyte",
                       "Late_erythroid_me3" = "Late_erythroid")

me3$Var1 <- as.character(replacement_encode[as.character(me3$Var1)])
me3$Var2 <- as.character(replacement_sc[as.character(me3$Var2)])

p1 <- ggplot(me3, aes(Var1, Var2, fill = `Pearson correlation`)) +
  geom_tile() +
  theme_classic() +
  scale_fill_distiller(palette = "RdYlBu") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  xlab("Bulk-cell ChIP-seq") +
  ylab("Single-cell multiplexed NTT-seq") +
  ggtitle("H3K27me3") +
  theme(legend.position = "none")

p2 <- ggplot(ac, aes(Var1, Var2, fill = `Pearson correlation`)) +
  geom_tile() +
  theme_classic() +
  scale_fill_distiller(palette = "RdYlBu") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  xlab("Bulk-cell ChIP-seq") +
  ylab("Single-cell multiplexed NTT-seq") +
  ggtitle("H3K27ac")

pp <- (p1 / p2)

ggsave(filename = "plots/pbmc/encode_cor_spearman.png", plot = pp, height = 6, width = 5)

# correlation between replicates
mat <- read.table("data/pbmc_protein/replicate_cor.tsv", sep = "\t", header = TRUE, comment.char = "")
mat <- mat[, 4:ncol(mat)] # remove region
mat <- as.matrix(mat)
mat[is.nan(mat)] <- 0
mat <- as.data.frame(mat)

colnames(mat) <- c(
  "Rep. 1 H3K27me3",
  "Rep. 1 H3K27ac",
  "Rep. 2 H3K27me3",
  "Rep. 2 H3K27ac"
)

mat <- mat[, c("Rep. 1 H3K27me3", "Rep. 2 H3K27me3", "Rep. 1 H3K27ac", "Rep. 2 H3K27ac")]

all_cor <- cor(mat, method = "pearson")
ac <- reshape2::melt(all_cor)
ac$`Pearson correlation` <- ac$value

p1 <- ggplot(ac, aes(Var1, Var2, fill = `Pearson correlation`)) +
  geom_tile() +
  theme_classic() +
  scale_fill_distiller(palette = "RdBu") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  xlab("") +
  ylab("")

ggsave("plots/pbmc/replicate_correlation.png", plot = p1, height = 3, width = 5)
