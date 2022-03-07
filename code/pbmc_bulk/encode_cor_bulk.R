library(ggplot2)
library(patchwork)
if (!requireNamespace("reshape2", quietly = TRUE)) {
  install.packages("reshape2")
}


mat <- read.table("data/pbmc_bulk/encode_cor.tsv", sep = "\t", header = TRUE, comment.char = "")
mat <- mat[, 4:ncol(mat)] # remove region
mat <- as.matrix(mat)
mat[is.nan(mat)] <- 0
mat <- as.data.frame(mat)

colnames(mat) <- c(
  "sc-H3K27me3", "sc-H3K27ac",
  "bulk-H3K27me3 (mono)", "bulk-H3K27me3 (multi)",
  "bulk-H3K27ac (mono)", "bulk-H3K27ac (multi)",
  "H3K27me3 (ChIP)", "H3K27ac (ChIP)",
  "H3K27ac (CT-pro)", "H3K27me3 (CT-pro)"
)

mat <- mat[, c("sc-H3K27me3", "H3K27me3 (CT-pro)", "sc-H3K27ac", "H3K27ac (CT-pro)")]

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

ggsave("plots/pbmc/encode_cor_all.png", plot = p1, height = 3, width = 5)