library(variancePartition)
library(limma)

metadata <- read.table(file = "../../metadata_combined.txt", header = T, tryLogical = F, sep = "\t")
count_matrix <- 