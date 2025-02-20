library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(GenomicRanges)
library(Seurat)
library(Signac)
library(Rsamtools)
library(rtracklayer)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(MAST)
library(plotly)
library(parallel)
library(presto)
library(tidyr)
library(abind)


MATCH_TUMOR_CNC <- data.frame(tumor_type = c("BRCA_NON_BASAL","BRCA_BASAL", "PDAC", "CRC", 
                                             "UCEC", "OV", "HNSCC", "CESC", "SKCM", "ccRCC", "GBM"),
                              CNC = c("Luminal mature", "Luminal progenitor", "Ductal like 2", 
                                      "Distal stem cells", "Secretory endometrial epithelial cells", 
                                      "Secretory endometrial epithelial cells", "Normal squamous cells", 
                                      "Normal squamous cells", "Melanocytes", "Proximal tubule cells", 
                                      "Oligodendrocyte precursor cells"))

GENOME_ASSEMBLY <- BSgenome.Hsapiens.UCSC.hg38
WINDOW_SIZE <-  50000
COHORT_IDS_FILE <- "BRCA_NON_BASAL_IDS"
CANCER_TYPE <- "BRCA_NON_BASAL"
CNC <- MATCH_TUMOR_CNC[MATCH_TUMOR_CNC$tumor_type == CANCER_TYPE, "CNC"]
WORKERS <- 20
  
sample_ids <- readLines(COHORT_IDS_FILE)

chrom_sizes <- seqlengths(GENOME_ASSEMBLY)
chrom_sizes_canon <- chrom_sizes[1:24]
chrom_cum_start <- c(0, cumsum(as.numeric(chrom_sizes_canon[-length(chrom_sizes_canon)])))
names(chrom_cum_start) <- names(chrom_sizes_canon)

windows <- define_windows(window_size = WINDOW_SIZE, 
                          chrom_sizes_canon = chrom_sizes_canon)

################ COHORT INSPECTION ###############

# load all the .rds objects

HT088B1_S1H1 <- readRDS("HT088B1-S1H1/BRCA_HT088B1-S1H1.rds")
HT088B1_S1H2 <- readRDS("HT088B1-S1H2/BRCA_HT088B1-S1H2.rds")
HT128 <- readRDS("HT128B1-S1H4/BRCA_HT128B1-S1H4.rds")
HT137 <- readRDS("HT137B1-S1H7/BRCA_HT137B1-S1H7.rds")
HT163 <- readRDS("HT163B1-S1H6/BRCA_HT163B1-S1H6.rds")
HT214 <- readRDS("HT214B1-S1H2/BRCA_HT214B1-S1H2.rds")
HT217 <- readRDS("HT217B1-S1H1/BRCA_HT217B1-S1H1.rds")
HT235 <- readRDS("HT235B1-S1H1/BRCA_HT235B1-S1H1.rds")
HT243 <- readRDS("HT243B1-S1H4/BRCA_HT243B1-S1H4.rds")
HT263 <- readRDS("HT263B1-S1H1/BRCA_HT263B1-S1H1.rds")
HT297 <- readRDS("HT297B1-S1H1/BRCA_HT297B1-S1H1.rds")
HT305 <- readRDS("HT305B1-S1H1/BRCA_HT305B1-S1H1.rds")
HT497 <- readRDS("HT497B1-S1H1/BRCA_HT497B1-S1H1.rds")
HT514 <- readRDS("HT514B1-S1H3/BRCA_HT514B1-S1H3.rds")
HT545 <- readRDS("HT545B1-S1H1/BRCA_HT545B1-S1H1.rds")

# make a list out of them

list_of_rds <- list(HT088B1_S1H1, HT088B1_S1H2, HT128, HT137, HT163, HT214, HT217,
                    HT235, HT243, HT263, HT297, HT305, HT497, HT514, HT545)

# already inspected peaks, they are all the same


# get number of cells
num_of_cells <- unlist(lapply(X = list_of_rds, FUN = function(x){dim(x)[2]}))

# get number of tumoral cells
num_of_tum_cells <- unlist(lapply(X = list_of_rds, 
                                  FUN = function(x){table(x$cell_type)["Tumor"]}), 
                           use.names = F)

# get number of CNCs
num_of_CNCs <- unlist(lapply(X = list_of_rds, 
                                  FUN = function(x){table(x$cell_type)[CNC]}), 
                           use.names = F)
num_of_CNCs <- ifelse(is.na(num_of_CNCs), 0, num_of_CNCs)


cohort_df <- data.frame(sample_id = sample_ids,
                        n_cells = num_of_cells,
                        n_tum_cells = num_of_tum_cells,
                        n_CNCs = num_of_CNCs)

# put everything in a df
cohort_df <- pivot_longer(cohort_df, cols = -sample_id, names_to = "Category", values_to = "Value")

# make a barplot to visualize
png(filename = "barplot_inspection.png", width = 4000, height = 2500, res = 400)
barplot_inspection <- ggplot(cohort_df, aes(x = sample_id, y = Value, fill = Category)) +
                              geom_col(position = position_dodge(width = 0.9)) +
                              geom_text(aes(label = Value), 
                                        position = position_dodge(width = 0.9), 
                                        vjust = -0.5, size = 2.5) +
                              theme_minimal() +
                              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                              labs(x = "Sample", y = "Value", fill = "Category")
barplot_inspection
dev.off()

# HT088B1-S1H1 is not worth merging with HT088B1-S1H2

# DISCARD: HT088B1-S1H1, HT128, HT163, HT214, 
         # HT217, HT235, HT297, HT497, HT545

# KEEP: HT088B1-S1H2, HT137, HT243, HT263, HT305, HT514

# Final cohort will be of n = 6 patients

list_of_rds <- list(HT088B1_S1H2, HT137, HT243, HT263, HT305, HT514)
sample_ids <- c("HT088B1-S1H2", "HT137B1-S1H7", "HT243B1-S1H4", 
                "HT263B1-S1H1", "HT305B1-S1H1", "HT514B1-S1H3")
n_of_samples <- length(sample_ids)
new_assay_name <- "bin_level_counts"
test.use = "wilcox"


################ COHORT PROCESSING ###############


markers_list <- mclapply(X = seq_along(1:length(sample_ids)), 
                         FUN = function(x){
  
  obj <- list_of_rds[[x]]
  sample_id <- sample_ids[x]
  markers <- PROCESS_ONE_SAMPLE(obj = obj,
                                sample_id = sample_id,
                                CNC = CNC, 
                                windows = windows, 
                                new_assay_name = new_assay_name, 
                                test.use = test.use)
  return(markers)
}, mc.set.seed = T, mc.cores = WORKERS)

################ PROCESS MARKERS ###############

DACRs_list <- lapply(X = markers_list, FUN = get_DACRs_coords)
DACRs_list <- lapply(X = DACRs_list, FUN = function(x){
  x$significant <- ifelse(x$pval_adj < 0.05, 1, 0)
  return(x)
})

coords <- DACRs_list[[1]][,1:3]
DACRs_tensor <- abind(lapply(DACRs_list, as.matrix), along = 3)

n_of_significant_per_bin <- apply(X = DACRs_tensor, 
                                  MARGIN = 1, 
                                  FUN = function(x){sum(as.numeric(x[6,]))})

n_of_significant_per_bin <- n_of_significant_per_bin / n_of_samples

log2FC_slice <- DACRs_tensor[,4,]
mode(log2FC_slice) <- "numeric"

mean_log2FC <- apply(X = log2FC_slice, MARGIN = 1, FUN = mean)

weighted_log2FC <- mean_log2FC * n_of_significant_per_bin # non so se ha senso
collapsed_DACRs_df <- cbind(coords, n_of_significant_per_bin, mean_log2FC, weighted_log2FC)

collapsed_DACRs_df <- process_DACRs_coords(collapsed_DACRs_df, chrom_cum_start)

png(filename = "manhattan_chr17.png", width = 4000, height = 2500, res = 400)
manhattan_plot_chr17 <- plot_log2FC_collapsed(collapsed_DACRs_df[collapsed_DACRs_df$chr == "chr17",], thr = 0.05)
manhattan_plot_chr17
dev.off()
