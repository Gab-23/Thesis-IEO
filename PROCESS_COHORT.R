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
library(future)
library(future.apply)


MATCH_TUMOR_CNC <- data.frame(tumor_type = c("BRCA_NON_BASAL","BRCA_BASAL", "PDAC", "CRC", 
                                             "UCEC", "OV", "HNSCC", "CESC", "SKCM", "ccRCC", "GBM"),
                              CNC = c("Luminal mature", "Luminal progenitor", "Ductal like 2", 
                                      "Distal stem cells", "Secretory endometrial epithelial cells", 
                                      "Secretory endometrial epithelial cells", "Normal squamous cells", 
                                      "Normal squamous cells", "Melanocytes", "Proximal tubule cells", 
                                      "Oligodendrocyte precursor cells"))

GENOME_ASSEMBLY <- BSgenome.Hsapiens.UCSC.hg38
WINDOW_SIZE <-  50000
COHORT_IDS_FILE <- "../BRCA_NON_BASAL_IDS.txt"
CANCER_TYPE <- "BRCA_NON_BASAL"
CNC <- MATCH_TUMOR_CNC[MATCH_TUMOR_CNC$tumor_type == CANCER_TYPE, "CNC"]
  
sample_ids <- readLines(COHORT_IDS_FILE)

chrom_sizes <- seqlengths(GENOME_ASSEMBLY)
chrom_sizes_canon <- chrom_sizes[1:24]
chrom_cum_start <- c(0, cumsum(as.numeric(chrom_sizes_canon[-length(chrom_sizes_canon)])))
names(chrom_cum_start) <- names(chrom_sizes_canon)

windows <- define_windows(window_size = WINDOW_SIZE, 
                          chrom_sizes_canon = chrom_sizes_canon)

# check if duplicated
# split duplicated from unique
# process duplicated --> first recompute the matrix --> then merge


# parallel processing
plan(multisession, workers = 4)
markers_list <- future_lapply(sample_ids, function(sample_id){
  markers <- PROCESS_ONE_SAMPLE(sample_id = sample_id, 
                                CNC = CNC, 
                                windows = windows, 
                                new_assay_name = "bin_level_counts", 
                                test.use = "LR")
  return(markers)
}, future.seed = TRUE)

