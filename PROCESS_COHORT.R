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

GENOME_ASSEMBLY <- BSgenome.Hsapiens.UCSC.hg38
MATCH_TUMOR_CNC <- data.frame(tumor_type = c("BRCA_NON_BASAL","BRCA_BASAL", "PDAC", "CRC", 
                                             "UCEC", "OV", "HNSCC", "CESC", "SKCM", "ccRCC", "GBM"),
                              CNC = c("Luminal mature", "Luminal progenitor", "Ductal like 2", 
                                      "Distal stem cells", "Secretory endometrial epithelial cells", 
                                      "Secretory endometrial epithelial cells", "Normal squamous cells", 
                                      "Normal squamous cells", "Melanocytes", "Proximal tubule cells", 
                                      "Oligodendrocyte precursor cells"))

WINDOW_SIZE <-  100000 # user defined
COHORT_IDS_FILE <- "BRCA_NON_BASAL_IDS.txt" # user defined
CANCER_TYPE <- "BRCA_NON_BASAL" # user defined

sample_ids <- readLines(COHORT_IDS_FILE)
CNC <- MATCH_TUMOR_CNC[MATCH_TUMOR_CNC$tumor_type == CANCER_TYPE, "CNC"]
  
chrom_sizes <- seqlengths(GENOME_ASSEMBLY)
chrom_sizes_canon <- chrom_sizes[1:24]
chrom_cum_start <- c(0, cumsum(as.numeric(chrom_sizes_canon[-length(chrom_sizes_canon)])))
names(chrom_cum_start) <- names(chrom_sizes_canon)

windows <- define_windows(window_size = WINDOW_SIZE, 
                          chrom_sizes_canon = chrom_sizes_canon)

markers_list <- list()
plotting_data_list <- list()
for (idx in seq_along(1:length(sample_ids))){
  
  sample_id <- sample_ids[idx]
  markers <- PROCESS_ONE_SAMPLE(sample_id = sample_id, 
                                CNC = CNC, 
                                windows = windows, 
                                new_assay_name = "bin_level_counts", 
                                test.use = "LR")
  DACRs_coords <- get_DACRs_coords(markers)
  
  markers_list[[idx]] <- markers
}

