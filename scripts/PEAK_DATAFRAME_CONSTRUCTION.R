library(GenomicFeatures)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(dplyr)

peaks_path <- "/home/ieo7429/Desktop/THESIS_GAB/bulk_ATAC/GSE231481/peaks/DA_output_peaks.tsv"
gene_list_path <- "/home/ieo7429/Desktop/THESIS_GAB/tables/cancerGeneList.tsv"
SCREEN_path <- "/home/ieo7429/Desktop/THESIS_GAB/tables/GRCh37-cCREs.bed"

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes <- genes(txdb)

peaks_DA_analysis <- read.table(peaks_path, header = T, sep = "\t", row.names = 1)
peaks_DA_analysis$peakID <- paste0("peak_",1:nrow(peaks_DA_analysis))
peaks_DA_analysis <- peaks_DA_analysis[peaks_DA_analysis$fdr <= 0.05,]
peaks_DA_analysis_gr <- GRanges(peaks_DA_analysis)

gene_list <- read.delim(gene_list_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
gene_list <- gene_list %>% dplyr::select(`Entrez Gene ID`, `Is Oncogene`, `Is Tumor Suppressor Gene`)
gene_list$`Is Tumor Suppressor Gene` <- ifelse(gene_list$`Is Tumor Suppressor Gene` == "Yes", 1, 0)
gene_list$`Is Oncogene` <- ifelse(gene_list$`Is Oncogene` == "Yes", 1, 0)

genes <- genes[mcols(genes)$gene_id %in% gene_list$`Entrez Gene ID`]

nearest_OncoKBGene <- nearest(peaks_DA_analysis_gr, genes)
nearest_distance_OncoKBGene <- distance(peaks_DA_analysis_gr, genes[nearest_OncoKBGene])

map_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = names(genes), columns = c("SYMBOL", "GENENAME"), keytype = "ENTREZID")

SCREEN_annot <- read.table(SCREEN_path, header = F, sep = "\t", col.names = c("chr","start","end","annot"))
SCREEN_annot_gr <- GRanges(SCREEN_annot)

hits_SCREEN <- findOverlaps(query = peaks_DA_analysis_gr, subject = SCREEN_annot_gr)
peaks_hits <- peaks_DA_analysis_gr[queryHits(hits_SCREEN)]
SCREEN_hits <- SCREEN_annot_gr[subjectHits(hits_SCREEN)]

SCREEN_overlaps_summary <- data.frame(peakID = mcols(peaks_hits)$peakID,
                                      annot  = mcols(SCREEN_hits)$annot
                                      )

non_overlapping_peaks <- peaks_DA_analysis[which(!peaks_DA_analysis$peakID %in% SCREEN_overlaps_summary$peakID), "peakID"]
NA_vec <- rep(x = NA, times = length(non_overlapping_peaks))

SCREEN_overlaps_summary <- rbind(SCREEN_overlaps_summary, data.frame(peakID = non_overlapping_peaks, annot = NA_vec))
SCREEN_overlaps_summary$annot_clean <- ifelse(is.na(SCREEN_overlaps_summary$annot), "NONE", SCREEN_overlaps_summary$annot)
dummy_matrix <- model.matrix(~ annot_clean - 1, data = SCREEN_overlaps_summary)
dummy_matrix <- dummy_matrix[,colnames(dummy_matrix) != "annot_cleanNONE", drop = FALSE]

binary_peaks <- data.frame(cbind(SCREEN_overlaps_summary$peakID, dummy_matrix))
colnames(binary_peaks) <- c("peakID", 
                            gsub(pattern = "annot_clean", 
                                 replacement = "", 
                                 x = colnames(dummy_matrix)))

binary_peaks[,-1] <- lapply(X = binary_peaks[,-1], FUN = as.numeric)

binary_peaks <- binary_peaks %>%
                  group_by(peakID) %>%
                  summarise(across(everything(), sum), .groups = "drop")

nearest_df <- data.frame(peakID = peaks_DA_analysis$peakID,
                         chr = peaks_DA_analysis$chr,
                         start = peaks_DA_analysis$start, 
                         end = peaks_DA_analysis$end,
                         logFC = peaks_DA_analysis$logFC,
                         entrez_id = map_ids[nearest_OncoKBGene, "ENTREZID"],
                         gene_symbol = map_ids[nearest_OncoKBGene, "SYMBOL"],
                         dist.to.closest.OncoKBGene = nearest_distance_OncoKBGene
                         )

nearest_df_with_annot <- merge(x = nearest_df, by.x = "entrez_id",
                               y = gene_list,  by.y = "Entrez Gene ID", 
                               all.x = TRUE)

nearest_df_with_annot$entrez_id <- NULL

final_tab <- merge(x = nearest_df_with_annot, y = binary_peaks, by = "peakID")

colnames(final_tab) <- c("peakID", "chr", "start", "end", 
                         "logFC", 
                         "closest.OncoKBGene", "dist.to.closest.OncoKBGene", "is_og", "is_tsg",
                         "n_overlaps_CA", "n_overlaps_CA-CTCF", "n_overlaps_CA-H3K4Me3", "n_overlaps_CA-TF", 
                         "n_overlaps_dELS", "n_overlaps_pELS", "n_overlaps_PLS", "n_overlaps_TF")









