source("/home/ieo7429/Desktop/THESIS_GAB/scripts/HELPER_FUNCTIONS.R") # load helpers
import_libraries() # import libraries

load_variables(genome_assembly = BSgenome.Hsapiens.UCSC.hg19, # load variables
               genome_assembly_code = "hg19",
               window_size = 100000,
               files_outdir = "/home/ieo7429/Desktop/THESIS_GAB/outfiles/",
               tables_dir = "/home/ieo7429/Desktop/THESIS_GAB/tables/",
               backbone_path = paste0(TABLES_DIR,"backbone_", WINDOW_SIZE_STRING, ".tsv"),
               bins_path_cooler = "/home/ieo7429/Desktop/HiC_pipeline_processing/cool/MCF10A_merged/MCF10A_merged_bins_hg19_0.1Mbp",
               hic_matrix_cooler = "/home/ieo7429/Desktop/HiC_pipeline_processing/cool/MCF10A_merged/MCF10A_merged_zscore_hg19_0.1Mbp.cool",
               hic_matrix_name = "/home/ieo7429/Desktop/HiC_pipeline_processing/cool/MCF10A_merged/hic_mat_clean_with_bins")

backbone <- read.table(backbone_path, header = T, sep = "\t") # load backbone
bins_cooler <- read.table(bins_path_cooler, header = F, sep = "\t", col.names = c("chr","start","end","weight"))
bins_cooler_with_ids <- merge(x = bins_cooler, y = backbone, by = c("chr", "start", "end"))

hic_matrix <- read.table(hic_matrix_cooler, header = F, sep = "\t")

if (exists("bins_cooler_with_ids")) {
  hic_matrix <- cbind(bins_cooler_with_ids$bin, hic_matrix)
  colnames(hic_matrix) <- c("bin", bins_cooler_with_ids$bin)
}

all_NaNs_mask_rows <- apply(X = hic_matrix[,-1], MARGIN = 1, FUN = function(x){!all(is.na(x))})
all_NaNs_mask_cols <- apply(X = hic_matrix, MARGIN = 2, FUN = function(x){!all(is.na(x))})
hic_matrix_noNANs <- hic_matrix[all_NaNs_mask_rows, all_NaNs_mask_cols]

var_array <- apply(X = hic_matrix_noNANs[,-1], MARGIN = 2, FUN = function(x){var(x)})

k = 200
top_k_var <- names(sort(var_array, decreasing = TRUE)[1:k])

hic_mat_final <- hic_matrix_noNANs[,top_k_var]
hic_mat_final <- cbind(hic_matrix_noNANs$bin, hic_mat_final)
colnames(hic_mat_final)[1] <- "bin"

write.table(hic_mat_final, 
            file = paste0(hic_matrix_name, SUFFIX, ".tsv"), 
            sep = "\t", row.names = F, col.names = T)
