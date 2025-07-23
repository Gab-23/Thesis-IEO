################ LOAD LIBRARIES ################

source("/home/ieo7429/Desktop/THESIS_GAB/scripts/HELPER_FUNCTIONS.R")
import_libraries()

################

################ INITIALIZE VARIABLES #############

load_variables(tables_dir = "/home/ieo7429/Desktop/THESIS_GAB/tables/",
               files_outdir = "/home/ieo7429/Desktop/THESIS_GAB/outfiles/",
               genome_assembly_str = "BSgenome.Hsapiens.UCSC.hg19", 
               genome_assembly = BSgenome.Hsapiens.UCSC.hg19, 
               genome_assembly_code = "hg19", 
               window_size = 100000,
               cancer_type = "BRCA",
               backbone_path = paste0(TABLES_DIR,"backbone_", WINDOW_SIZE_STRING,".tsv"),
               
               training_dataset_path = paste0(TABLES_DIR, 
                                              "training_dataset_wo_target_w_CNA_", 
                                              WINDOW_SIZE_STRING, "_", GENOME_ASSEMBLY_CODE, ".RData"),
               
               final_table_outpath = paste0(FILES_OUTDIR, 
                                            "bulk_ML_tables_regression_classification_BRCA_", 
                                            WINDOW_SIZE_STRING, "_with_HIC_and_Repliseq.RData")
               )


################

################ GENERATE FINAL DATASETS WITH TARGET VARIABLE ################ 

backbone <- read.table(file = backbone_path, header = T, sep = "\t")
load(training_dataset_path)

table_target_no_processing_path <- "/home/ieo7429/Desktop/THESIS_GAB/bulk_ATAC/GSE231481/peaks/DA_output_backbone_0.1Mbp.tsv"
table_target_regrout_path <- "/home/ieo7429/Desktop/THESIS_GAB/bulk_ATAC/GSE231481/peaks/regrout_DA_output_backbone_0.1Mbp.tsv"

table_target_no_processing <- read.table(table_target_no_processing_path, header = T, sep = "\t")
table_target_regrout <- read.table(table_target_regrout_path, header = T, sep = "\t")

table_target_no_processing <- merge(x = backbone, y = table_target_no_processing, by = c("chr", "start", "end"), sort = F) 
table_target_no_processing <- table_target_no_processing %>% 
                                dplyr::select(bin, logFC, bool_diff_acc)

log2k_1 <- log2(1); table_target_no_processing$sign_mean_log2FC_1 <- ifelse(table_target_no_processing$logFC < -log2k_1, -1, 
                                                                      ifelse(table_target_no_processing$logFC > log2k_1, 1, 0))

log2k_2 <- log2(2); table_target_no_processing$sign_mean_log2FC_2 <- ifelse(table_target_no_processing$logFC < -log2k_2, -1, 
                                                                      ifelse(table_target_no_processing$logFC > log2k_2, 1, 0))

log2k_3 <- log2(3); table_target_no_processing$sign_mean_log2FC_3 <- ifelse(table_target_no_processing$logFC < -log2k_3, -1, 
                                                                      ifelse(table_target_no_processing$logFC > log2k_3, 1, 0))

training_dataset_subset <- ML_dataset_with_target[ML_dataset_with_target$Type == CANCER_TYPE,]

hic_matrix_path <- paste0("/home/ieo7429/Desktop/HiC_pipeline_processing/cool/MCF10A_merged/hic_mat_clean_with_bins_", WINDOW_SIZE_STRING, "_hg19.tsv")
hic_matrix <- read.table(hic_matrix_path, header = T, sep = "\t", check.names = F)

repliseq_df_path <- paste0("/home/ieo7429/Desktop/RepliSeq_pipeline_processing/repliseq_postproc_data_", WINDOW_SIZE_STRING, ".tsv")
repliseq_df <- read.table(file = repliseq_df_path, header = T, sep = "\t")
repliseq_df <- repliseq_df %>% dplyr::select(bin, log2smoothed, replication_class)

training_dataset_subset <- merge(x = training_dataset_subset, y = hic_matrix, by = "bin")
training_dataset_subset <- merge(x = training_dataset_subset, y = repliseq_df, by = "bin")

training_dataset_subset_no_processing <- merge(x = table_target_no_processing, y = training_dataset_subset, by = "bin")
training_dataset_subset_regrout <- merge(x = table_target_regrout, y = training_dataset_subset, by = "bin")

ML_tables <- list(no_processing = training_dataset_subset_no_processing,
                  regrout = training_dataset_subset_regrout)

save(ML_tables, file = final_table_outpath)

################

print("Done!")
quit(save = "no", status = 0)


