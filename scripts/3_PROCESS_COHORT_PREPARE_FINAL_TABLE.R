################ LOAD LIBRARIES ################

source("/home/ieo7429/Desktop/THESIS_GAB/scripts/HELPER_FUNCTIONS.R")
import_libraries()

################

################ LOAD HELPERS ################

files_paths <- c("/home/ieo7429/Desktop/THESIS_GAB/scripts/H1_EPIANEUFINDER_BIN_FILTERING.R",
                 "/home/ieo7429/Desktop/THESIS_GAB/scripts/H1_NO_PROCESSING_VANILLA.R",
                 "/home/ieo7429/Desktop/THESIS_GAB/scripts/H2.1_REGRESS_OUT_COHORT.R",
                 "/home/ieo7429/Desktop/THESIS_GAB/scripts/LIFT_FILE.R",
                 "/home/ieo7429/Desktop/THESIS_GAB/scripts/my_makeWindows.R")

import_scripts(files_paths = files_paths)

################

################ INITIALIZE VARIABLES #############

load_variables(tables_dir = "/home/ieo7429/Desktop/THESIS_GAB/tables/",
               plots_outdir = "/home/ieo7429/Desktop/THESIS_GAB/plots/", 
               files_outdir = "/home/ieo7429/Desktop/THESIS_GAB/outfiles/",
               genome_assembly_str = "BSgenome.Hsapiens.UCSC.hg19", 
               genome_assembly = BSgenome.Hsapiens.UCSC.hg19, 
               genome_assembly_code = "hg19", 
               window_size = 1000000,
               test.use = "LR",
               cancer_type = "BRCA",
               table_no_processing_path = paste0(FILES_OUTDIR, "no_processing/collapsed_table_no_processing_", SUFFIX, "_LR.tsv"),
               table_regrout_path = paste0(FILES_OUTDIR, "regrout/collapsed_table_regrout_", SUFFIX, "_LR.tsv"),
               table_bin_filtering_path = paste0(FILES_OUTDIR, "bin_filter/collapsed_table_bin_filtering_", SUFFIX, "_LR.tsv"),
               table_no_processing_name = "target_no_processing_",
               table_regrout_name = "target_regrout_",
               table_bin_filtering_name = "target_bin_filtering_",
               backbone_path = paste0(TABLES_DIR,"backbone_", WINDOW_SIZE_STRING,".tsv"),
               training_dataset_path = paste0(TABLES_DIR, "training_dataset_wo_target_w_CNA_1Mbp_hg19.RData"),
               table_target_no_processing_path = paste0(FILES_OUTDIR, "target_variables/target_no_processing_", SUFFIX, "_", test.use,".tsv"),
               table_target_bin_filtering_path = paste0(FILES_OUTDIR, "target_variables/target_bin_filtering_", SUFFIX, "_", test.use,".tsv"),
               table_target_regrout_path = paste0(FILES_OUTDIR, "target_variables/target_regrout_", SUFFIX, "_", test.use,".tsv"),
               final_table_outpath = paste0(FILES_OUTDIR, "ML_tables_regression_classification_BRCA_v2.RData"))


################

################ GENERATE FINAL DATASETS WITH TARGET VARIABLE ################ 

backbone <- read.table(file = backbone_path, header = T, sep = "\t")
load(training_dataset_path)

BRCA_table_no_processing <- read.table(table_no_processing_path, header = T, sep = "\t")
BRCA_table_regrout <- read.table(table_regrout_path, header = T, sep = "\t")
BRCA_table_bin_filtering <- read.table(table_bin_filtering_path, header = T, sep = "\t")

to_keep <- c(1,2,3,4,5,6,7,8,9,12,13) # regrout has already been processed
BRCA_table_no_processing <- BRCA_table_no_processing[,to_keep]
BRCA_table_bin_filtering<- BRCA_table_bin_filtering[,to_keep]

BRCA_table_no_processing <- merge(x = backbone, y = BRCA_table_no_processing, by = c("chr", "start", "end")); BRCA_table_no_processing[,c(1,2,3)] <- NULL
BRCA_table_bin_filtering <- merge(x = backbone, y = BRCA_table_bin_filtering, by = c("chr", "start", "end")); BRCA_table_bin_filtering[,c(1,2,3)] <- NULL

write.table(BRCA_table_no_processing, file = paste0(FILES_OUTDIR, table_no_processing_name, SUFFIX, "_", test.use, ".tsv"), sep = "\t", row.names = F, col.names = T)
write.table(BRCA_table_regrout, file = paste0(FILES_OUTDIR, table_regrout_name, SUFFIX, "_", test.use, ".tsv"), sep = "\t", row.names = F, col.names = T)
write.table(BRCA_table_bin_filtering, file = paste0(FILES_OUTDIR, table_bin_filtering_name, SUFFIX, "_", test.use, ".tsv"), sep = "\t", row.names = F, col.names = T)

# conclusions: 
# - will still use linear model to regress out since it smooths better the effect, yet not perfectly
# - will use normal average instead of weighted average, since it preserves better features of BRCA_NON_BASAL

table_target_no_processing <- read.table(table_target_no_processing_path, header = T, sep = "\t")
table_target_regrout <- read.table(table_target_regrout_path, header = T, sep = "\t")
table_target_bin_filtering <- read.table(table_target_bin_filtering_path, header = T, sep = "\t")

training_dataset_subset <- ML_dataset_with_target[ML_dataset_with_target$Type == CANCER_TYPE,]

training_dataset_subset_no_processing <- merge(x = table_target_no_processing, y = training_dataset_subset, by = "bin")
training_dataset_subset_bin_filtering <- merge(x = table_target_bin_filtering, y = training_dataset_subset, by = "bin")
training_dataset_subset_regrout <- merge(x = table_target_regrout, y = training_dataset_subset, by = "bin")
  
ML_tables <- list(no_processing = training_dataset_subset_no_processing,
                  regrout = training_dataset_subset_regrout,
                  bin_filtering = training_dataset_subset_bin_filtering)

save(ML_tables, 
     file = final_table_outpath)

################
