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
               samples_folder_path = "/home/ieo7429/Desktop/THESIS_GAB/samples/BRCA_NON_BASAL/",
               match_tumor_cnc_path = paste0(TABLES_DIR, "MATCH_TUMOR_CNC"), 
               cohort_ids_file_path = paste0(TABLES_DIR, "BRCA_NON_BASAL_IDS"), 
               blacklist_path = paste0(TABLES_DIR,"hg19_blacklist.v2.bed"), 
               unlifted_peaks_path = paste0(TABLES_DIR,"pancan_peaks"), 
               chain_path = paste0(TABLES_DIR,"hg38ToHg19.over.chain"), 
               genome_assembly_str = "BSgenome.Hsapiens.UCSC.hg19", 
               genome_assembly = BSgenome.Hsapiens.UCSC.hg19, 
               genome_assembly_code = "hg19", 
               window_size = 1000000, 
               cancer_type = "BRCA", 
               workers = 25, 
               alpha = 0.05, 
               barplot_title = "barplot_inspection_", 
               markers_list_name = "markers_list_", 
               collapsed_DACRs_df_name = "collapsed_table_", 
               whole_genome_plot_name = "manhattan_whole_genome_", 
               mode = "no_processing_",
               lift_files = FALSE,
               make_plots = TRUE,
               test.use = "LR",
               exclude = c("chrY","chrM"),
               new_assay_name = "bin_level_counts",
               backbone_path = paste0(TABLES_DIR,"backbone_", WINDOW_SIZE_STRING,".tsv"),
               simple_avg = TRUE)


################

################ PROCESS MARKERS ###############

# load markers list for BRCA_BASAL and BRCA_NON_BASAL and process them

load(paste0(FILES_OUTDIR, "no_processing/markers_list_no_processing_BRCA_BASAL_1Mbp_hg19_LR.RData"))
markers_list_BASAL <- markers_list

load(paste0(FILES_OUTDIR, "no_processing/markers_list_no_processing_BRCA_NON_BASAL_1Mbp_hg19_LR.RData"))
markers_list_NON_BASAL <- markers_list

if (BOOL_SIMPLE_AVG) {
  
  markers_list_BASAL <- lapply(X = markers_list_BASAL, 
                               FUN = function(df){df$subtype <- "BRCA" 
                               return(df)})
  
  markers_list_NON_BASAL <- lapply(X = markers_list_NON_BASAL, 
                                   FUN = function(df){df$subtype <- "BRCA" 
                                   return(df)})
  
} else {
  
  markers_list_BASAL <- lapply(X = markers_list_BASAL, 
                               FUN = function(df){df$subtype <- "BRCA_BASAL"
                               return(df)})
  
  markers_list_NON_BASAL <- lapply(X = markers_list_NON_BASAL, 
                                   FUN = function(df){df$subtype <- "BRCA_NON_BASAL"
                                   return(df)})
}

markers_list_BRCA <- append(markers_list_BASAL, markers_list_NON_BASAL)
n_of_samples <- length(markers_list_BRCA)

collapsed_DACRs_df <- process_markers_list(markers_list = markers_list_BRCA,
                                           chrom_cum_start = chrom_cum_start,
                                           n_of_samples = n_of_samples,
                                           alpha = ALPHA,
                                           alpha_classification = ALPHA) # collapse n markers data.frame into one summary table

write.table(collapsed_DACRs_df, # save the collapsed df as a table
            file = paste0(FILES_OUTDIR, collapsed_DACRs_df_name, SUFFIX, "_", test.use, ".tsv"), 
            sep = "\t", row.names = F, col.names = T)

################

################ PLOT VARIATION OF CHROMATIN ACCESSIBILITY ###############

if (BOOL_MAKE_PLOTS) {

  manhattan_plot_whole <- plot_log2FC_collapsed(collapsed_DACRs_df, n_of_samples, ALPHA, CNC, test.use, mode = "mean")
  
  png(filename = paste0(PLOTS_OUTDIR, whole_genome_plot_name, SUFFIX, "_", test.use, ".png"), width = 21000, height = 7000, res = 800)
  print(manhattan_plot_whole)
  dev.off()
  
}

################

