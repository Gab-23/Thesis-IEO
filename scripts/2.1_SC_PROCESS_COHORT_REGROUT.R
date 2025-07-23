################ LOAD LIBRARIES ################

source("/home/ieo7429/Desktop/THESIS_GAB/scripts/HELPER_FUNCTIONS.R")
import_libraries()

################

################ LOAD HELPERS ################

files_paths <- c("/home/ieo7429/Desktop/THESIS_GAB/scripts/H1_EPIANEUFINDER_BIN_FILTERING.R",
                 "/home/ieo7429/Desktop/THESIS_GAB/scripts/H1_NO_PROCESSING_VANILLA.R",
                 "/home/ieo7429/Desktop/THESIS_GAB/scripts/H1_LIFT_FILE.R",
                 "/home/ieo7429/Desktop/THESIS_GAB/scripts/H2.1_REGRESS_OUT_COHORT.R",
                 "/home/ieo7429/Desktop/THESIS_GAB/scripts/my_makeWindows.R")

import_scripts(files_paths = files_paths)

################

################ INITIALIZE VARIABLES #############

load_variables(tables_dir = "/home/ieo7429/Desktop/THESIS_GAB/tables/",
               files_outdir = "/home/ieo7429/Desktop/THESIS_GAB/outfiles/",
               plots_outdir = "/home/ieo7429/Desktop/THESIS_GAB/plots/",
               match_tumor_cnc_path = paste0(TABLES_DIR, "MATCH_TUMOR_CNC"),
               genome_assembly_str = "BSgenome.Hsapiens.UCSC.hg19", 
               genome_assembly = BSgenome.Hsapiens.UCSC.hg19, 
               genome_assembly_code = "hg19", 
               window_size = 100000, 
               cancer_type = "BRCA", 
               accessibility_table_path = paste0(FILES_OUTDIR, "no_processing/collapsed_table_no_processing_BRCA_0.1Mbp_hg19_LR.tsv"),
               mode = "regrout_",
               make_plots = TRUE,
               backbone_path = paste0(TABLES_DIR,"backbone_", WINDOW_SIZE_STRING,".tsv"),
               ampl_del_table_path = paste0(TABLES_DIR, "training_dataset_wo_target_w_CNA_", WINDOW_SIZE_STRING, "_", GENOME_ASSEMBLY_CODE, ".RData"),
               corrected_df_name = "collapsed_table_regrout_",
               test.use = "LR",
               correlation_plots_name = "regrout_corplots_",
               alpha = 0.05)


################

################ [OPTIONAL] REGRESS OUT CNA EFFECT ################

if (MODE == "regrout_") {
  
  regrout_results_list <- REGRESS_OUT_COHORT(accessibility_table_path = accessibility_table_path, 
                                             ampl_del_table_path = ampl_del_table_path, 
                                             backbone_path = backbone_path,
                                             cancer_type = CANCER_TYPE,
                                             mode = "lm")
  
  corrected_data_frame <- regrout_results_list$results

  if (BOOL_MAKE_PLOTS) {
    correlation_plot <- regrout_results_list$plots
    png(filename = paste0(PLOTS_OUTDIR, correlation_plots_name, SUFFIX, "_", test.use, ".png"), width = 6000, height = 4000, res = 400)
    plot(correlation_plot)
    dev.off()
  }
}  

################

################ PLOT VARIATION OF CHROMATIN ACCESSIBILITY ###############

if (BOOL_MAKE_PLOTS) {
  
  manhattan_plot_whole <- plot_log2FC_collapsed(corrected_data_frame, n_samples = 8, ALPHA, CNC, test.use, mode = "mean")
  
  png(filename = paste0(PLOTS_OUTDIR, whole_genome_plot_name, SUFFIX, "_", test.use, ".png"), width = 21000, height = 7000, res = 800)
  print(manhattan_plot_whole)
  dev.off()
  
}

################

################ SAVE DATA ################ 

if (MODE == "regrout_") {
  
  corrected_data_frame$chr <- NULL # set chr to NULL
  corrected_data_frame$sd_log2FC <- NULL
  corrected_data_frame$perc_of_significant_per_bin <- NULL
  corrected_data_frame$cum_start <- NULL
  corrected_data_frame$cum_end <- NULL
  corrected_data_frame$absolute_midpoints <- NULL
  
  write.table(corrected_data_frame, # save already corrected file
              file = paste0(FILES_OUTDIR, corrected_df_name, SUFFIX, "_", test.use, ".tsv"), 
              sep = "\t", 
              row.names = F, col.names = T)

}

################

print("Done!")
quit(save = "no", status = 0)



