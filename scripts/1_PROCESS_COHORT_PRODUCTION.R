################ LOAD LIBRARIES ################

source("/home/ieo7429/Desktop/THESIS_GAB/scripts/HELPER_FUNCTIONS.R")
import_libraries()

################

################ LOAD HELPERS ################

files_paths <- c("/home/ieo7429/Desktop/THESIS_GAB/scripts/H1_EPIANEUFINDER_BIN_FILTERING.R",
                 "/home/ieo7429/Desktop/THESIS_GAB/scripts/H1_NO_PROCESSING_VANILLA.R",
                 "/home/ieo7429/Desktop/THESIS_GAB/scripts/H2.1_REGRESS_OUT_COHORT.R",
                 "/home/ieo7429/Desktop/THESIS_GAB/scripts/H1_LIFT_FILE.R",
                 "/home/ieo7429/Desktop/THESIS_GAB/scripts/my_makeWindows.R",
                 "/home/ieo7429/Desktop/THESIS_GAB/scripts/obsolete/PROCESS_ONE_SAMPLE.R")

import_scripts(files_paths = files_paths)

################

################ INITIALIZE VARIABLES #############

load_variables(tables_dir = "/home/ieo7429/Desktop/THESIS_GAB/tables/",
               scripts_dir = "/home/ieo7429/Desktop/THESIS_GAB/scripts/",
               plots_outdir = "/home/ieo7429/Desktop/THESIS_GAB/plots/", 
               files_outdir = "/home/ieo7429/Desktop/THESIS_GAB/outfiles/",
               samples_folder_path = "/home/ieo7429/Desktop/THESIS_GAB/samples/BRCA_BASAL/",
               match_tumor_cnc_path = paste0(TABLES_DIR, "MATCH_TUMOR_CNC"), 
               cohort_ids_file_path = paste0(TABLES_DIR, "BRCA_BASAL_IDS"), 
               blacklist_path = paste0(TABLES_DIR,"hg19_blacklist.v2.bed"), 
               chain_path = paste0(TABLES_DIR,"hg38ToHg19.over.chain"),
               unlifted_peaks_path = paste0(TABLES_DIR,"pancan_peaks"),
               lifted_peaks_path = paste0(TABLES_DIR,"pancan_peaks_lifted"),
               bash_script = paste0(SCRIPTS_DIR,"filter_fragments.sh"),
               conda_env = "GAB_ENV",
               complementary = "no",
               genome_assembly_str = "BSgenome.Hsapiens.UCSC.hg19", 
               genome_assembly = BSgenome.Hsapiens.UCSC.hg19, 
               genome_assembly_code = "hg19", 
               window_size = 1000000, 
               cancer_type = "BRCA_BASAL", 
               workers = 80, 
               alpha = 0.05, 
               barplot_title = "barplot_inspection_", 
               markers_list_name = "markers_list_", 
               collapsed_DACRs_df_name = "collapsed_table_", 
               whole_genome_plot_name = "manhattan_whole_genome_", 
               mode = "no_processing_",
               with_epianeufinder = TRUE,
               lift_files = TRUE,
               make_plots = TRUE,
               test.use = "LR",
               exclude = c("chrY","chrM"),
               new_assay_name = "bin_level_counts",
               backbone_path = paste0(TABLES_DIR,"backbone_",WINDOW_SIZE_STRING,".tsv"))

################

################ COHORT INSPECTION ###############

# read SeuratObject file names
setwd(SAMPLES_FOLDER_PATH)
files <- list.files(recursive = TRUE, pattern = "HT*/*")
files_name <- grep(pattern = "rds", files, value = TRUE)

# make a list of SeuratObjects
list_of_rds <- lapply(files_name, readRDS)
names(list_of_rds) <- sample_ids

# already inspected peaks, they are all the same

# make a barplot to visualize
if (BOOL_MAKE_PLOTS) {
  
  barplot_inspection <- plot_barplot_inspection(list_of_rds, sample_ids, CNC, CANCER_TYPE)

  png(filename = paste0(PLOTS_OUTDIR, barplot_title, CANCER_TYPE, ".png"), width = 7000, height = 6000, res = 800)
  print(barplot_inspection)
  dev.off()
  
  rm(barplot_inspection)
}

################

################ CONSIDERATIONS AFTER COHORT INSPECTION ################ 

# Observations FOR BRCA_NON_BASAL: 

# HT088B1-S1H1 is not worth merging with HT088B1-S1H2
# DISCARD: HT088B1-S1H1, HT128, HT163, HT214, HT217, HT235, HT297, HT497, HT545
# KEEP: HT088B1-S1H2, HT137, HT243, HT263, HT305, HT514

# Final cohort will be of n = 6 patients

# subset the list
if (CANCER_TYPE == "BRCA_NON_BASAL") {
  sample_ids <- c("HT088B1-S1H2", "HT137B1-S1H7", "HT243B1-S1H4", 
                  "HT263B1-S1H1", "HT305B1-S1H1", "HT514B1-S1H3")
} else if (CANCER_TYPE == "BRCA_BASAL"){

# Observations FOR BRCA_BASAL: 

# DISCARD: HT029, HT035, HT1408-06, HT141, HT271, HT378B1-S1H1, HT378B1-S1H2, HT384B1-S1H1
# KEEP: HT206, HT517

# Final cohort will be of n = 2 patients


# AFTER PROCESSING OF THE TWO COHORTS, THEY WILL BE MERGED GENERATING A UNIQUE BRCA_COHORT

# subset the list
sample_ids <- c("HT206B1-S1H4", "HT517B1-S1H1")}

list_of_rds <- list_of_rds[sample_ids]
names(list_of_rds) <- sample_ids
gc()

################

################ LIFTING OVER GENOMIC COORDINATES ################

if (BOOL_LIFT_FILES) {
  # lifting peaks
  LIFT_FILE(file_path = UNLIFTED_PEAKS_PATH, 
            chain_path = CHAIN_PATH, 
            mode = "peaks")
  
  
  # lifting fragments file
  sapply(X = sample_ids, FUN = function(ID){
    
    selected_fragment_file <- FRAGMENT_FILES_PATHS[which(startsWith(x = FRAGMENT_FILES_PATHS, prefix = ID))] 
    
    LIFT_FILE(file_path = selected_fragment_file, 
              chain_path = CHAIN_PATH, 
              mode = "fragments")
    
  })
}

################

################ FILTERING FRAGMENTS FILE ################

params <- c(CANCER_TYPE, LIFTED_PEAKS_PATH, CONDA_ENV, COMPLEMENTARY)

system2(BASH_SCRIPT, params)

################

################ ADJUST EPIANEUFINDER TOOL SOURCE CODE ################

# assign the source code the modified version of makeWindows in order to account for correct windows definition
assignInNamespace("makeWindows", 
                  my_makeWindows, 
                  asNamespace("epiAneufinder"))

################

################ CONSIDERATIONS ABOUT COPY NUMBER ALTERATIONS ################ 

# It appears that CNAs are correlated to Chromatin accessibility.
# So there are different approaches to fix it:

#  - Decorrelate CNAs from Chromatin Accessibility:
#    - Epianeufinder to classify + Quantile Mapping to adjust the distributions (look at previous sections)
#    - Standard processing + Regress out CNA effect (look at next sections)

#  - Epianeufinder to classify + NAs imputation (look at next sections)
#  - Epianeufinder to classify + bin removal (look at next sections)

################

################ SET UP PROCESSING FUNCTIONS ################ 

if (CANCER_TYPE == "BRCA_NON_BASAL") {
  minFrags_vec <- c(1000, 1000, 1000, 1000, 1000, 2000) # for BRCA_NON_BASAL
} else if (CANCER_TYPE == "BRCA_BASAL"){
  minFrags_vec <- c(1000, 1000)}

if (MODE == "no_processing_" && BOOL_WITH_EPIANEUFINDER){
  PROCESSING_FUNCTION <- NO_PROCESSING_VANILLA # no correction
} else if (MODE == "no_processing_" && !BOOL_WITH_EPIANEUFINDER){
  PROCESSING_FUNCTION <- PROCESS_ONE_SAMPLE
} else if (MODE == "bin_filtering_"){
  PROCESSING_FUNCTION <- EPIANEUFINDER_BIN_FILTERING # bin filtering
}

################ 

################ COHORT PROCESSING WITH EPIANEUFINDER ###############
if (BOOL_WITH_EPIANEUFINDER) {
  
  n_of_samples <- length(sample_ids)
  
  markers_list <- list()
  for (idx in seq_along(1:n_of_samples)) {
    
    obj <- list_of_rds[[idx]]
    sample_id <- sample_ids[[idx]]
    minFrags <- minFrags_vec[[idx]]
    
    file_names_in_sample_folder <- list.files(path = sample_id)
    fragments_mask_filtered <- grep(pattern = "*lifted_filtered.tsv.gz", file_names_in_sample_folder)
    fragments_path_filtered <- file_names_in_sample_folder[fragments_mask_filtered]
    fragments_path_filtered <- paste0(sample_id, "/", fragments_path_filtered)
    
    output_path <- paste0(sample_id,"/")
    
    epianeufinder_corrected_matrix <- PROCESSING_FUNCTION(obj = obj, # input Seurat object
                                                          input_path = fragments_path_filtered, # input path of the fragments file
                                                          output_path = output_path, # where to store the output
                                                          blacklist = BLACKLIST_PATH, # ENCODE blacklisted regions
                                                          window_size = WINDOW_SIZE, # window size
                                                          genome = GENOME_ASSEMBLY_STR , # genome assembly to use
                                                          exclude = exclude, # chromosomes to exclude
                                                          reuse.existing = FALSE, # resume a previously started epiAneufinder run. Bad if parameters change
                                                          ncores = WORKERS, # parallel processing
                                                          minFrags = minFrags, # minimum number of fragments to keep a cell
                                                          minsizeCNV = 1, # minimum bin size of a CNA event
                                                          k = 4, # how many breakpoints (2^k)
                                                          plotKaryo = FALSE,  # make the karyoplot
                                                          folder_name = paste0("epiAneufinder_results_", WINDOW_SIZE_STRING))
    
    # fix colnames before creating chromatin assay
    
    new_assay_name_norm <- paste0(new_assay_name,"_norm")
    obj[[new_assay_name_norm]] <- CreateAssayObject(counts = epianeufinder_corrected_matrix)
    
    obj <- RunTFIDF(obj, assay = new_assay_name_norm) # normalize
    print("I just normalized the counts!")
    
    print("Beginning FindMarkers() execution!")
    markers <- FindMarkers(obj, # perform DA analysis
                           group.by = "cell_type",
                           assay = new_assay_name_norm,
                           test.use = test.use,
                           ident.1 = "Tumor",
                           ident.2 = CNC,
                           logfc.threshold = 0, 
                           min.pct = 0)
    
    print("Done!")
    
    markers_list[[idx]] <- markers # add markers to the list
  }
} else {
################

################ COHORT PROCESSING WITHOUT EPIANEUFINDER ###############

  n_of_samples <- length(sample_ids)
  
  markers_list <- list()
  for (idx in seq_along(1:n_of_samples)) {
    
    obj <- list_of_rds[[idx]]
    sample_id <- sample_ids[[idx]]
    
    file_names_in_sample_folder <- list.files(path = sample_id)
    fragments_mask_filtered <- grep(pattern = "*lifted_filtered.tsv.gz", file_names_in_sample_folder)
    fragments_path_filtered <- file_names_in_sample_folder[fragments_mask_filtered]
    fragments_path_filtered <- paste0(sample_id, "/", fragments_path_filtered)
    
    corrected_matrix <- PROCESSING_FUNCTION(obj = obj, # input Seurat object
                                            input_path = fragments_path_filtered, # input path of the fragments file
                                            blacklist = BLACKLIST_PATH, # ENCODE blacklisted regions
                                            genome_assembly = GENOME_ASSEMBLY_STR,
                                            window_size = WINDOW_SIZE, # window size
                                            exclude = exclude)
    
    # fix colnames before creating chromatin assay
    
    new_assay_name_norm <- paste0(new_assay_name,"_norm")
    obj[[new_assay_name_norm]] <- CreateAssayObject(counts = corrected_matrix)
    
    obj <- RunTFIDF(obj, assay = new_assay_name_norm) # normalize
    print("I just normalized the counts!")
    
    print("Beginning FindMarkers() execution!")
    markers <- FindMarkers(obj, # perform DA analysis
                           group.by = "cell_type",
                           assay = new_assay_name_norm,
                           test.use = test.use,
                           ident.1 = "Tumor",
                           ident.2 = CNC,
                           logfc.threshold = 0, 
                           min.pct = 0)
    
    print("Done!")
    
    markers_list[[idx]] <- markers # add markers to the list
  }
}
################


################ SAVE DATA AND QUIT ################ 

rm(list_of_rds, sample_ids, new_assay_name)

save(markers_list, 
     file = paste0(FILES_OUTDIR, markers_list_name, SUFFIX, "_", test.use, ".RData")) # save as .RData, load it using load()

print("Done!")
quit(save = "no", status = 0)

################








