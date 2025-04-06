source("/home/ieo7429/Desktop/THESIS_GAB/scripts/HELPER_FUNCTIONS.R") # load helpers
import_libraries() # import libraries

load_variables(genome_assembly = BSgenome.Hsapiens.UCSC.hg19,# load variables
               genome_assembly_code = "hg19",
               window_size = 100000,
               seg_file_dir = "/home/ieo7429/Desktop/segmentation_files_all_cohorts/",
               tables_dir = "/home/ieo7429/Desktop/THESIS_GAB/tables/",
               backbone_path = paste0(TABLES_DIR,"backbone_", WINDOW_SIZE_STRING, ".tsv"),
               training_dataset_path = paste0(TABLES_DIR, "Feature_tables_for_ML_models_11cohorts_IntINSIDER.RData"),
               centromere_table_path = paste0(TABLES_DIR, "centomere.tsv")
               )

backbone <- read.table(backbone_path, header = T, sep = "\t") # load backbone

CANCER_TYPES <- c("BRCA","COADREAD","ESCA","GBMLGG","KIRC","KIRP","LUSC","LUAD","OV","PAAD","STAD") # define 11 cancer cohorts
thr_signal <- 0.2 # define threshold signal

load(training_dataset_path) # load training dataset with all windows
Feature.tables.for.ML <- Feature.tables.for.ML[[WINDOW_SIZE_STRING]] # take the windows you want

merged_frequencies_list <- list() # initialize list
for (CANCER_TYPE in CANCER_TYPES) { # iterate over different cohorts
  
  file_name <- paste0(CANCER_TYPE, "_hg19_FireBrowse_totalCopyNumber_classifiedALL.tsv") # take the file_name
  
  tab <- read.table(paste0(SEG_FILE_DIR, file_name), header = T, sep = "\t") # read it
  tab$classified <- NULL # remove useless columns
  tab$Num_Probes <- NULL
  tab$Chromosome <- paste0("chr", tab$Chromosome) # change chr formulation
  tab <- tab[,c(1,7,8,9,2,3,4,10,5,6)] # reorded columns
  
  n_patients <- length(unique(tab$patient_id)) # take the number of patients
  
  tab_amplification <- tab[tab$Segment_Mean > thr_signal,] # table for amplification according to thr_signal
  tab_deletion <- tab[tab$Segment_Mean < -thr_signal,] # table for deletion according to thr_signal
  
  # define granges
  granges_backbone <- GRanges(seqnames = backbone$chr, ranges = IRanges(start = backbone$start, end = backbone$end), bin = backbone$bin)
  granges_amplification <- GRanges(seqnames = tab_amplification$Chromosome, ranges = IRanges(start = tab_amplification$Start, end = tab_amplification$End))
  granges_deletion <- GRanges(seqnames = tab_deletion$Chromosome, ranges = IRanges(start = tab_deletion$Start, end = tab_deletion$End))
  
  # find overlaps (Default: any overlap counts)
  overlaps_amplification<- findOverlaps(query = granges_amplification, subject = granges_backbone, )
  overlaps_deletion<- findOverlaps(query = granges_deletion, subject = granges_backbone)
  
  # convert overlaps to data.frames
  overlaps_amplification_df <- as.data.frame(overlaps_amplification); colnames(overlaps_amplification_df) <- c("tab_amplification", "backbone") 
  overlaps_deletion_df <- as.data.frame(overlaps_deletion); colnames(overlaps_deletion_df) <- c("tab_deletion", "backbone")
  
  # define amplification and deletion frequency tables
  ampl_frequencies <- overlaps_amplification_df %>% group_by(backbone) %>% summarise(ampl_score = (n() / n_patients)) %>% ungroup()
  del_frequencies <- overlaps_deletion_df %>% group_by(backbone) %>% summarise(del_score = (n() / n_patients)) %>% ungroup()
  
  # merge everything based on backbone region
  merged_frequencies <- merge(x = ampl_frequencies, y = del_frequencies, by = "backbone")
  merged_frequencies$bin <- backbone$bin[merged_frequencies$backbone]
  merged_frequencies$Type <- rep(CANCER_TYPE, nrow(merged_frequencies))
  merged_frequencies_list[[CANCER_TYPE]] <- merged_frequencies[,c(5,4,2,3)]
  print(paste0("Finished execution of ",CANCER_TYPE))
  
}; print("Done!")

pancan_merged_frequencies <- merged_frequencies_list %>% bind_rows()
ML_dataset_with_target <- merge(x = pancan_merged_frequencies, y = Feature.tables.for.ML, by = c("Type", "bin"))

chrom_sizes_canon_df <- data.frame(chr = names(chrom_sizes_canon),
                                   Chromosome_Length = chrom_sizes_canon,
                                   row.names = NULL)

centromere_table <- read.table(centromere_table_path, header = T, sep = "")
centromere_table$start <- NULL; centromere_table$end <- NULL
centromere_table$chr <- paste0("chr", centromere_table$chr)

colnames(centromere_table) <- c("chr", "Centromere_Length", "Centromere_Type", "Centromere")

chrom_information <- merge(x = chrom_sizes_canon_df, y = centromere_table, by = "chr")
chrom_information <- merge(x = backbone, y = chrom_information, by = "chr")
chrom_information$chr <- NULL; chrom_information$start <- NULL; chrom_information$end <- NULL; chrom_information$Centromere <- NULL

ML_dataset_with_target <- merge(x = ML_dataset_with_target, y = chrom_information, by = "bin")

save(ML_dataset_with_target, 
     file = paste0(TABLES_DIR, 
                   "training_dataset_wo_target_w_CNA_", 
                   WINDOW_SIZE_STRING, "_", GENOME_ASSEMBLY_CODE, ".RData"))










