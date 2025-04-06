import_scripts <- function(files_paths){
  for (file in files_paths) {
    source(file)
  }
} # source all the helper scripts

import_libraries <- function(){
  
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(BSgenome.Hsapiens.UCSC.hg19)
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
  library(parallel)
  library(presto)
  library(tidyr)
  library(abind)
  library(ggnewscale)
  library(ggrepel)
  library(qmap)
  library(epiAneufinder)
  library(impute)
  library(rtracklayer)
  library(data.table)
  library(Hmisc)
  library(rCGH)
  
} # function to import all needed libraries

load_variables <- function(match_tumor_cnc_path = NULL, # load all needed variables
                           cohort_ids_file_path = NULL, 
                           blacklist_path = NULL, 
                           samples_folder_path = NULL, 
                           unlifted_peaks_path = NULL, 
                           chain_path = NULL, 
                           plots_outdir = NULL, 
                           files_outdir = NULL,
                           seg_file_dir = NULL,
                           tables_dir = NULL,
                           genome_assembly_str = NULL, 
                           genome_assembly = NULL, 
                           genome_assembly_code = NULL, 
                           window_size = NULL, 
                           cancer_type = NULL, 
                           workers = NULL, 
                           alpha = NULL,
                           barplot_title = NULL,
                           markers_list_name = NULL,
                           collapsed_DACRs_df_name = NULL,
                           whole_genome_plot_name = NULL, 
                           mode = NULL,
                           lift_files = NULL,
                           make_plots = NULL,
                           simple_avg = NULL,
                           test.use = NULL,
                           exclude = NULL,
                           new_assay_name = NULL,
                           corrected_df_name = NULL,
                           correlation_plots_name = NULL,
                           accessibility_table_path = NULL,
                           ampl_del_table_path = NULL,
                           backbone_path = NULL,
                           table_no_processing_path = NULL,
                           table_regrout_path = NULL,
                           table_bin_filtering_path = NULL,
                           table_no_processing_name = NULL,
                           table_regrout_name = NULL,
                           table_bin_filtering_name = NULL,
                           ML.Tables_path = NULL,
                           table_target_no_processing_path = NULL,
                           table_target_regrout_path = NULL,
                           table_target_bin_filtering_path = NULL,
                           final_table_outpath = NULL,
                           training_dataset_path = NULL,
                           centromere_table_path = NULL){
  
  PLOTS_OUTDIR <<- plots_outdir # plots directories
  FILES_OUTDIR <<- files_outdir # files directories
  TABLES_DIR <<- tables_dir # tables directory
  SEG_FILE_DIR <<- seg_file_dir # segmentation files directory
  
  MATCH_TUMOR_CNC_PATH <<- match_tumor_cnc_path # matching tumor type and corresponding CNC
  COHORT_IDS_FILE_PATH <<- cohort_ids_file_path # cohort members' IDs file
  BLACKLIST_PATH <<- blacklist_path # ENCODE blacklisted regions
  SAMPLES_FOLDER_PATH <<- samples_folder_path # samples directories location
  
  if (!is.null(SAMPLES_FOLDER_PATH)) {
    FRAGMENT_FILES_PATHS <<- list.files(path = SAMPLES_FOLDER_PATH, 
                                     pattern = "*fragments.tsv.gz", 
                                     recursive = TRUE)
  }
  
  
  UNLIFTED_PEAKS_PATH <<- unlifted_peaks_path
  
  CHAIN_PATH <<- chain_path
  
  # read tables
  if (!is.null(MATCH_TUMOR_CNC_PATH)) {
    MATCH_TUMOR_CNC <<- read.table(MATCH_TUMOR_CNC_PATH, header = T, sep = "\t") 
  }
  
  # define the genome assembly
  GENOME_ASSEMBLY_STR <<- genome_assembly_str
  GENOME_ASSEMBLY <<- genome_assembly
  GENOME_ASSEMBLY_CODE <<- genome_assembly_code
  
  # define the window size
  WINDOW_SIZE <<-  window_size
  
  if (!is.null(WINDOW_SIZE)) {
    WINDOW_SIZE_STRING <<- retrieve_window_size_string(WINDOW_SIZE)
  }
  
  # define cancer type and CNC
  CANCER_TYPE <<- cancer_type
  if (!is.null(MATCH_TUMOR_CNC_PATH)) {
    CNC <<- MATCH_TUMOR_CNC[grep(pattern = CANCER_TYPE, MATCH_TUMOR_CNC$tumor_type), "CNC"]
  if (length(CNC) > 1) {
    CNC <<- paste0(CNC[1], " - ", CNC[2])  
    }
  }
  
  # how many cores to use during parallel executed operations?
  WORKERS <<- workers
  
  # significance level to evaluate DACRs
  ALPHA <<- alpha
  
  # suffix for files and plots
  if (!is.null(WINDOW_SIZE)) {
    SUFFIX <<- paste0(CANCER_TYPE, "_", WINDOW_SIZE_STRING, "_", GENOME_ASSEMBLY_CODE)
  }
  
  # titles of files and plots
  
  if (!is.null(mode)) {
    MODE <<- mode
    barplot_title <<- barplot_title
    markers_list_name <<- paste0(markers_list_name, mode)
    collapsed_DACRs_df_name <<- paste0(collapsed_DACRs_df_name, mode)
    whole_genome_plot_name <<- paste0(whole_genome_plot_name, mode)
    corrected_df_name <<- corrected_df_name
    correlation_plots_name <<- correlation_plots_name
    accessibility_table_path <<- accessibility_table_path
    ampl_del_table_path <<- ampl_del_table_path
  }
  
  if(!is.null(COHORT_IDS_FILE_PATH)){
    sample_ids <<- readLines(COHORT_IDS_FILE_PATH) # read the cohort IDs
  }
  
  if (!is.null(GENOME_ASSEMBLY)) {
    chrom_sizes <- seqlengths(GENOME_ASSEMBLY) # get chromosome sizes
    chrom_sizes_canon <<- chrom_sizes[1:22] # subset to canon chromosomes
    chrom_cum_start <<- c(0, cumsum(as.numeric(chrom_sizes_canon[-length(chrom_sizes_canon)]))) # get cumulative starts
    names(chrom_cum_start) <<- names(chrom_sizes_canon)
    
    BOOL_LIFT_FILES <<- lift_files
    BOOL_MAKE_PLOTS <<- make_plots
    BOOL_SIMPLE_AVG <<- simple_avg
    
    test.use <<- test.use
    exclude <<- exclude
    new_assay_name <<- new_assay_name
    
    backbone_path <<- backbone_path
    
    table_no_processing_path <<- table_no_processing_path
    table_regrout_path <<- table_regrout_path
    table_bin_filtering_path <<- table_bin_filtering_path
    
    table_no_processing_name <<- table_no_processing_name
    table_regrout_name <<- table_regrout_name
    table_bin_filtering_name <<- table_bin_filtering_name
    
    ML.Tables_path <<- ML.Tables_path
    
    table_target_no_processing_path <<- table_target_no_processing_path
    table_target_regrout_path <<- table_target_regrout_path
    table_target_bin_filtering_path <<- table_target_bin_filtering_path
    
    final_table_outpath <<- final_table_outpath
    
    training_dataset_path <<- training_dataset_path
    
    centromere_table_path <<- centromere_table_path
    
  }
}

retrieve_window_size_string <- function(window_size){
  if (window_size == 1000000) {
    return("1Mbp")
  }
  
  if (window_size == 100000) {
    return("0.1Mbp")
  }
  
  if (window_size == 50000) {
    return("50Kbp")
  }
} # convert the window size integer to a string representation

retrieve_rownames_of_GRanges_object <- function(GRanges_object){
  names_GRanges_object <- paste(seqnames(GRanges_object), 
                                start(GRanges_object), 
                                end(GRanges_object),
                                sep = "-")
} # retrieve string representation of GRanges object

define_windows <- function(window_size, chrom_sizes_canon){
  
  division_factor <- round(chrom_sizes_canon / window_size)
  per_chromosome_bin_size <- floor(chrom_sizes_canon / division_factor)
  small_windows_idxs <- which(per_chromosome_bin_size < window_size)
  
  per_chromosome_bin_size[small_windows_idxs] <- floor(chrom_sizes_canon[small_windows_idxs] / 
                                                         (division_factor[small_windows_idxs] - 1))
  
  windows <- lapply(X = seq_along(1:length(chrom_sizes_canon)), FUN = function(idx){
    
    chrom_windows <- tileGenome(seqlengths = chrom_sizes_canon[idx],
                                tilewidth = per_chromosome_bin_size[idx],
                                cut.last.tile.in.chrom = T)
    
    bin_widths <- unique(chrom_windows@ranges@width)
    
    if (length(bin_widths) > 1) {
      value_bin_to_remove <- min(bin_widths)
      chrom_windows <- chrom_windows[chrom_windows@ranges@width != value_bin_to_remove]
    }
    
    return(chrom_windows)
    
  })
  
  windows <- GRangesList(windows)
  windows <- unlist(windows)
  genome(windows) <- "hg19"
  names(windows) <- retrieve_rownames_of_GRanges_object(windows)
  
  return(windows)
} # given a window size and chromosome size, take windows as GRanges

generate_backbone <- function(windows){
  
  internal_idx <- 1
  curr_chromosome <- "chr1"
  
  bin_ids <- c()
  for (idx in seq_along(1:length(windows))) {
    
    range <- windows[idx]
    
    chr <- as.character(seqnames(range))
    chr_order <- which(chr == levels(seqnames(range)))
    
    if (chr != curr_chromosome) {
      internal_idx <- 1
    }
    
    bin_ids <- c(bin_ids, paste0(chr_order,"_",internal_idx))
    
    internal_idx <- internal_idx + 1
    curr_chromosome <- chr
    
  }
  coords_df <- as.data.frame(windows)
  coords_df$bin <- bin_ids
  coords_df$width <- NULL
  coords_df$strand <- NULL
  rownames(coords_df) <- NULL
  colnames(coords_df) <- c("chr", "start", "end", "bin")
  return(coords_df)
} # starting from a GRanges object, generate a backbone file, associating genomic coordinates and bin_id

filter_windows <- function(windows, original_peaks){
  
  n_of_overlaps <- countOverlaps(query = windows,
                                 subject = original_peaks, 
                                 type = "any")
  windows_filtered <- windows[n_of_overlaps > 0]
  return(windows_filtered)
} # keep only those windows that overlap at least one peak

retrieve_original_barcodes <- function(seurat_obj){
  sample_level_names <- colnames(seurat_obj)
  sample_names_split <- do.call(what = rbind, strsplit(sample_level_names, split = "_"))
  raw_names <- sample_names_split[,3]
  return(raw_names)
} # obsolete, but from Seurat Object take the original barcodes (without sample name)

handle_fragments <- function(fragments_path, seurat_obj, original_barcodes){
  
  file_name_parts <- strsplit(fragments_path, "\\.")[[1]][1:2] # divide the fragments path into parts
  fragments_bgz_path <- paste0(file_name_parts[1], ".", file_name_parts[2],".","bgz") # generate .bgz path
  fragments_idx_path <- paste0(fragments_bgz_path,".tbi") # generate .tbi path
  
  if (!(file.exists(fragments_path))) { # check for existence of files / presence of .bgz file
    stop("fragments file does not exist")
  } else if (file.exists(fragments_bgz_path)) {
    message("fragments file is already bgzipped")
  } else {bgzip(fragments_path)} # if not exist bgzip the file
  
  if (!(file.exists(fragments_bgz_path))) { # check for existence of  files / presence of .tbi file
    stop("bgzipped fragments file does not exist")
  } else if (file.exists(fragments_idx_path)) {
    message("fragments file index already exists")  
  } else {indexTabix(fragments_bgz_path, format = "bed")} # if not exist index it
  
  fragment_obj <- CreateFragmentObject(path = fragments_bgz_path, cells = original_barcodes) # create the fragments object
  return(fragment_obj)
} # function to deal with the fragments file

resize_count_matrix <- function(fragment_obj, windows, cells_barcodes, sample_level_names){ # reassign fragments from file to a defined cell of windows
  
  resized_count_matrix <- FeatureMatrix(fragments = fragment_obj,
                                        features = windows, 
                                        cells = cells_barcodes)
  
  colnames(resized_count_matrix) <- sample_level_names
  return(resized_count_matrix)
}

get_DACRs_coords <- function(markers){ 
  
  DACRs <- rownames(markers)
  DACRs_coords <- as.data.frame(do.call(what = rbind, sapply(X = DACRs, FUN = strsplit, split = "-")))
  colnames(DACRs_coords) <- c("chr","start","end")
  DACRs_coords$chr <- factor(DACRs_coords$chr, levels = paste0("chr",c(as.character(1:22), "X","Y")))
  DACRs_coords$start <- as.integer(DACRs_coords$start)
  DACRs_coords$end <- as.integer(DACRs_coords$end)
  DACRs_coords$log2FC <- markers$avg_log2FC
  DACRs_coords$pval_adj <- markers$p_val_adj
  DACRs_coords$subtype <- markers$subtype
  DACRs_coords <- DACRs_coords %>% arrange(chr, start)
  
  return(DACRs_coords)
} # process the markers data.frame

process_DACRs_coords <- function(DACRs_coords, chrom_cum_start){
  DACRs_coords <- DACRs_coords %>%
    dplyr::mutate(
      cum_start = start + chrom_cum_start[chr], # define cumulative coordinates
      cum_end = end + chrom_cum_start[chr]
    )
  
  DACRs_coords$absolute_midpoints <- (DACRs_coords$cum_start + DACRs_coords$cum_end) / 2 # define cumulative middle point
  return(DACRs_coords)
} # process the data frame for plotting purposes

process_markers_list <- function(markers_list, chrom_cum_start, n_of_samples, alpha, alpha_classification){
  
  import_libraries()
  
  markers_list_rownames <- lapply(markers_list, rownames) # keep only common regions if different regions are present
  common_rownames <- Reduce(intersect, markers_list_rownames)
  markers_list <- lapply(X = markers_list, FUN = function(x){
    x <- x[common_rownames,]
  })
  
  DACRs_list <- lapply(X = markers_list, FUN = get_DACRs_coords) # apply get_DACRs_coords to all elements
  DACRs_list <- lapply(X = DACRs_list, FUN = function(x){
    x$significant <- ifelse(x$pval_adj < alpha, 1, 0)
    return(x)
  })
  
  coords <- DACRs_list[[1]][,1:3] # take the coordinates
  DACRs_tensor <- abind(lapply(DACRs_list, as.matrix), along = 3) # generate a tensor with all tables
  
  n_of_significant_per_bin <- apply(X = DACRs_tensor, # take the number of patients for which a particular change is significant
                                    MARGIN = 1, 
                                    FUN = function(x){sum(as.numeric(x[7,]))})
  
  perc_of_significant_per_bin <- n_of_significant_per_bin / n_of_samples # take the % of patients for which a change is significant
  
  log2FC_slice <- DACRs_tensor[,4,] # take the log2FC slice n (regions) * 1 (column) * m (samples) 3D slice
  pval_slice <- DACRs_tensor[,5,] # take the adj pval slice n (regions) * 1 (column) * m (samples) 3D slice
  
  mode(log2FC_slice) <- "numeric" # set it to numeric
  mode(pval_slice) <- "numeric"
  
  subtype_vec <- table(DACRs_tensor[1,6,]) # take the cohort subtypes
  
  if (length(subtype_vec) > 1){ # if more than 1 subtype
    
    initial_weight <- min(subtype_vec) / max(subtype_vec) # 2 / 6 = 0.33333
    weight <- c(rep(1, min(subtype_vec)), rep(initial_weight, max(subtype_vec))) # set up the vector
    print("Performing weighted average and sd!")
    
  } else {
    weight <- rep(1, subtype_vec) # if only one subtype, then all weights are 1
    print("Performing simple average and sd!")
  }
  
  weight_slice <- t(replicate(dim(DACRs_tensor)[1], weight)) # retrieve the weight slice
  DACRs_tensor <- abind(DACRs_tensor, weight_slice, along = 2) # add it to the tensor
  
  mean_log2FC <- sapply(X = seq_along(1:dim(DACRs_tensor)[1]), # compute the weighted average
                       FUN = function(x){
                              log2FC_vec <- log2FC_slice[x,]
                              weight_vec <- weight_slice[x,]
                              weighted_avg <- wtd.mean(x = log2FC_vec, 
                                                       weights = weight_vec)
  })
  
  sd_log2FC <- sapply(X = seq_along(1:dim(DACRs_tensor)[1]), # compute the weighted average
                        FUN = function(x){
                          log2FC_vec <- log2FC_slice[x,]
                          weight_vec <- weight_slice[x,]
                          weighted_sd <- sqrt(wtd.var(x = log2FC_vec, 
                                                   weights = weight_vec))
                        })
  
  joint_probability <- sapply(X = seq_along(1:dim(DACRs_tensor)[1]), # compute the weighted average
                              FUN = function(x){
                                      pval_vec <- pval_slice[x,]
                                      fischer_method_joined <- -2 * sum(log(pval_vec)) # fischer method tells that the sum of ln() of a 
                                      aggregated_pvalue <- pchisq(fischer_method_joined, # set of independent pvalues follow a chi-squared distribution with k (num of pvalues) dof
                                                                  df = 2 * length(pval_vec), 
                                                                  lower.tail = FALSE)
                                    })
  
  bool_diff_acc <- ifelse(joint_probability < alpha_classification, 1, 0) 
  # if the joint adjusted pvalues are significant, this mean that the region is a DACR (1), otherwise not (0)
  # the joint differential accessibility would be evaluated using a lower threshold than 0.05
  
  log2k_1 <- log2(1) 
  log2k_2 <- log2(2)
  log2k_3 <- log2(3)
  
  sign_mean_log2FC_1 <- ifelse(mean_log2FC < -log2k_1, -1, ifelse(mean_log2FC > log2k_1, 1, 0))
  sign_mean_log2FC_2 <- ifelse(mean_log2FC < -log2k_2, -1, ifelse(mean_log2FC > log2k_2, 1, 0))
  sign_mean_log2FC_3 <- ifelse(mean_log2FC < -log2k_3, -1, ifelse(mean_log2FC > log2k_3, 1, 0))
  
  weighted_log2FC <- mean_log2FC * perc_of_significant_per_bin # the weighted log2FC for a region is the percentage of patients with a significant change * the value of mean log2FC
  weighted_coef_of_var_log2FC <- (mean_log2FC / sd_log2FC) * perc_of_significant_per_bin
  
  # generate the collapsed data.frame
  collapsed_DACRs_df <- cbind(coords, # coords
                              joint_probability, bool_diff_acc, # classification related variables
                              sign_mean_log2FC_1, sign_mean_log2FC_2, sign_mean_log2FC_3, 
                              mean_log2FC, sd_log2FC, perc_of_significant_per_bin, # regression related variables
                              weighted_log2FC, weighted_coef_of_var_log2FC)
  
  collapsed_DACRs_df <- process_DACRs_coords(collapsed_DACRs_df, chrom_cum_start) # call process_DACRs_coords to continue data.frame processing 
  
  return(collapsed_DACRs_df)
} # process the markers list to obtain a collapsed version over n patients

plot_log2FC_collapsed <- function(DACRs_df, n_samples, alpha, cnc, test, mode) {
  
  if (mode == "mean") {
    p <- ggplot(DACRs_df, aes(x = absolute_midpoints, y = mean_log2FC, 
                              color = mean_log2FC,
                              size = perc_of_significant_per_bin))
    
    important_points <- DACRs_df %>% arrange(desc(abs(weighted_log2FC))) %>% slice_head(prop = 0.2)
    
    ax_str <- "mean_log2FC"
    
  } else if (mode == "weighted") {
    p <- ggplot(DACRs_df, aes(x = absolute_midpoints, y = weighted_log2FC, 
                              color = weighted_log2FC,
                              size = perc_of_significant_per_bin))
    
    important_points <- DACRs_df %>% arrange(desc(abs(weighted_log2FC))) %>% slice_head(prop = 0.2)
    
    ax_str <- "weighted_log2FC"
    
  } else if (mode == "coef_of_var") {
    p <- ggplot(DACRs_df, aes(x = absolute_midpoints, y = weighted_coef_of_var_log2FC, 
                              color = weighted_coef_of_var_log2FC,
                              size = perc_of_significant_per_bin))
    
    important_points <- DACRs_df %>% arrange(desc(abs(weighted_coef_of_var_log2FC))) %>% slice_head(prop = 0.2)
    
    ax_str <- "weighted_coef_of_var_log2FC"
  }
  
  
  chr_midpoints <- DACRs_df %>%
    group_by(chr) %>%
    summarise(midpoint = mean(absolute_midpoints, na.rm = TRUE))
  
  chr_labels <- unique(DACRs_df$chr)
  if (length(chr_labels) > 1) {
    title_chr <- "[WHOLE GENOME]"
  } else {
    title_chr <- paste0("[", chr_labels, "]")
  }
  
  p <- p +
    geom_point() +
    geom_line(data = important_points, color = "black", alpha = 0.8, linewidth = 1) + 
    scale_color_gradientn(colours = c("darkblue", "red", "gold"), values = c(0,0.5,1)) +
    scale_size_continuous(range = c(0.5,6)) +
    theme_minimal() +
    labs(x = "", 
         y = ax_str, 
         title = paste(paste0("Manhattan Plot [Tumor vs ", cnc, "]"), title_chr), 
         subtitle = paste(paste0("Genomic Coordinates vs ", ax_str), 
                          "[ n =", n_samples, "]",
                          "[ \u03B1 =", alpha, "]",
                          "[ test :", test, "]")) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.title.y = element_text(size = 15),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 15), 
      axis.text.y = element_text(size = 15)
    ) +
    scale_x_continuous(
      breaks = chr_midpoints$midpoint,
      labels = chr_midpoints$chr
    ) +
    ylim(-4,4)
  
  return(p)
} # plot chromatin accessibility variation for all samples

plot_barplot_inspection <- function(list_of_rds, sample_ids, cnc, cancer_type){

    num_of_cells <- unlist(lapply(X = list_of_rds, FUN = function(x){dim(x)[2]}))
    num_of_tum_cells <- unlist(lapply(X = list_of_rds, FUN = function(x){table(x$cell_type)["Tumor"]}), use.names = F)
    num_of_CNCs <- unlist(lapply(X = list_of_rds, FUN = function(x){table(x$cell_type)[CNC]}), use.names = F)
    num_of_CNCs <- ifelse(is.na(num_of_CNCs), 0, num_of_CNCs)
  
    cohort_df <- data.frame(sample_id = sample_ids,
                            n_cells = num_of_cells, 
                            n_tum_cells = num_of_tum_cells,
                            n_CNCs = num_of_CNCs)
  
    cohort_df <- pivot_longer(cohort_df, cols = -sample_id, 
                              names_to = "Category", values_to = "Value")
  
    p <- ggplot(cohort_df, aes(x = sample_id, y = Value, fill = Category)) +
         geom_col(position = position_dodge(width = 0.9)) +
         geom_text(aes(label = Value), 
                position = position_dodge(width = 0.9), 
                vjust = -0.5, size = 2.5) +
         theme_minimal() +
         theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
         labs(x = "Sample", 
              y = "Value", 
              fill = "Category",
              title = paste("Barplot of Cell Counts in", cancer_type, "cohort" ))
    
} # plot the barplot for cohort inspection

plot_CNA_acc_correlation <- function(df_before, df_after, mode){
  
  if (mode == "lm") {
    formula <- as.formula(y ~ x)
  } else if (mode == "poly") {
    formula <- as.formula(y ~ poly(x, 2))
  }
  
  ampl_acc_test <- cor.test(df_before$ampl_score, df_before$mean_log2FC, method = "spearman") # perform correlation test. spearman because lack of normality for both variables
  cor_ampl_acc <- round(ampl_acc_test$estimate, 2) # take correlation value
  pval_ampl_acc <- round(ampl_acc_test$p.value, 5) # take pvalue
  
  ampl_acc <- ggplot(df_before, aes(x = ampl_score, y = mean_log2FC, colour = chr)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_smooth(method = "lm", formula = formula, se = FALSE, color = "blue") +
    annotate("text", x = max(df_before$ampl_score) * 0.7, y = max(df_before$ampl_score) * 0.9, 
             label = paste("Correlation: ", cor_ampl_acc, "\n",
                           "P-value: ", pval_ampl_acc), 
             color = "black", size = 5, fontface = "bold") +
    labs(title = "[BEFORE NORMALIZATION]",
         x = "Amplification score",
         y = "Chromatin accessibility [mean Log2FC]") +
    theme_minimal()
  
  ampl_acc_regrout_test <- cor.test(df_after$ampl_score, df_after$mean_log2FC, method = "spearman")
  cor_ampl_acc_regrout <- round(ampl_acc_regrout_test$estimate, 2)
  pval_ampl_acc_regrout <- round(ampl_acc_regrout_test$p.value, 5)
  
  ampl_acc_regrout <- ggplot(df_after, aes(x = ampl_score, y = mean_log2FC, colour = chr)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_smooth(method = "lm", formula = formula, se = FALSE, color = "blue") +
    annotate("text", x = max(df_after$ampl_score) * 0.7, y = max(df_after$ampl_score) * 0.9, 
             label = paste("Correlation: ", cor_ampl_acc_regrout, "\n",
                           "P-value: ", pval_ampl_acc_regrout), 
             color = "black", size = 5, fontface = "bold") +
    labs(title = "[REGRESS OUT AMPLIFICATION + DELETION EFFECT]",
         x = "Amplification score",
         y = "Chromatin accessibility [mean Log2FC]") +
    theme_minimal()
  
  del_acc_test <- cor.test(df_before$del_score, df_before$mean_log2FC, method = "spearman") # perform correlation test. spearman because lack of normality for both variables
  cor_del_acc <- round(del_acc_test$estimate, 2) # take correlation value
  pval_del_acc <- round(del_acc_test$p.value, 5) # take pvalue
  
  del_acc <- ggplot(df_before, aes(x = del_score, y = mean_log2FC, colour = chr)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_smooth(method = "lm", formula = formula, se = FALSE, color = "blue") +
    annotate("text", x = max(df_before$del_score) * 0.7, y = max(df_before$del_score) * 0.9, 
             label = paste("Correlation: ", cor_del_acc, "\n",
                           "P-value: ", pval_del_acc), 
             color = "black", size = 5, fontface = "bold") +
    labs(title = "[BEFORE NORMALIZATION]",
         x = "Deletion score",
         y = "Chromatin accessibility [mean Log2FC]") +
    theme_minimal()
  
  del_acc_regrout_test <- cor.test(df_after$del_score, df_after$mean_log2FC, method = "spearman")
  cor_del_acc_regrout <- round(del_acc_regrout_test$estimate, 2)
  pval_del_acc_regrout <- round(del_acc_regrout_test$p.value, 5)
  
  del_acc_regrout <- ggplot(df_after, aes(x = del_score, y = mean_log2FC, colour = chr)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_smooth(method = "lm", formula = formula, se = FALSE, color = "blue") +
    annotate("text", x = max(df_after$del_score) * 0.7, y = max(df_after$del_score) * 0.9, 
             label = paste("Correlation: ", cor_del_acc_regrout, "\n",
                           "P-value: ", pval_del_acc_regrout), 
             color = "black", size = 5, fontface = "bold") +
    labs(title = "[REGRESS OUT AMPLIFICATION + DELETION EFFECT]",
         x = "Deletion score",
         y = "Chromatin accessibility [mean Log2FC]") +
    theme_minimal()
  
  gridplot <- gridExtra::grid.arrange(ampl_acc, ampl_acc_regrout, 
                                      del_acc, del_acc_regrout,
                                      nrow = 2, ncol = 2)

  return(gridplot)
  
} # plot correlation between chromatin accessibility and CNA before and after out regression step

plot_CNA_acc_manhattan <- function(merged_ampl_acc){  
  
  merged_ampl_acc$bin <- factor(merged_ampl_acc$bin, levels = unique(merged_ampl_acc$bin))
  
  cor_matrix <- merged_ampl_acc %>%
    dplyr::select(ampl_score, del_score, mean_log2FC) %>%
    cor(method = "spearman", use = "complete.obs")
  
  data_long <- merged_ampl_acc %>%
    dplyr::select(bin, absolute_midpoints, ampl_score, del_score, mean_log2FC) %>%
    pivot_longer(cols = c(ampl_score, del_score, mean_log2FC), 
                 names_to = "Score Type", values_to = "value")  # Renaming for legend
  
  color_map <- c("ampl_score" = "red", "del_score" = "blue", "mean_log2FC" = "purple")
  
  cor_labels <- data_long %>%
    group_by(`Score Type`) %>%
    reframe(
      x_pos = max(absolute_midpoints) * 0.95,
      y_pos = max(value, na.rm = TRUE) * 0.9,
      cor_text = case_when(
        `Score Type` == "ampl_score" ~ paste0("Corr: ", round(cor_matrix["ampl_score", "del_score"], 2), " (del) | ",
                                              round(cor_matrix["ampl_score", "mean_log2FC"], 2), " (mean_log2FC)"),
        `Score Type` == "del_score" ~ paste0("Corr: ", round(cor_matrix["del_score", "ampl_score"], 2), " (ampl) | ",
                                             round(cor_matrix["del_score", "mean_log2FC"], 2), " (mean_log2FC)"),
        `Score Type` == "mean_log2FC" ~ paste0("Corr: ", round(cor_matrix["mean_log2FC", "ampl_score"], 2), " (ampl) | ",
                                               round(cor_matrix["mean_log2FC", "del_score"], 2), " (del)")
      )
    )
  
  cor_labels <- unique(cor_labels)
  
  ggplot(data_long, aes(x = absolute_midpoints, y = value, color = `Score Type`)) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_color_manual(values = color_map, labels = c("Amplification Score", "Deletion Score", "mean_log2FC")) +
    facet_wrap(~`Score Type`, ncol = 1, scales = "free_y") + 
    geom_text(data = cor_labels, aes(x = x_pos, y = y_pos, label = cor_text), 
              inherit.aes = FALSE, size = 4, hjust = 1, vjust = 1) +
    theme_minimal() +
    labs(title = "CNA Scores vs Chromatin Accessibility Variation", x = "Genomic Coordinates", y = "Score", color = "Score Type") + 
    theme(
      strip.text = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 12),  
      legend.text = element_text(size = 10)
    )
} # plot a stacked manhattan plot with amplification/deletion scores vs chromatin accessibility variation

