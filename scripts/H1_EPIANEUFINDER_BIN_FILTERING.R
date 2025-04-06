EPIANEUFINDER_BIN_FILTERING <- function(obj, input_path, output_path, 
                                     blacklist, window_size = 1000000,
                                     genome = "BSgenome.Hsapiens.UCSC.hg19", 
                                     exclude = c('chrY','chrM'), reuse.existing = FALSE,
                                     ncores = 30, minFrags = 1000, 
                                     minsizeCNV = 0, k = 4, plotKaryo = FALSE,
                                     folder_name){
  
  epiAneufinder_dir_path <- paste0(output_path, folder_name)
  
  if (dir.exists(epiAneufinder_dir_path)) {
    print("epiAneufinder results already found, skipping execution!")
  } else {
    
    epiAneufinder_dir_path <- paste0(output_path, "epiAneufinder_results")
    
    print("Beginning epiAneufinder execution!")
    
    epiAneufinder(input = input_path,
                  outdir = output_path, 
                  blacklist = blacklist,
                  windowSize = window_size, 
                  genome = genome,
                  exclude = exclude, 
                  reuse.existing = reuse.existing,
                  ncores = ncores, 
                  minFrags = minFrags,
                  minsizeCNV = minsizeCNV, 
                  k = k, 
                  plotKaryo = plotKaryo)
    
    print("Finished epiAneufinder execution!")
    
  }
  
  all_cells <- obj$Original_barcode
  non_tumor_cells <- obj@meta.data[obj$cell_type != "Tumor",]$Original_barcode
  
  print("Reading epiAneufinder table and count matrix!")
  epianeufinder_count_matrix_path <- paste0(epiAneufinder_dir_path, "/count_summary.rds")
  epianeufinder_table_path <- paste0(epiAneufinder_dir_path, "/results_table.tsv")
  
  epianeufinder_count_matrix <- readRDS(epianeufinder_count_matrix_path)
  epianeufinder_count_matrix <- epianeufinder_count_matrix@assays@data@listData$counts
  
  print("Processing epiAneufinder table and count matrix!")
  epianeufinder_table <- read.table(epianeufinder_table_path)
  rownames(epianeufinder_table) <- paste0(epianeufinder_table$seq, "-",
                                          epianeufinder_table$start, "-", epianeufinder_table$end)
  
  epianeufinder_table <- epianeufinder_table[,-c(1,2,3)]
  colnames(epianeufinder_table) <- gsub("cell-", "", 
                                        gsub("\\.", "-", 
                                             colnames(epianeufinder_table)) )
  
  epianeufinder_table[, colnames(epianeufinder_table) %in% non_tumor_cells] <- 1 # assume all non tumor cells as diploid
  
  n_cells <- dim(epianeufinder_count_matrix)[2] # use same approach defined in source code to discard low quality bins
  blacklist_thr <- 0.85
  bins_to_discard <- which((rowSums(epianeufinder_count_matrix == 0)) > (n_cells * blacklist_thr))
  
  if (length(bins_to_discard) > 0) {
    epianeufinder_count_matrix <- epianeufinder_count_matrix[-bins_to_discard,]
    print("Discarded excess bins!")
  }
  
  rownames(epianeufinder_count_matrix) <- rownames(epianeufinder_table)
  
  epianeufinder_table <- epianeufinder_table[, colnames(epianeufinder_table) %in% all_cells]
  epianeufinder_count_matrix <- epianeufinder_count_matrix[, colnames(epianeufinder_count_matrix) %in% all_cells]
  
  ################### BIN FILTERING #########################
  
  # RATIONALE: Remove those bins that have frequency of CNA > thr
  
  print("Beginning bin filtering!")
  
  CNA_mask <- epianeufinder_table != 1 # generate a logical mask to identify CNA events
  bin_level_frequency <- rowSums(CNA_mask) / dim(CNA_mask)[2] # get the frequency of cells having non copy number neutral status for a bin
  
  
  q1 <- quantile(bin_level_frequency, 0.25) # retrieve quantiles
  q3 <- quantile(bin_level_frequency, 0.75)
  
  bin_thr <- q3 + (1.5 * (q3 - q1)) # set the threshold for the bin
  
  bin_level_frequency_filtered <- bin_level_frequency[bin_level_frequency <= bin_thr] # filter the frequency vector
  filtered_bins <- names(bin_level_frequency_filtered) # take the filtered bins
  
  epianeufinder_table_filtered <- epianeufinder_table[filtered_bins,] # subset the table
  epianeufinder_count_matrix_filtered <- epianeufinder_count_matrix[filtered_bins,] # subset the matrix
  
  epianeufinder_count_matrix_filtered <- as(epianeufinder_count_matrix_filtered, "CsparseMatrix")
  epianeufinder_count_matrix_filtered <- epianeufinder_count_matrix_filtered[, order(colnames(epianeufinder_count_matrix_filtered))] # fix order of columns in CNA_norm_matrix
  
  obj_metadata <- obj@meta.data[order(obj$Original_barcode), ] # fix order of rows in obj@meta.data
  
  
  if (all(colnames(epianeufinder_count_matrix_filtered) == obj$Original_barcode[obj$Original_barcode %in% colnames(epianeufinder_count_matrix_filtered)])) { # now order is the same
    colnames(epianeufinder_count_matrix_filtered) <- names(obj$Original_barcode)[obj$Original_barcode %in% colnames(epianeufinder_count_matrix_filtered)] # we can assign correct barcode name
    print("Done!")
  }
  
  return(epianeufinder_count_matrix_filtered)
}
