EPIANEUFINDER_QUANTILE_MAPPING <- function(obj, input_path, output_path, 
                                     blacklist, window_size = 1000000,
                                     genome = "BSgenome.Hsapiens.UCSC.hg38", 
                                     exclude = c('chrY','chrM'), reuse.existing = FALSE,
                                     ncores = 30, minFrags = 1000, 
                                     minsizeCNV = 0, k = 4, plotKaryo = FALSE,
                                     method = "RQUANT"){
  
  epiAneufinder_dir_path <- paste0(output_path, "epiAneufinder_results")
  
  if (dir.exists(epiAneufinder_dir_path)) {
    print("epiAneufinder results already found, skipping execution!")
  } else {
    
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
  epianeufinder_count_matrix_path <- paste0(output_path, "/epiAneufinder_results/count_summary.rds")
  epianeufinder_table_path <- paste0(output_path, "/epiAneufinder_results/results_table.tsv")
  
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
  
  ################### CORRECTION #########################
  
  # RATIONALE: Leverage a classification tool to build random samples of bins in cells
  # treat them as random sample and apply Quantile Mapping approaches 
  # to correct depleted/amplified regions
  
  
  print("Beginning CNA correction!")
  
  zeros_mask <- epianeufinder_table == 0 # get the mask
  zeros_idx <- which(zeros_mask, arr.ind = TRUE) # get the indexes
  
  ones_mask <- epianeufinder_table == 1
  ones_idx <- which(ones_mask, arr.ind = TRUE)
  
  twos_mask <- epianeufinder_table == 2
  twos_idx <- which(twos_mask, arr.ind = TRUE)
  
  epianeufinder_count_matrix <- as.matrix(epianeufinder_count_matrix)
  epianeufinder_table <- as.matrix(epianeufinder_table)
  
  zeros <- epianeufinder_count_matrix[zeros_mask] # extract the elements
  ones <- epianeufinder_count_matrix[ones_mask]
  twos <- epianeufinder_count_matrix[twos_mask]
  
  # Quantile Mapping to correct samples of an observed distribution to one of a target distribution
  
  fit_zeros <- fitQmap(obs = ones, mod = zeros, method = method) # fit the model
  zeros_mapped <- doQmap(x = zeros, fit_zeros) # generate the corrected values
  
  fit_twos <- fitQmap(obs = ones, mod = twos, method = method)
  twos_mapped <- doQmap(x = twos, fit_twos)
  
  print("Processing CNA normalized matrix!")
  
  CNA_norm_matrix <- epianeufinder_count_matrix # copy the matrix and assign normalized values
  CNA_norm_matrix[ones_idx] <- ones 
  CNA_norm_matrix[zeros_idx] <- zeros_mapped
  CNA_norm_matrix[twos_idx] <- twos_mapped
  
  CNA_norm_matrix <- as(CNA_norm_matrix, "CsparseMatrix")
  CNA_norm_matrix <- CNA_norm_matrix[, order(colnames(CNA_norm_matrix))] # fix order of columns in CNA_norm_matrix
  
  obj_metadata <- obj@meta.data[order(obj$Original_barcode), ] # fix order of rows in obj@meta.data
  
  
  if (all(colnames(CNA_norm_matrix) == obj$Original_barcode[obj$Original_barcode %in% colnames(CNA_norm_matrix)])) { # now order is the same
    colnames(CNA_norm_matrix) <- names(obj$Original_barcode)[obj$Original_barcode %in% colnames(CNA_norm_matrix)] # we can assign correct barcode name
    print("Done!")
  }
  
  return(CNA_norm_matrix)
}
