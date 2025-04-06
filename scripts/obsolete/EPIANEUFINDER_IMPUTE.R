EPIANEUFINDER_IMPUTE <- function(obj, input_path, output_path, 
                                        blacklist, window_size = 1000000,
                                        genome = "BSgenome.Hsapiens.UCSC.hg38", 
                                        exclude = c('chrY','chrM'), reuse.existing = FALSE,
                                        ncores = 30, minFrags = 1000, 
                                        minsizeCNV = 0, k = 4, plotKaryo = FALSE){
  
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
  tumor_cells <- all_cells[!(all_cells %in% non_tumor_cells)]
  
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
  
  ################ COHORT PROCESSING WITH EPIANEUFINDER (NAs IMPUTING) ############
  
  CNA_mask <- epianeufinder_table != 1
  table_CNA_idx <- which(CNA_mask, arr.ind = T)
  
  epianeufinder_count_matrix[table_CNA_idx] <- NA
  
  imputed_matrix <- impute::impute.knn(data = epianeufinder_count_matrix, k = round(sqrt(dim(obj)[1])))
  epianeufinder_count_matrix <- imputed_matrix$data
  
  # range <- seq_along(1:nrow(epianeufinder_count_matrix))
  
  # for (idx in range) {
    # row <- epianeufinder_count_matrix[idx, ]
    # row_no_NAs <- row[!is.na(row)]
    
    # median_row <- mean(row_no_NAs)
    # row[which(is.na(row))] <- median_row
    
    # epianeufinder_count_matrix[idx, ] <- row  
  # }
  
  epianeufinder_count_matrix <- as(epianeufinder_count_matrix, "CsparseMatrix")
  epianeufinder_count_matrix <- epianeufinder_count_matrix[, order(colnames(epianeufinder_count_matrix))] # fix order of columns in CNA_norm_matrix
  
  obj_metadata <- obj@meta.data[order(obj$Original_barcode), ] # fix order of rows in obj@meta.data
  
  
  if (all(colnames(epianeufinder_count_matrix) == obj$Original_barcode[obj$Original_barcode %in% colnames(epianeufinder_count_matrix)])) { # now order is the same
    colnames(epianeufinder_count_matrix) <- names(obj$Original_barcode)[obj$Original_barcode %in% colnames(epianeufinder_count_matrix)] # we can assign correct barcode name
    print("Done!")
  }
  
  return(epianeufinder_count_matrix)
  
}
 