NO_PROCESSING_VANILLA <- function(obj, input_path, output_path, 
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
  
  print("Reading epiAneufinder count matrix!")
  epianeufinder_assay_path <- paste0(epiAneufinder_dir_path, "/count_summary.rds")
  epianeufinder_assay_obj <- readRDS(epianeufinder_assay_path)
  
  epianeufinder_ranges <- epianeufinder_assay_obj@rowRanges
  epianeufinder_ranges <- paste0(epianeufinder_ranges$wSeq, "-", 
                                 epianeufinder_ranges$wStart, "-", 
                                 epianeufinder_ranges$wEnd)
  
  epianeufinder_count_matrix <- epianeufinder_assay_obj@assays@data@listData$counts
  
  print("Processing epiAneufinder count matrix!")

  epianeufinder_count_matrix <- epianeufinder_count_matrix[, colnames(epianeufinder_count_matrix) %in% all_cells]
  epianeufinder_count_matrix <- as(epianeufinder_count_matrix, "CsparseMatrix")
  epianeufinder_count_matrix <- epianeufinder_count_matrix[, order(colnames(epianeufinder_count_matrix))] # fix order of columns in CNA_norm_matrix
  
  obj_metadata <- obj@meta.data[order(obj$Original_barcode), ] # fix order of rows in obj@meta.data
  
  if (all(colnames(epianeufinder_count_matrix) == obj$Original_barcode[obj$Original_barcode %in% colnames(epianeufinder_count_matrix)])) { # now order is the same
    colnames(epianeufinder_count_matrix) <- names(obj$Original_barcode)[obj$Original_barcode %in% colnames(epianeufinder_count_matrix)] # we can assign correct barcode name
    rownames(epianeufinder_count_matrix) <- epianeufinder_ranges
    print("Done!")
  }
  
  return(epianeufinder_count_matrix)
}
