PROCESS_ONE_SAMPLE <- function(obj,
                               sample_id,
                               CNC, 
                               windows, 
                               new_assay_name = "bin_level_counts", 
                               test.use = "LR"){
  
  # LIST FILES IN SAMPLE DIRECTORY
  file_names_in_sample_folder <- list.files(path = sample_id)
  
  fragments_mask <- grep(pattern = "*.tsv.gz", file_names_in_sample_folder)
  
  if ((length(fragments_mask) == 0)) {
    stop("fragments.tsv.gz file does not exist!")    
  }
  
  # DEFINE PATHS
  fragments_path <- paste0(sample_id,"/",file_names_in_sample_folder[fragments_mask])

  # GET RAW AND SAMPLE-LEVEL CELL BARCODES
  raw_names <- retrieve_original_barcodes(seurat_obj = obj)
  sample_level_names <- colnames(obj)
  
  # RETRIEVE ORIGINAL PEAKS
  original_peak_ranges <- obj@assays$pancan@ranges
  
  # FILTER WINDOWS NOT OVERLAPPING WITH PEAKS
  filtered_windows <- filter_windows(windows = windows,
                                     original_peaks = original_peak_ranges)
  print("I just filtered the windows!")
  
  # LOAD THE FRAGMENTS OBJECT
  fragment_obj <- handle_fragments(fragments_path = fragments_path, 
                                   seurat_obj = obj, 
                                   original_barcodes = raw_names)
  print("I just took care of the fragments object!")
  
  # RECOMPUTE THE COUNT MATRIX
  resized_count_matrix <- resize_count_matrix(fragment_obj = fragment_obj, 
                                              windows = filtered_windows, 
                                              cells_barcodes = raw_names, 
                                              sample_level_names = sample_level_names)
  print("I just recomputed the count matrix!")
  
  # CREATE ASSAY OBJECT AND SET IT AS DEFAULT
  obj[[new_assay_name]] <- CreateAssayObject(counts = resized_count_matrix)
  DefaultAssay(obj) <- new_assay_name
  
  # PERFORM TFIDF NORMALIZATION ON NEW COUNT MATRIX
  obj <- RunTFIDF(obj, 
                  assay = new_assay_name)
  print("I just normalized the counts!")
  
  
  # RUN FindMarkers() FUNCTION
  print("Beginning FindMarkers() execution!")
  markers <- FindMarkers(obj,
                         group.by = "cell_type",
                         assay = new_assay_name,
                         test.use = test.use,
                         ident.1 = "Tumor",
                         ident.2 = CNC,
                         logfc.threshold = 0, 
                         min.pct = 0)
  print("Done!")
  
  return(markers)
}


