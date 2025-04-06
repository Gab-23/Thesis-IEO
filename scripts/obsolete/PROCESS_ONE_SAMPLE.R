PROCESS_ONE_SAMPLE <- function(obj,
                               sample_id,
                               CNC, 
                               windows, 
                               new_assay_name = "bin_level_counts", 
                               test.use = "LR"){
  
  # LIST FILES IN SAMPLE DIRECTORY
  file_names_in_sample_folder <- list.files(path = sample_id)
  
  fragments_mask_filtered <- grep(pattern = "*_filtered.tsv.gz", file_names_in_sample_folder)
  
  if ((length(fragments_mask_filtered) == 0)) {
    stop("filtered_fragments.tsv.gz file does not exist!")    
  }
  
  
  # DEFINE PATHS
  fragments_path_filtered <- paste0(sample_id,"/",file_names_in_sample_folder[fragments_mask_filtered])
  
  
  # GET RAW AND SAMPLE-LEVEL CELL BARCODES
  raw_names <- retrieve_original_barcodes(seurat_obj = obj)
  sample_level_names <- colnames(obj)
  
  # RETRIEVE ORIGINAL PEAKS
  original_peak_ranges <- obj@assays$pancan@ranges
  
  # FILTER WINDOWS NOT OVERLAPPING WITH PEAKS
  filtered_windows <- filter_windows(windows = windows,
                                     original_peaks = original_peak_ranges)
  print("I just filtered the windows!")
  
  # LOAD THE FRAGMENTS OBJECTS
  fragment_obj_filtered <- handle_fragments(fragments_path = fragments_path_filtered, 
                                            seurat_obj = obj, 
                                            original_barcodes = raw_names)
  
  
  print("I just took care of the fragment objects!")
  
  # RECOMPUTE THE COUNT MATRIX
  
  resized_count_matrix_filtered <- resize_count_matrix(fragment_obj = fragment_obj_filtered, 
                                                       windows = filtered_windows, 
                                                       cells_barcodes = raw_names, 
                                                       sample_level_names = sample_level_names)
  
  print("I just recomputed the count matrices!")
  
  # CREATE ASSAY OBJECT AND SET IT AS DEFAULT
  
  new_assay_name_filtered <- paste0(new_assay_name, "_filtered")
  
  obj[[new_assay_name_filtered]] <- CreateAssayObject(counts = resized_count_matrix_filtered)
  
  DefaultAssay(obj) <- new_assay_name_filtered
  
  # PERFORM TFIDF NORMALIZATION ON NEW COUNT MATRIX
  # USE scRNA seq to normalize ??? Matrix is less sparse now
  
  obj <- RunTFIDF(obj, assay = new_assay_name_filtered)
  
  print("I just normalized the counts!")
  
  
  # RUN FindMarkers() FUNCTION
  
  print("Beginning FindMarkers() execution!")
  markers <- FindMarkers(obj,
                         group.by = "cell_type",
                         assay = new_assay_name_filtered,
                         test.use = test.use,
                         ident.1 = "Tumor",
                         ident.2 = CNC,
                         logfc.threshold = 0, 
                         min.pct = 0)
  
  print("Done!")
  return(markers)
}


