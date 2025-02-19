PROCESS_ONE_SAMPLE <- function(sample_id, 
                               CNC, 
                               windows, 
                               new_assay_name = "bin_level_counts", 
                               test.use = "LR"){
  
  # LIST FILES IN SAMPLE DIRECTORY
  file_names_in_sample_folder <- list.files(path = sample_id)
  
  # DEFINE MASK TO FIND WANTED FILE
  obj_mask <- grep(pattern = "*.rds", file_names_in_sample_folder)
  fragments_mask <- grep(pattern = "*.tsv.gz", file_names_in_sample_folder)
  
  if ((length(obj_mask) == 0) || (length(fragments_mask) == 0)) {
    stop("Either fragments.tsv.gz file or .rds file do not exist!")    
  }
  
  # DEFINE PATHS
  obj_path <- paste0(sample_id,"/",file_names_in_sample_folder[obj_mask])
  fragments_path <- paste0(sample_id,"/",file_names_in_sample_folder[fragments_mask])

  # READ RDS OBJECT
  obj <- readRDS(obj_path)
  
  # GET RAW AND SAMPLE-LEVEL CELL BARCODES
  raw_names <- retrieve_original_barcodes(seurat_obj = obj)
  sample_level_names <- colnames(obj)
  
  # RETRIEVE ORIGINAL PEAKS
  original_peak_ranges <- obj@assays$pancan@ranges
  
  # FILTER WINDOWS NOT OVERLAPPING WITH PEAKS
  filtered_windows <- filter_windows(windows = windows,
                                     original_peaks = original_peak_ranges)
  # LOAD THE FRAGMENTS OBJECT
  fragment_obj <- handle_fragments(fragments_path = fragments_path, 
                                   seurat_obj = obj, 
                                   original_barcodes = raw_names)
  
  # RECOMPUTE THE COUNT MATRIX
  resized_count_matrix <- resize_count_matrix(fragment_obj = fragment_obj, 
                                              windows = filtered_windows, 
                                              cells_barcodes = raw_names, 
                                              sample_level_names = sample_level_names)
  
  # CREATE ASSAY OBJECT AND SET IT AS DEFAULT
  obj[[new_assay_name]] <- CreateAssayObject(counts = resized_count_matrix)
  DefaultAssay(obj) <- new_assay_name
  
  # PERFORM TFIDF NORMALIZATION ON NEW COUNT MATRIX
  obj <- RunTFIDF(obj, 
                  assay = new_assay_name)
  
  
  # RUN FindMarkers() FUNCTION
  markers <- FindMarkers(obj,
                         group.by = "cell_type",
                         assay = new_assay_name,
                         test.use = test.use,
                         ident.1 = "Tumor",
                         ident.2 = CNC,
                         logfc.threshold = 0, 
                         min.pct = 0)
  
  return(markers)
}
