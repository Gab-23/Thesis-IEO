PROCESS_ONE_SAMPLE <- function(obj,
                               input_path,
                               blacklist_path, 
                               genome_assembly_str,
                               window_size = 1000000,
                               exclude = c('chrY','chrM')){

  # GET RAW AND SAMPLE-LEVEL CELL BARCODES
  raw_names <- retrieve_original_barcodes(seurat_obj = obj)
  sample_level_names <- colnames(obj)
  
  # GET BLACKLISTED GRANGES
  
  blacklist <- GenomicRanges::GRanges(geneviewer::read_bed(blacklist_path))
  
  # DEFINE WINDOWS
  windows <- my_makeWindows(genome = genome_assembly_str,
                            blacklist = blacklist,
                            windowSize = window_size,
                            exclude = exclude)
  
  keep <- which(windows$N < 0.001)
  windows <- windows[keep]
  
  # LOAD THE FRAGMENTS OBJECTS
  fragment_obj <- handle_fragments(fragments_path = input_path, 
                                   seurat_obj = obj, 
                                   original_barcodes = raw_names)
  
  print("I just took care of the fragment objects!")
  
  # RECOMPUTE THE COUNT MATRIX
  
  resized_count_matrix <- resize_count_matrix(fragment_obj = fragment_obj, 
                                              windows = windows, 
                                              cells_barcodes = raw_names, 
                                              sample_level_names = sample_level_names)
  
  print("I just recomputed the count matrices!")
  
  return(resized_count_matrix)
}


