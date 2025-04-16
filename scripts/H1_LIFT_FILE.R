LIFT_FILE <- function(file_path, chain_path, mode){
  
  options(scipen=999) # remove scientific notation so that big numbers are not written like 1e+08 
  
  file_path_parts <- strsplit(file_path, split = "\\.")[[1]] # take the parts composing the file path
  file_tab <- data.table::fread(input = file_path) # read the table in efficient way
  chain <- import.chain(chain_path) # import conversion chain
  
  granges_unlifted <- GRanges(seqnames = file_tab[[1]],  # with data.table it's possible to call columns like this
                              ranges = IRanges(start = file_tab[[2]],
                                               end = file_tab[[3]]))
  
  if (mode == "fragments") {
    
    colnames(file_tab) <- c("chr", "start", "end", "barcode", "pcr") # set the colnames
    values(granges_unlifted) <- DataFrame(barcode = file_tab$barcode, pcr = file_tab$pcr) # set the metadata
    
    outfile_name <- paste0(paste0(file_path_parts[1], "_lifted"), # set the output filename
                           ".", file_path_parts[2],
                           ".", file_path_parts[3])
    
    lifted_granges <- liftOver(x = granges_unlifted, chain = chain) # lift coordinates
    mappings_to_discard <- which(elementNROWS(lifted_granges) != 1) # check if there are unmapped regions (GRanges of length 0) or multimapping regions (GRanges of length > 1)
    
    lifted_granges_filtered <- lifted_granges[-mappings_to_discard]
    lifted_granges_filtered_unlisted <- unlist(lifted_granges_filtered) # get the GRanges object, note that when unlisting, mapping failures (encoded as empty GRanges) will disappear
    
    lifted_granges_filtered_unlisted_sorted <- sort(lifted_granges_filtered_unlisted)
    invalid_idx <- which(start(lifted_granges_filtered_unlisted_sorted) > end(lifted_granges_filtered_unlisted_sorted))
    
    if (length(invalid_idx) > 0){
      lifted_granges_filtered_unlisted_sorted_clean <- lifted_granges_filtered_unlisted_sorted[-invalid_idx]
    } else {
      lifted_granges_filtered_unlisted_sorted_clean <- lifted_granges_filtered_unlisted_sorted
    } 
    
    failed_or_split_fragments <- (length(granges_unlifted@ranges@width) - # display number of discarded regions
                                    length(lifted_granges_filtered_unlisted_sorted_clean@ranges@width))
    
    summary_width_unlifted <- summary(granges_unlifted@ranges@width) # get distribution summaries
    summary_width_lifted <- summary(lifted_granges_filtered_unlisted_sorted_clean@ranges@width)
    
    print(paste0("# DISCARDED FRAGMENTS DURING LIFTING PROCESS: ", failed_or_split_fragments))
    cat("\n")
    
    print("PRINTING WIDTH DISTRIBUTION SUMMARY BEFORE LIFTING: ")
    print(summary_width_unlifted)
    cat("\n")
    print("PRINTING WIDTH DISTRIBUTION SUMMARY AFTER LIFTING: ")
    print(summary_width_lifted)
    
    df <- as.data.frame(lifted_granges_filtered_unlisted_sorted_clean) # create a df for .bed file writing
    df$strand <- NULL # remove useless columns
    df$width <- NULL
    
    write.table(x = df, 
                file = gzfile(outfile_name), # gzip the file directly
                sep="\t", 
                quote=FALSE, row.names=FALSE, col.names=FALSE)
    
    
  } else if (mode == "peaks") {
    
    colnames(file_tab) <- c("chr", "start", "end")
    outfile_name <- paste0(file_path_parts, "_lifted")
    
    lifted_granges <- liftOver(x = granges_unlifted, chain = chain)
    mappings_to_discard <- which(elementNROWS(lifted_granges) != 1)
    
    lifted_granges_filtered <- lifted_granges[-mappings_to_discard]
    lifted_granges_filtered_unlisted <- unlist(lifted_granges_filtered)
    
    widths_to_discard <- which(lifted_granges_filtered_unlisted@ranges@width < 500) # sanity check: original peaks were 500 bp long. If < 500, then I discard (arbitrary threshold)
    
    lifted_granges_filtered_unlisted_clean <- lifted_granges_filtered_unlisted[-widths_to_discard]
    
    failed_or_split_fragments <- (length(granges_unlifted@ranges@width)
                                  - length(lifted_granges_filtered_unlisted@ranges@width))
    
    cleaned_fragments <- (length(lifted_granges_filtered_unlisted@ranges@width)
                                  - length(lifted_granges_filtered_unlisted_clean@ranges@width))
    
    summary_width_unlifted <- summary(granges_unlifted@ranges@width)
    summary_width_lifted <- summary(lifted_granges_filtered_unlisted@ranges@width)
    summary_width_cleaned <- summary(lifted_granges_filtered_unlisted_clean@ranges@width)
    
    print(paste0("# DISCARDED PEAKS DURING LIFTING PROCESS: ", failed_or_split_fragments))
    print(paste0("# DISCARDED PEAKS DURING CLEANING PROCESS: ", cleaned_fragments))
    cat("\n")
    
    print("PRINTING WIDTH DISTRIBUTION SUMMARY BEFORE LIFTING: ")
    print(summary_width_unlifted)
    cat("\n")
    print("PRINTING WIDTH DISTRIBUTION SUMMARY AFTER LIFTING: ")
    print(summary_width_lifted)
    cat("\n")
    print("PRINTING WIDTH DISTRIBUTION SUMMARY AFTER CLEANING: ")
    print(summary_width_cleaned)
    
    df <- as.data.frame(lifted_granges_filtered_unlisted_clean)
    df <- df[,c(1,2,3)]
    
    write.table(x = df, 
                file = outfile_name, 
                sep="\t", 
                quote=FALSE, row.names=FALSE, col.names=FALSE)
    
  }
  
  options(scipen = 0) # set option for scientific notation to default
  
  cat("\n")
  print("Done!")
}
