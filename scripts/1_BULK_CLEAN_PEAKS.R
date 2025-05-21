library(GenomicRanges)

setwd("/home/ieo7429/Desktop/THESIS_GAB/bulk_ATAC/")

summit_paths <- list.files(path = ".", recursive = T, pattern = "summits") # get the summits

outfile_paths <- unlist(lapply(X = summit_paths, FUN = function(x){

    parts <- strsplit(x = x, split = "\\.")[[1]]
    outfile_name <- paste0(parts[1], "_expanded.", parts[2])
  })) # get the output file names

lapply(X = seq_along(along.with = summit_paths), FUN = function(x){
  
  curr_summit <- summit_paths[x]
  curr_expanded <- outfile_paths[x]
  
  params <- c(curr_summit, "/home/ieo7429/Desktop/hg19/hg19.chrom.sizes.canon",
            5000, curr_expanded)

  system2("/home/ieo7429/Desktop/THESIS_GAB/scripts/expand_summits.sh", params)
  
}) # expand the summits

total_peaks_paths <- list.files(path = ".", recursive = T, pattern = "total") # get the total peaks paths

total_peaks_list <- lapply(X = total_peaks_paths, FUN = function(x){
  tab <- read.table(file = x, 
                    header = F, 
                    sep = "\t", 
                    col.names = c("chr", "start", "end",
                                  "name", "score"))
  
  tab$sample_name <- unlist(lapply(X = strsplit(x = tab$name, split = "_"), FUN = function(x){x[[1]]}))
  
  tab <- GRanges(tab)

  }) # conver them to GRanges

get_peaks_to_keep <- function(peaks_granges_sorted, n_workers, stringency = 0.5){ 
  
  n_of_samples_tot <- length(unique(peaks_granges_sorted$sample_name)) # number of total samples composing te consensus
  required_overlaps <- n_of_samples_tot * stringency # stringency tells us how conservative we are, default: a region is kept if present in half of the samples
  
  overlapping_ranges <- reduce(x = peaks_granges_sorted, with.revmap = TRUE) # get the overlapping peaks with the indexes
  revmap <- mcols(overlapping_ranges)$revmap # extract the indexes
  
  to_keep_outer <- unlist(parallel::mclapply(X = revmap, # iterate over the indexes
                                             FUN = function(indexes){
                                               
                                               if (length(indexes) <= 1) { 
                                                 # if no overlaps, discard
                                               } else {
                                                 
                                                 sample_names <- unique(mcols(peaks_granges_sorted[indexes])$sample_name) # take the unique sample names
                                                 n_of_samples_inner <- length(sample_names) # take the number of samples originating the overlapping peaks
                                                 if (n_of_samples_inner >= required_overlaps) { # if the number of samples is exceeding the threshold, we keep the most significant peak among the overlappign ones
                                                   qvals <- mcols(peaks_granges_sorted[indexes])$score 
                                                   to_keep_inner_idx <- indexes[which.max(qvals)]
                                                   to_keep_inner <- mcols(peaks_granges_sorted[to_keep_inner_idx])$name
                                                   return(to_keep_inner)
                                                   
                                                 } else {
                                                }
                                               }
                                              }, mc.cores = n_workers))
  
  return(to_keep_outer)
  
}

n_workers <- 50

clean_granges_df_list <- lapply(X = seq_along(along.with = total_peaks_list), FUN = function(x){
  
  peaks_granges <- total_peaks_list[[x]]
  
  ranking <- order(mcols(peaks_granges)$score, decreasing = TRUE) # sort the peaks
  peaks_granges_sorted <- peaks_granges[ranking]
  
  to_keep <- get_peaks_to_keep(peaks_granges_sorted = peaks_granges_sorted, 
                               n_workers = n_workers, 
                               stringency = 0.5)
  
  peaks_granges_clean <- peaks_granges[peaks_granges$name %in% to_keep]
  peaks_granges_clean_df <- as.data.frame(peaks_granges_clean) # convert to df

})

names <- c("NAT", "tumor")

lapply(X = seq_along(along.with = clean_granges_df_list), # save them
       FUN = function(x){
         
         filename <- paste0("/home/ieo7429/Desktop/THESIS_GAB/bulk_ATAC/GSE231481/peaks/", names[x], "_clean_peaks.bed")
         
         write.table(x = clean_granges_df_list[[x]][,c(1,2,3,6,7,8)], 
            file = filename, 
            row.names = F, col.names = F, quote = F, sep = "\t")
         
       })

# then take bedtools intersect -a -b -v
          # bedtools intersect -b -a -v
          # cat the two files ---> XOR peak set
# which is by definition a set of non overlapping, 10 kbp centered peaks

