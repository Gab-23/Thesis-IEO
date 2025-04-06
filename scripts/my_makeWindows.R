my_makeWindows <- function (genome, blacklist, windowSize, exclude = NULL) {
  
  print("HI GABRIELE! YOU ARE RUNNING THE MODIFIED VERSION OF THIS FUNCTION")
  
  # load the helper functions
  source("/home/ieo7429/Desktop/THESIS_GAB/scripts/HELPER_FUNCTIONS.R")
  print("Loaded helper function")
  
  import_libraries()
  print("Loaded helper libraries")
  
  genome <- getFromNamespace(genome, ns = genome)
  
  chrom_sizes <- seqlengths(genome) # get chromosome sizes
  chrom_sizes_canon <- chrom_sizes[1:22] # subset to canon chromosomes
  print("Defined chromosome sizes")
  
  windows <- define_windows(window_size = windowSize, 
                            chrom_sizes_canon = chrom_sizes_canon)
  print("Defined windows")
  
  windows <- GenomeInfoDb::keepStandardChromosomes(windows, 
                                                   pruning.mode = "coarse")
  windows <- dropSeqlevels(windows, exclude, pruning.mode = "coarse")
  mcols(windows)$wSeq <- as.character(seqnames(windows))
  mcols(windows)$wStart <- BiocGenerics::start(windows)
  mcols(windows)$wEnd <- BiocGenerics::end(windows)
  message("Subtracting Blacklist...")
  overlaps <- findOverlaps(windows, blacklist)
  idx <- setdiff(1:length(windows), S4Vectors::queryHits(overlaps))
  windowsBL <- windows[idx]
  names(windowsBL) <- paste0("w", seq_along(windowsBL))
  mcols(windowsBL)$name <- names(windowsBL)
  message("Adding Nucleotide Information...")
  windowSplit <- split(windowsBL, as.character(seqnames(windowsBL)))
  windowNuc <- lapply(seq_along(windowSplit), function(x) {
    message(sprintf("%s of %s", x, length(windowSplit)))
    chrSeq <- Biostrings::getSeq(genome, names(windowSplit)[x])
    grx <- windowSplit[[x]]
    aFreq <- Biostrings::alphabetFrequency(Biostrings::Views(chrSeq, 
                                                             ranges(grx)))
    mcols(grx)$GC <- rowSums(aFreq[, c("G", "C"), drop = FALSE])/rowSums(aFreq)
    mcols(grx)$AT <- rowSums(aFreq[, c("A", "T"), drop = FALSE])/rowSums(aFreq)
    return(grx)
  }) %>% GRangesList %>% unlist %>% sortSeqlevels %>% sort
  windowNuc$N <- 1 - (windowNuc$GC + windowNuc$AT)
  print("Finished making windows successfully")
  windowNuc
}