retrieve_rownames_of_GRanges_object <- function(GRanges_object){
  names_GRanges_object <- paste(seqnames(GRanges_object), 
                                start(GRanges_object), 
                                end(GRanges_object),
                                sep = "-")
}
define_windows <- function(window_size, chrom_sizes_canon){
  
  windows <- tileGenome(seqlengths = chrom_sizes_canon,
                        tilewidth = window_size,
                        cut.last.tile.in.chrom = T)
  
  genome(windows) <- "hg38"
  names(windows) <- retrieve_rownames_of_GRanges_object(windows)
  return(windows)
}
filter_windows <- function(windows, original_peaks){
  
  n_of_overlaps <- countOverlaps(query = windows,
                                 subject = original_peaks, 
                                 type = "any")
  windows_filtered <- windows[n_of_overlaps > 0]
  return(windows_filtered)
}
retrieve_original_barcodes <- function(seurat_obj){
  sample_level_names <- colnames(seurat_obj)
  sample_names_split <- do.call(what = rbind, strsplit(sample_level_names, split = "_"))
  raw_names <- sample_names_split[,3]
  return(raw_names)
}
handle_fragments <- function(fragments_path, seurat_obj, original_barcodes){
  
  file_name_parts <- strsplit(fragments_path, "\\.")[[1]][1:2]
  fragments_bgz_path <- paste0(file_name_parts[1], ".", file_name_parts[2],".","bgz")
  fragments_idx_path <- paste0(fragments_bgz_path,".tbi")
  
  if (!(file.exists(fragments_path))) {
    stop("fragments file does not exist")
  } else if (file.exists(fragments_bgz_path)) {
    message("fragments file is already bgzipped")
  } else {bgzip(fragments_path)}
  
  if (!(file.exists(fragments_bgz_path))) {
    stop("bgzipped fragments file does not exist")
  } else if (file.exists(fragments_idx_path)) {
    message("fragments file index already exists")  
  } else {indexTabix(fragments_bgz_path, format = "bed")}
  
  fragment_obj <- CreateFragmentObject(path = fragments_bgz_path, cells = original_barcodes)
  return(fragment_obj)
}
resize_count_matrix <- function(fragment_obj, windows, cells_barcodes, sample_level_names){
  
  resized_count_matrix <- FeatureMatrix(fragments = fragment_obj,
                                        features = windows, 
                                        cells = cells_barcodes)
  
  colnames(resized_count_matrix) <- sample_level_names
  return(resized_count_matrix)
}
get_DACRs_coords <- function(markers){
  
  DACRs <- rownames(markers)
  DACRs_coords <- as.data.frame(do.call(what = rbind, sapply(X = DACRs, FUN = strsplit, split = "-")))
  colnames(DACRs_coords) <- c("chr","start","end")
  DACRs_coords$chr <- factor(DACRs_coords$chr, levels = paste0("chr",c(as.character(1:22), "X","Y")))
  DACRs_coords$start <- as.integer(DACRs_coords$start)
  DACRs_coords$end <- as.integer(DACRs_coords$end)
  DACRs_coords$log2FC <- markers$avg_log2FC
  DACRs_coords$pval_adj <- markers$p_val_adj
  DACRs_coords <- DACRs_coords %>% arrange(chr, start)
  
  return(DACRs_coords)
}
process_DACRs_coords <- function(DACRs_coords, chrom_cum_start){
  DACRs_coords <- DACRs_coords %>%
    mutate(
      cum_start = start + chrom_cum_start[chr],
      cum_end = end + chrom_cum_start[chr]
    )
  
  DACRs_coords$absolute_midpoints <- (DACRs_coords$cum_start + DACRs_coords$cum_end) / 2 # define cumulative middle point
  return(DACRs_coords)
}
plot_log2FC <- function(DACRs_df, thr = 0.05) {
  
  chr_midpoints <- DACRs_df %>%
    group_by(chr) %>%
    summarise(midpoint = mean(absolute_midpoints, na.rm = TRUE))
  
  chr_labels <- unique(DACRs_df$chr)
  if (length(chr_labels) > 1) {
    title_chr <- "[WHOLE GENOME]"
  } else {
    title_chr <- paste0("[", chr_labels, "]")
  }
  
  p <- ggplot(DACRs_df, aes(x = absolute_midpoints, y = log2FC)) +
    geom_point(color = ifelse(DACRs_df$pval_adj < thr, "red", "blue"), 
               alpha = ifelse(DACRs_df$pval_adj < thr, 1, 0.5), 
               size = ifelse(DACRs_df$pval_adj < thr, 3, 1.5)) +
    theme_minimal() +
    labs(x = "", 
         y = "log2FC", 
         title = paste("Manhattan Plot", title_chr), 
         subtitle = "Genomic Coordinates vs Log2FC") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 15), 
      axis.text.y = element_text(size = 15)
    ) +
    scale_x_continuous(
      breaks = chr_midpoints$midpoint,
      labels = chr_midpoints$chr
    ) +
    ylim(-(max(abs(DACRs_df$log2FC))),max(abs(DACRs_df$log2FC)))
  
  return(p)
}

plot_log2FC_collapsed <- function(DACRs_df, thr = 0.05) {
  
  chr_midpoints <- DACRs_df %>%
    group_by(chr) %>%
    summarise(midpoint = mean(absolute_midpoints, na.rm = TRUE))
  
  chr_labels <- unique(DACRs_df$chr)
  if (length(chr_labels) > 1) {
    title_chr <- "[WHOLE GENOME]"
  } else {
    title_chr <- paste0("[", chr_labels, "]")
  }
  
  p <- ggplot(DACRs_df, aes(x = absolute_midpoints, y = mean_log2FC, 
                            color = mean_log2FC,
                            size = n_of_significant_per_bin)) +
    geom_point() +
    scale_color_gradient(low = "blue", high = "red") +
    scale_size_continuous(range = c(0.5,6)) +
    theme_minimal() +
    labs(x = "", 
         y = "mean log2FC", 
         title = paste("Manhattan Plot", title_chr), 
         subtitle = "Genomic Coordinates vs mean Log2FC [Tumor VS CNCs]") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.title.y = element_text(size = 15),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 15), 
      axis.text.y = element_text(size = 15)
    ) +
    scale_x_continuous(
      breaks = chr_midpoints$midpoint,
      labels = chr_midpoints$chr
    ) +
    ylim(-(max(abs(DACRs_df$mean_log2FC))),
         max(abs(DACRs_df$mean_log2FC)))
  
  return(p)
}
