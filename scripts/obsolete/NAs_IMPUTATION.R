################ COHORT PROCESSING WITH EPIANEUFINDER (NAs IMPUTING) ############

epianeufinder_count_matrix <- readRDS("HT088B1-S1H2/epiAneufinder_results/count_summary.rds") # + process it
epianeufinder_table <- read.table("HT088B1-S1H2/epiAneufinder_results/results_table.tsv") # + process it

CNA_mask <- epianeufinder_table != 1 # compute a binary mask for CNA events
table_CNA_idx <- which(CNA_mask, arr.ind = T) # take the CNA indexes

epianeufinder_count_matrix[table_CNA_idx] <- NA # set them as NAs

range <- seq_along(1:nrow(epianeufinder_count_matrix)) # take the iteration range

for (idx in range) { # iterate each row
  row <- epianeufinder_count_matrix[idx, ] # take the i-th row
  row_no_NAs <- row[!is.na(row)] # take the non NAs from the row
  
  cells_selection <- sample(x = row_no_NAs[tumor_cells %in% names(row_no_NAs)], 
                            size = length(non_tumor_cells))
  
  median_row <- median(cells_selection) # compute the median without NAs
  row[which(is.na(row))] <- median_row # assign the median to NAs
  
  epianeufinder_count_matrix[idx, ] <- row  # update the matrix
}