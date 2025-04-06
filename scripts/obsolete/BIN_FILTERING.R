################### FILTERING ########################

# RATIONALE: find those bins with frequency of NON-copy number neutral > threshold

epianeufinder_count_matrix <- readRDS("HT088B1-S1H2/epiAneufinder_results/count_summary.rds") # + process it
epianeufinder_table <- read.table("HT088B1-S1H2/epiAneufinder_results/results_table.tsv") # + process it

CNA_mask <- epianeufinder_table != 1 # generate a logical mask to identify CNA events
bin_level_frequency <- rowSums(CNA_mask) / dim(CNA_mask)[2] # get the frequency of cells having non copy number neutral status for a bin


q1 <- quantile(bin_level_frequency, 0.25) # retrieve quantiles
q3 <- quantile(bin_level_frequency, 0.75)

bin_thr <- q3 + (1.5 * (q3 - q1)) # set the threshold for the bin

bin_level_frequency_filtered <- bin_level_frequency[bin_level_frequency <= bin_thr] # filter the frequency vector
filtered_bins <- names(bin_level_frequency_filtered) # take the filtered bins

par(mfrow = c(3,2)) # plot
hist(bin_level_frequency, breaks = 100)
hist(bin_level_frequency_filtered, breaks = 100)
plot(density(bin_level_frequency))
plot(density(bin_level_frequency_filtered))
boxplot(bin_level_frequency)
boxplot(bin_level_frequency_filtered)
par(mfrow = c(1,1))

epianeufinder_table_filtered <- epianeufinder_table[filtered_bins,] # subset the table
epianeufinder_count_matrix_filtered <- epianeufinder_count_matrix[filtered_bins,] # subset the matrix
