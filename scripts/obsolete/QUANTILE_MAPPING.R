################### QUANTILE MAPPING #########################

# RATIONALE: Leverage a classification tool to build random samples of bins in cells
# treat them as random sample and apply Quantile Mapping approaches to correct depleted/amplified regions

library(epiAneufinder)
library(qmap)

epiAneufinder(input="./HT137B1-S1H7/HT137B1-S1H7_fragments_filtered.tsv.gz",
              outdir="../../samples/BRCA_NON_BASAL/HT137B1-S1H7/", 
              blacklist="../../tables/hg38-blacklist.v2.bed",
              windowSize=1e6, 
              genome="BSgenome.Hsapiens.UCSC.hg38",
              exclude=c('chrY','chrM'), reuse.existing=FALSE,
              title_karyo="Karyogram of sample data", 
              ncores=30, minFrags=1000,
              minsizeCNV=0, k=4, plotKaryo=FALSE)

# outputs a table with classes and one with normalized gc corrected counts.

obj <- readRDS("./HT088B1-S1H2/BRCA_HT088B1-S1H2.rds") # load the Seurat Object
non_tumor_cells <- obj@meta.data[obj$cell_type != "Tumor",]$Original_barcode
# process and obtain the resized count matrix
resized_count_matrix_filtered <- as.matrix(resized_count_matrix_filtered) # convert to dense matrix
colnames(resized_count_matrix_filtered) <- retrieve_original_barcodes(obj) # set the right colnames

epianeufinder_count_matrix <- readRDS("./HT088B1-S1H2/epiAneufinder_results/count_summary.rds") # take the raw counts
epianeufinder_count_matrix <- epianeufinder_count_matrix@assays@data@listData$counts

epianeufinder_table <- read.table("./HT088B1-S1H2/epiAneufinder_results/results_table.tsv") # load the class table
rownames(epianeufinder_table) <- paste0(epianeufinder_table$seq, "-",  # set the correct rownames
                                        epianeufinder_table$start, "-", 
                                        epianeufinder_table$end)
epianeufinder_table <- epianeufinder_table[,-c(1,2,3)] # remove coordinates
colnames(epianeufinder_table) <- gsub("cell-", "", gsub("\\.", "-", colnames(epianeufinder_table)) ) # set correct colnames
epianeufinder_table[,non_tumor_cells] <- 1

n_cells <- dim(epianeufinder_count_matrix)[2] # use same approach defined in source code to discard low quality bins
blacklist_thr <- 0.85
bins_to_discard <- which((rowSums(epianeufinder_count_matrix == 0)) > (n_cells * blacklist_thr))
epianeufinder_count_matrix <- epianeufinder_count_matrix[-bins_to_discard,] # remove those bins (in this case 2)
rownames(epianeufinder_count_matrix) <- rownames(epianeufinder_table) # assign rownames to counts

# now the Seurat matrix needs to be subsetted for the bins
resized_count_matrix_filtered <- resized_count_matrix_filtered[rownames(resized_count_matrix_filtered) %in% 
                                                                 rownames(epianeufinder_count_matrix),]

# the two epianeufinder matrices have same rows but need to be subsetted for the columns

epianeufinder_table <- epianeufinder_table[, colnames(epianeufinder_table) %in% colnames(resized_count_matrix_filtered)]
epianeufinder_count_matrix <- epianeufinder_count_matrix[, colnames(epianeufinder_count_matrix) %in% colnames(resized_count_matrix_filtered)]

epianeufinder_table <- epianeufinder_table[,colnames(resized_count_matrix_filtered)]
epianeufinder_count_matrix <- epianeufinder_count_matrix[,colnames(resized_count_matrix_filtered)]

# check for the matrix dimensions
dim(epianeufinder_table)
dim(epianeufinder_count_matrix)
dim(resized_count_matrix_filtered)

# inspect any differences --> only 60 cases, with off by 1. Why?
differences <- which(resized_count_matrix_filtered != epianeufinder_count_matrix, arr.ind = TRUE)
resized_count_matrix_filtered[differences]
epianeufinder_count_matrix[differences]
cor.test(resized_count_matrix_filtered, epianeufinder_count_matrix)

zeros_mask <- epianeufinder_table == 0 # get the mask
zeros_idx <- which(zeros_mask, arr.ind = TRUE) # get the indexes

ones_mask <- epianeufinder_table == 1
ones_idx <- which(ones_mask, arr.ind = TRUE)

twos_mask <- epianeufinder_table == 2
twos_idx <- which(twos_mask, arr.ind = TRUE)

epianeufinder_count_matrix <- as.matrix(epianeufinder_count_matrix)
epianeufinder_table <- as.matrix(epianeufinder_table)

zeros <- epianeufinder_count_matrix[zeros_mask] # extract the elements
ones <- epianeufinder_count_matrix[ones_mask]
twos <- epianeufinder_count_matrix[twos_mask]

zeros_str <- rep('0', length(zeros)) # set labels for plotting
ones_str <- rep('1', length(ones))
twos_str <- rep('2', length(twos))

df <- data.frame(labels <- c(zeros_str, ones_str, twos_str), # build the df
                 val = c(zeros, ones, twos), stringsAsFactors = TRUE)

p_boxplot <- ggplot(df, aes(x = labels, y = val)) + # first boxplot
  geom_boxplot(fill = c("purple", "green", "red")) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", aes(group = 1)) +
  theme_minimal() +
  ggtitle("Boxplots of Vectors with Different Sizes") +
  theme(legend.position = "none") +
  ylim(0,30)

quant_zeros <- quantile(zeros, seq(0,1,0.05)) # create quantiles per class
quant_ones <- quantile(ones, seq(0,1,0.05))
quant_twos <- quantile(twos, seq(0,1,0.05))

quant_df <- data.frame(quant_zeros = quant_zeros, # build the df
                       quant_ones = quant_ones,
                       quant_twos = quant_twos,
                       quantile = seq(0,1,0.05))

p_lineplot <- ggplot(quant_df, aes(x = quantile)) + # first lineplot
  geom_line(aes(y = quant_zeros, color = "quant_zeros"), size = 1, alpha = 0.5) +
  geom_line(aes(y = quant_ones, color = "quant_ones"), size = 1) +
  geom_line(aes(y = quant_twos, color = "quant_twos"), size = 1, alpha = 0.5) +
  labs(title = "Quantile Plot for Different Groups",
       x = "Quantile",
       y = "Value",
       color = "Group") +
  theme_minimal() +
  scale_color_manual(values = c("quant_zeros" = "red", 
                                "quant_ones" = "blue", 
                                "quant_twos" = "green"))

# Quantile Mapping to correct samples of an observed distribution to one of a target distribution

fit_zeros <- fitQmap(obs = ones, mod = zeros, method = "RQUANT") # fit the model
zeros_mapped <- doQmap(x = zeros, fit_zeros) # generate the corrected values

fit_twos <- fitQmap(obs = ones, mod = twos, method = "RQUANT")
twos_mapped <- doQmap(x = twos, fit_twos)

df_mapped <- data.frame(labels <- c(zeros_str, ones_str, twos_str), # build the df
                        val = c(zeros_mapped, ones, twos_mapped), stringsAsFactors = TRUE)

p_boxplot_mapped <- ggplot(df_mapped, aes(x = labels, y = val)) + # second boxplot
  geom_boxplot(fill = c("purple", "green", "red")) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", aes(group = 1)) +
  theme_minimal() +
  ggtitle("Boxplots of Vectors with Different Sizes") +
  theme(legend.position = "none") +
  ylim(0,30)

quant_zeros_mapped <- quantile(zeros_mapped, seq(0,1,0.05)) # generate corresponding quantiles
quant_twos_mapped <- quantile(twos_mapped, seq(0,1,0.05))

quant_df_mapped <- data.frame(quant_zeros_mapped = quant_zeros_mapped, # build the df
                              quant_ones = quant_ones,
                              quant_twos_mapped = quant_twos_mapped,
                              quantile = seq(0,1,0.05))

p_lineplot_mapped <- ggplot(quant_df, aes(x = quantile)) + # second lineplot
  geom_line(aes(y = quant_ones, color = "quant_ones"), size = 1) +
  geom_line(aes(y = quant_zeros_mapped, color = "quant_zeros_mapped"), size = 1, alpha = 0.5) +
  geom_line(aes(y = quant_twos_mapped, color = "quant_twos_mapped"), size = 1, alpha = 0.5) +
  labs(title = "Quantile Plot for Different Groups",
       x = "Quantile",
       y = "Value",
       color = "Group") +
  theme_minimal() +
  scale_color_manual(values = c("quant_ones" = "blue",
                                "quant_zeros_mapped" = "red",
                                "quant_twos_mapped" = "green"))

png("/home/ieo7429/Desktop/THESIS_GAB/plots/Boxplot_QM_CORRECTION.png", width = 5000, height = 7500, res = 800)
gridExtra::grid.arrange(p_boxplot, p_boxplot_mapped) # plot boxplots together
dev.off()

# retrieve summary stats
summary(zeros)
summary(zeros_mapped)
summary(ones)
summary(twos)
summary(twos_mapped)

png("/home/ieo7429/Desktop/THESIS_GAB/plots/Lineplot_QM_CORRECTION.png", width = 5000, height = 7500, res = 800)
gridExtra::grid.arrange(p_lineplot, p_lineplot_mapped) # plot lineplots together
dev.off()

png("/home/ieo7429/Desktop/THESIS_GAB/plots/HIST_QM_CORRECTION.png", width = 5000, height = 7500, res = 800)
par(mfrow = c(3,2)) # plot histograms together
hist(zeros, breaks = 150, xlim = c(0,150))
hist(zeros_mapped, breaks = 150, xlim = c(0,150))
hist(ones, breaks = 150, xlim = c(0,150))
hist(ones, breaks = 150, xlim = c(0,150))
hist(twos, breaks = 150, xlim = c(0,150))
hist(twos_mapped, breaks = 150, xlim = c(0,150))
dev.off()
par(mfrow = c(1,1))

png("/home/ieo7429/Desktop/THESIS_GAB/plots/ECDF_QM_CORRECTION.png", width = 5000, height = 7500, res = 800)
par(mfrow = c(3,2)) # plot ECDFs together
plot(ecdf(zeros))
plot(ecdf(zeros_mapped))
plot(ecdf(ones))
plot(ecdf(ones))
plot(ecdf(twos))
plot(ecdf(twos_mapped))
dev.off()
par(mfrow = c(1,1))

CNA_norm_matrix <- epianeufinder_count_matrix # copy the matrix and assign normalized values
CNA_norm_matrix[ones_idx] <- ones 
CNA_norm_matrix[zeros_idx] <- zeros_mapped
CNA_norm_matrix[twos_idx] <- twos_mapped

CNA_norm_matrix <- as(CNA_norm_matrix, "CsparseMatrix")

if (all(colnames(CNA_norm_matrix) == obj$Original_barcode)) {
  colnames(CNA_norm_matrix) <- names(obj$Original_barcode)
  print("Done!")
}

new_assay_name_norm <- paste0("bin_level_counts_","norm")
obj[[new_assay_name_norm]] <- CreateAssayObject(counts = CNA_norm_matrix)

