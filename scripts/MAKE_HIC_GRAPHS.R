hic_data_path <- "/home/ieo7429/Downloads/HiC/Hi-C_MCF7_MCF10A_processed_HiCfiles/Heatmaps/allxalll/1Mb_resolution/HiCStein-MCF10a-WT__hg19__genome__C-1000000-iced.matrix"
hic_data <- read.delim(file = hic_data_path, header = TRUE, row.names = 1, sep = "\t")

actual_rownames <- unlist(lapply(X = rownames(hic_data), 
                                 FUN = function(x){strsplit(x, 
                                                            split = "\\|")}[[1]][[3]]))

rownames(hic_data) <- actual_rownames
colnames(hic_data) <- rownames(hic_data)

all_nans_rows <- apply(X = hic_data, MARGIN = 1, FUN = function(x){all(is.nan(x))})
all_nans_cols <- apply(X = hic_data, MARGIN = 2, FUN = function(x){all(is.nan(x))})

hic_data_clean <- hic_data[!all_nans_rows, !all_nans_cols]
hic_data_clean <- 1 / (hic_data_clean + 0.0001)
hic_data_clean <- as.matrix(hic_data_clean)

annot_table <- read.table("/home/ieo7429/Desktop/THESIS_GAB/tables/ANNOTATED_WINDOWS_hg19.tsv", 
                          header = TRUE, sep = "\t")

rownames(annot_table) <- paste0(annot_table$seqnames, ":", 
                                annot_table$start, "-", 
                                annot_table$end)

non_nas_nodes <- annot_table[!is.na(annot_table$type),]

g <- igraph::graph_from_adjacency_matrix(hic_data_clean, mode = "undirected", weighted = TRUE, diag = FALSE)

V(g)$type <- NA; V(g)$type[match(rownames(non_nas_nodes), V(g)$name)] <- non_nas_nodes$type
