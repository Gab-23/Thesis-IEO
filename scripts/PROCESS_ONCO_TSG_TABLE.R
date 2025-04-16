source("/home/ieo7429/Desktop/THESIS_GAB/scripts/HELPER_FUNCTIONS.R")
import_libraries()

unprocessed_table <- read.delim("/home/ieo7429/Desktop/THESIS_GAB/tables/cancerGeneList.tsv", 
                                header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)

unprocessed_table <- unprocessed_table[,c(1,7,8)]

to_discard <- which((unprocessed_table$Is.Oncogene == "Yes") & (unprocessed_table$Is.Tumor.Suppressor.Gene == "Yes"))

unprocessed_table_clean <- unprocessed_table[-to_discard,]

unprocessed_table_clean$type <- ifelse(unprocessed_table_clean$Is.Oncogene == "Yes", "OG", "TSG")
unprocessed_table_clean$Is.Oncogene <- NULL
unprocessed_table_clean$Is.Tumor.Suppressor.Gene <- NULL

ensembl37 <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", 
                              dataset = "hsapiens_gene_ensembl", 
                              host = "https://grch37.ensembl.org")


genes <- unprocessed_table_clean$Hugo.Symbol

coords <- biomaRt::getBM(attributes = c("hgnc_symbol", "chromosome_name", 
                                        "start_position", "end_position", 
                                        "strand", "gene_biotype"),
                         filters = "hgnc_symbol",
                         values = genes,
                         mart = ensembl37)

processed_data <- merge(x = unprocessed_table_clean, y = coords, by.x = "Hugo.Symbol", by.y = "hgnc_symbol")
processed_data <- processed_data[,c(1,2,3,4,5)]

accepted_chr <- c(as.character(seq(1:22)),"X","Y")

processed_data_clean <- processed_data[processed_data$chromosome_name %in% accepted_chr,]
processed_data_clean$chromosome_name <- paste0("chr",processed_data_clean$chromosome_name)
colnames(processed_data_clean) <- c("hugo_symbol", "type", "chr", "start", "end")

write.table(x = processed_data_clean, file = "/home/ieo7429/Desktop/THESIS_GAB/tables/ONCO_TSG_TABLE_hg19.tsv", 
            sep = "\t", col.names = TRUE, row.names = FALSE)


chrom_sizes_canon <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:22]
windows <- define_windows(1000000, chrom_sizes_canon)

processed_data_granges <- GRanges(seqnames = processed_data_clean$chr, 
                                             ranges = IRanges(processed_data_clean$start, 
                                                              processed_data_clean$end),
                                             name = processed_data_clean$hugo_symbol,
                                             type = processed_data_clean$type)

hits <- findOverlaps(query = processed_data_granges, subject = windows)

overlap_df <- data.frame(
  window_id = subjectHits(hits),
  gene_name = processed_data_granges$name[queryHits(hits)],
  gene_type = processed_data_granges$type[queryHits(hits)]
)

summary_df <- overlap_df %>% 
  group_by(window_id) %>%
  summarise(
    genes = paste(gene_name, collapse = ","),
    types = paste(gene_type, collapse = ","),
    n_of_genes = n(),
    n_of_TSGs = table(gene_type)["TSG"],
    n_of_OGs = table(gene_type)["OG"]
  )

annotated_windows_df <- as.data.frame(windows)
annotated_windows_df$window_id <- seq_len(nrow(annotated_windows_df))

final_df <- left_join(annotated_windows_df, summary_df, by = "window_id")
final_df[is.na(final_df)] <- 0

final_df$genes <- ifelse(final_df$genes == 0, NA, final_df$genes)
final_df$types <- ifelse(final_df$types == 0, NA, final_df$types)


final_df$width <- NULL
final_df$strand <- NULL
final_df$window_id <- NULL

colnames(final_df) <- c("chr", "start", "end", 
                        "genes", "types", 
                        "n_of_genes", "n_of_TSGs", "n_of_OGs")

write.table(final_df, 
            file = "/home/ieo7429/Desktop/THESIS_GAB/tables/ANNOTATED_WINDOWS_AGGREGATED_hg19.tsv", 
            sep = "\t", col.names = TRUE, row.names = FALSE)

