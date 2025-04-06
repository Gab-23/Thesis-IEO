TSG_ONCOGENES_DF <- read.table("../../TABLE_TSG_ONCO_hg38", header = F, sep = "\t")
TSG_ONCOGENES_DF$V3 <- c(rep("TSG", 21),rep("ONCO",(nrow(TSG_ONCOGENES_DF)-21)))
TSG_ONCOGENES_DF <- TSG_ONCOGENES_DF[-c(1,22),]
TSG_ONCOGENES_DF <- cbind(TSG_ONCOGENES_DF$V1, 
                          do.call(rbind, strsplit(TSG_ONCOGENES_DF$V2, split = "[:-]")), 
                          TSG_ONCOGENES_DF$V3)
TSG_ONCOGENES_DF <- as.data.frame(TSG_ONCOGENES_DF)
colnames(TSG_ONCOGENES_DF) <- c("gene_name", "chr", "start", "end", "type")

TSG_ONCOGENES_DF$start <- as.numeric(gsub(x = TSG_ONCOGENES_DF$start, 
                                          pattern = ",", 
                                          replacement = "")) 

TSG_ONCOGENES_DF$end <- as.numeric(gsub(x = TSG_ONCOGENES_DF$end, 
                                        pattern = ",", 
                                        replacement = ""))

TSG_ONCOGENES_DF$chr <- ifelse(TSG_ONCOGENES_DF$chr == "chrx", "chrX",TSG_ONCOGENES_DF$chr)

TSG_ONCOGENES_DF <- process_DACRs_coords(TSG_ONCOGENES_DF, 
                                         chrom_cum_start)

write.table(TSG_ONCOGENES_DF, "../../TABLE_TSG_ONCO_hg38_PROCESSED", sep = "\t", row.names = F, col.names = T)
