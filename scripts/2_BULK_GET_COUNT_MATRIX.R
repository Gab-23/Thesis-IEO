library(Rsubread)

convert_bed_to_saf <- function(input_bed, mode){ # the file is assumed to be a bed + last column for peakID or binID
  
  if (mode == "peaks") {
    input_bed$score <- NULL; input_bed$sample_name <- NULL  
  }
  
  input_bed$GeneID <- paste0(input_bed$Chr,":",input_bed$Start,"-",input_bed$End)
  
  output_saf <- input_bed[,c(4,1,2,3)]
  output_saf$Strand <- "*"
  
  order <- order(gsub(pattern = "chr", 
                      replacement = "", 
                      x = output_saf$Chr), 
                 output_saf$Start)
  
  output_saf <- output_saf[order, ]
  
  return(output_saf)
}

setwd("GSE231481/peaks/")
input_bed_path <- "../../../tables/backbone_1Mbp.tsv"
bam_pattern <- "*_shifted.bam"
mode <- "backbone_1Mbp"


if (grepl(pattern = "peaks", x = mode)) {
  
  input_bed <- read.table(file = input_bed_path, 
                          header = F, 
                          sep = "\t", 
                          col.names = c("Chr", 
                                        "Start", 
                                        "End", 
                                        "GeneID",
                                        "score",
                                        "sample_name"))
  
} else {
  
  input_bed <- read.table(file = input_bed_path, 
                        header = T, 
                        sep = "\t", 
                        col.names = c("Chr", 
                                      "Start", 
                                      "End", 
                                      "GeneID"))
}

output_saf <- convert_bed_to_saf(input_bed = input_bed, mode = mode)

bam_vector <- list.files(path = "../bam", 
                         recursive = T, full.names = T, 
                         pattern = bam_pattern)

count_matrix_output <- featureCounts(files = bam_vector, # BAM vector
                     annot.ext = output_saf, # the annotation
                     isGTFAnnotationFile = FALSE, # annotation is BED, not GTF
                     isPairedEnd = TRUE, # ATAC is paired end
                     nthreads = 10) # how many threads


count_matrix <- count_matrix_output[["counts"]]
metadata <- read.table(file = "../metadata.txt", header = T, sep = "\t")

metadata$sampleID <- unlist(lapply(X = metadata$sampleID, FUN = function(x){
  parts <- strsplit(x = x, split = "_")[[1]]
  new_colname <- paste0(parts[1],"_",parts[2])
}))

colnames(count_matrix) <- unlist(lapply(X = colnames(count_matrix), FUN = function(x){
  sample_name <- strsplit(x, split = "_")[[1]][1]
  sample_id <- metadata[metadata$sraID == sample_name, "sampleID"]
  return(sample_id)
}))

write.table(x = count_matrix, file = paste0("../count_matrices/", mode, "_count_matrix"), sep = "\t", col.names = TRUE)
