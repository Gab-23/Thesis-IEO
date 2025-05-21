library(Rsubread)

convert_bed_to_saf <- function(input_bed){ # the file is assumed to be a bed + last column for peakID or binID
  output_saf <- input_bed[,c(4,1,2,3)]
  output_saf$Strand <- "*"
  return(output_saf)
}

input_bed_path <- "total_peaks_clean.bed"
bam_pattern <- "*_shifted.bam"


input_bed <- read.table(file = input_bed_path, 
                        header = F, 
                        sep = "\t", 
                        col.names = c("Chr", "Start", "End", "GeneID"))


output_saf <- convert_bed_to_saf(input_bed = input_bed)

bam_vector <- list.files(path = "../bam", 
                         recursive = T, full.names = T, 
                         pattern = bam_pattern)

count_matrix <- featureCounts(files = bam_vector, # BAM vector
                     annot.ext = output_saf, # the annotation
                     isGTFAnnotationFile = F, # annotation is BED, not GTF
                     isPairedEnd = TRUE, # ATAC is paired end
                     nthreads = 10) # how many threads
