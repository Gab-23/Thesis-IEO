library(GenomicRanges)
library(zoo)

# load data
setwd("/home/ieo7429/Desktop/RepliSeq_pipeline_processing/covtracks/")

window_size <- "1Mbp"

k <- ifelse(window_size == "1Mbp",1,10)

backbone_path <- paste0("../../THESIS_GAB/tables/backbone_",window_size,".tsv")
E_bg_path <- paste0("SRR28552523_earlyS/SRR28552523_",window_size,"_covtrack")
L_bg_path <- paste0("SRR28552522_lateS/SRR28552522_",window_size,"_covtrack")
G1_bg_path <- paste0("SRR28552524_G1/SRR28552524_",window_size,"_covtrack")

backbone <- read.table(file = backbone_path, header = T, sep = "\t")

E <- read.table(file = E_bg_path, 
                header = F, 
                sep = "\t", 
                col.names = c("chr","start","end","binID","counts"))

L <- read.table(file = L_bg_path, 
                header = F, sep = "\t", 
                col.names = c("chr","start","end","binID","counts"))

G1 <- read.table(file = G1_bg_path, 
                 header = F, 
                 sep = "\t", 
                 col.names = c("chr","start","end","binID","counts"))

E$binID <- NULL; L$binID <- NULL; G1$binID <- NULL

# filter for 0 counts and common regions

E <- E[E$counts > 0, ]; L <- L[L$counts > 0, ]; G1 <- G1[G1$counts > 0, ] 

common_regions <- Reduce(intersect, list(rownames(E), rownames(L), rownames(G1)))

E_filt <- E[common_regions,]; L_filt <- L[common_regions,]; G1_filt <- G1[common_regions,]

total_reads_E <- sum(E_filt$counts)
total_reads_L <- sum(L_filt$counts)
total_reads_G1 <- sum(G1_filt$counts)

E_filt$bin_length <- E_filt$end - E_filt$start
L_filt$bin_length <- L_filt$end - L_filt$start
G1_filt$bin_length <- G1_filt$end - G1_filt$start

get_rpkm <- function(df, total_reads){
  df$rpkm <- df$counts / 
              ((df$bin_length / 1000) * 
                 (total_reads / 1000000))
  return(df)
}

E_filt <- get_rpkm(E_filt, total_reads_E)
L_filt <- get_rpkm(L_filt, total_reads_L)
G1_filt <- get_rpkm(G1_filt, total_reads_G1)

E_filt$E_G1_norm <- E_filt$rpkm / G1_filt$rpkm
L_filt$L_G1_norm <- L_filt$rpkm / G1_filt$rpkm
EL_ratio <- E_filt$E_G1_norm / L_filt$L_G1_norm
smoothed <- rollmean(EL_ratio, k = k, fill = NA, align = "center")

log2EL_ratio <- log2(EL_ratio)
log2smoothed <- log2(smoothed)

repliseq_df <- data.frame(chr = E_filt$chr,
                          start = E_filt$start,
                          end = E_filt$end,
                          E_rpkm = E_filt$rpkm,
                          L_rpkm = L_filt$rpkm,
                          G1_rpkm = G1_filt$rpkm,
                          E_G1_norm = E_filt$E_G1_norm,
                          L_G1_norm = L_filt$L_G1_norm,
                          EL_ratio = EL_ratio,
                          smoothed_ratio = smoothed,
                          log2EL_ratio = log2EL_ratio,
                          log2smoothed = log2smoothed,
                          replication_class = ifelse(log2smoothed >= 0, "E", "L"))

sub_df <- repliseq_df[repliseq_df$chr == "chr17", ]
x <- seq_along(sub_df$log2smoothed)
y <- sub_df$log2smoothed

plot(x, y, type = "n", 
     main = "Smoothed replication timing (1Mbp windows) (chr17)", 
     ylab = "Smoothed log2FC", 
     xlab = "Genomic bins")

for (i in 1:(length(x) - 1)) {
  
  x_poly <- c(x[i], x[i+1], x[i+1], x[i])
  y_poly <- c(0, 0, y[i+1], y[i])
  
  if (y[i] >= 0 && y[i+1] >= 0) {
    col_fill <- "lightblue"
  } else if (y[i] < 0 && y[i+1] < 0) {
    col_fill <- "lightcoral"
  } else {
    
    slope <- (y[i+1] - y[i]) / (x[i+1] - x[i])
    x0 <- x[i] - y[i] / slope
    
    polygon(c(x[i], x0, x[i]), c(y[i], 0, 0),
            col = ifelse(y[i] > 0, "lightblue", "lightcoral"), border = NA)
    
    polygon(c(x0, x[i+1], x[i+1]), c(0, y[i+1], 0),
            col = ifelse(y[i+1] > 0, "lightblue", "lightcoral"), border = NA)
    next
  }
  polygon(x_poly, y_poly, col = col_fill, border = NA)
}
lines(x, y, col = "black", lwd = 1.5)
abline(h = 0, col = "black", lty = 2)

repliseq_df$start <- repliseq_df$start + 1
repliseq_df <- merge(x = backbone, y = repliseq_df, by = c("chr","start","end"))
repliseq_df$chr <- NULL; repliseq_df$end <- NULL; repliseq_df$start <- NULL
write.table(x = repliseq_df, 
            file = paste0("../repliseq_postproc_data_",window_size,".tsv"), 
            sep = "\t", 
            col.names = T)



