# STARTING FROM FABIOS DATASETS, GENERATE THE TOTAL VERSION BY SUMMING AMPL AND DEL FREQUENCIES
# THE FEATURE MATRIX IS ALWAYS THE SAME BUT THE TARGETS CHANGE.
# AFTER DOING EVERYTHING I ATTACH THE TWO ADDITIONAL COLUMNS (CHROMOSOME LENGTH AND CENTROMERE LENGTH)

source("/home/ieo7429/Desktop/THESIS_GAB/scripts/HELPER_FUNCTIONS.R") # load helpers
import_libraries() # import libraries

load_variables(genome_assembly = BSgenome.Hsapiens.UCSC.hg19, # load genome info
               window_size = 100000)

backbone_name <- paste0("/home/ieo7429/Desktop/THESIS_GAB/tables/backbone_", WINDOW_SIZE_STRING, ".tsv")
backbone <- read.table(backbone_name, header = T, sep = "\t") # load backbone

wd_name <- paste0("/home/ieo7429/Desktop/THESIS_GAB/tables/fabios_tables/", WINDOW_SIZE_STRING)
setwd(wd_name)
files <- grep(pattern = WINDOW_SIZE_STRING, x = list.files("."), value = TRUE)[1:2] # probably won't need it

CANCER_TYPES <- c("BRCA")

for (file in files) { # take fabios files
  load(file)
  assign(file, ML.Tables)
  rm(ML.Tables)
}

# subset for current cohort

chr <- Length_Location_ampl_or_del_0.1Mbp_covThr_zero.RData$ampl_or_del$`Chromosome-level`; chr <- chr[chr$Type %in% CANCER_TYPES, ]
arm <- Length_Location_ampl_or_del_0.1Mbp_covThr_zero.RData$ampl_or_del$`Arm-level`; arm <- arm[arm$Type %in% CANCER_TYPES, ]
mid <- Length_ampl_or_del_0.1Mbp_covThr_zero.RData$ampl_or_del$`Mid-length`; mid <- mid[mid$Type %in% CANCER_TYPES, ]
small <- Length_ampl_or_del_0.1Mbp_covThr_zero.RData$ampl_or_del$`Small-scale`; small <- small[small$Type %in% CANCER_TYPES, ]

# make a list of data.frames

dataset_tables_different_levels <- list(chr_level = chr, 
                                        arm_level = arm,
                                        mid_level = mid, 
                                        small_level = small)

summed_scores <- dataset_tables_different_levels %>%
                 bind_rows() %>%  # rbind all rows
                 group_by(bin) %>% # group by bin ID
                 summarise(
                   ampl_score_tot = sum(ampl_score), # take the sum as the total score
                   del_score_tot = sum(del_score)) %>% # take the sum as the total score
                 ungroup()

tot <- dataset_tables_different_levels$chr_level # take one covariate dataset
tot <- merge(x = tot, y = summed_scores, by = "bin") # merge by bin
tot$ampl_score <- tot$ampl_score_tot; tot$ampl_score_tot <- NULL # reassign variable
tot$del_score <- tot$del_score_tot; tot$del_score_tot <- NULL # reassign variable

dataset_tables_different_levels$tot_level <- tot # add it to the list

chrom_information <- rCGH::hg19 # take chromosome information for hg19
chrom_information$cumlen <- NULL # remove cumulative length
chrom_information$cen_len <- (chrom_information$centromerEnd - chrom_information$centromerStart) # take centromere length (for hg18 to hg39 it is fixed 3Mbp)

centromere_table_from_fabio <- read.table("/home/ieo7429/Desktop/THESIS_GAB/tables/centomere.tsv", header = T, sep = "")
centromere_table_from_fabio$start <- NULL; centromere_table_from_fabio$end <- NULL

colnames(centromere_table_from_fabio) <- c("chr", "Centromere_Length_Fabio", "Centromere_Type", "Centromere")
colnames(chrom_information) <- c("chr", "Chromosome_Length", "Centromere_Start", "Centromere_End", "Centromere_Length")

chrom_information_updated <- merge(x = chrom_information, y = centromere_table_from_fabio, by = "chr")
chrom_information_updated$chr <- paste0("chr", chrom_information_updated$chr)

chrom_information_updated_with_backbone <- merge(x = backbone, y = chrom_information_updated, by = "chr")

chrom_information_updated_with_backbone$chr <- NULL
chrom_information_updated_with_backbone$start <- NULL
chrom_information_updated_with_backbone$end <- NULL

save(dataset_tables_different_levels, 
     file = "/home/ieo7429/Desktop/THESIS_GAB/tables/datasets_different_levels_total.RData")




