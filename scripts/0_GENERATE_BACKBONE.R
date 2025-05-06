source("/home/ieo7429/Desktop/THESIS_GAB/scripts/HELPER_FUNCTIONS.R") # load helpers
import_libraries() # import libraries

load_variables(genome_assembly = BSgenome.Hsapiens.UCSC.hg19, # load genome info
               window_size = 100000)

windows <- define_windows(window_size = WINDOW_SIZE, 
                          chrom_sizes_canon = chrom_sizes_canon)

# windows <- GenomicRanges::tileGenome(WINDOW_SIZE, 
#                                       chrom_size_canon, cut.tile.in.last.chrom = TRUE)

backbone <- generate_backbone(windows)

outfile_name <- paste0("/home/ieo7429/Desktop/THESIS_GAB/tables/backbone_", WINDOW_SIZE_STR, ".tsv")

write.table(x = backbone, 
            file = outfile_name, 
            row.names = F, col.names = T, 
            sep = "\t")

