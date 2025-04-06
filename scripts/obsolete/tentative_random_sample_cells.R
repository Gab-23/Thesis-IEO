B = 5
global_seed <- sample(1:1e6, 1)
local_seeds <- sample(1:1e6,B)

set.seed(global_seed)

obj <- readRDS("/home/ieo7429/Desktop/THESIS_GAB/samples/BRCA_NON_BASAL/HT088B1-S1H2/BRCA_HT088B1-S1H2.rds")

cell_type_1 <- "Tumor"

lapply(X = seq(1:B), 
       FUN = function(i){
         
         set.seed(local_seeds[i])
         non_tumor_cells <- obj@meta.data[obj@meta.data$cell_type != cell_type_1,]$Original_barcode
         num_non_tumor_cells <- length(non_tumor_cells)
         sampled_tumor_cells <- sample(x = obj@meta.data[obj@meta.data$cell_type == cell_type_1,]$Original_barcode , 
                                       size = num_non_tumor_cells * 1.5, 
                                       replace = F)
         
         barcodes_sample <- c(sampled_tumor_cells, non_tumor_cells)
         
         file_name <- paste0("/home/ieo7429/Desktop/THESIS_GAB/samples/barcodes_random_sample_",i)
         file <-  file(file_name, "wb")
         writeBin(paste(barcodes_sample, collapse="\n"), file)
         close(file)
         
         msg <- paste0("barcodes_random_sample_",i, " written")
         return(msg)
         })


