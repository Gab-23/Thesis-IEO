library("variancePartition")
library("BiocParallel")
library("limma")
library("edgeR")

setwd("/home/ieo7429/Desktop/THESIS_GAB/bulk_ATAC/GSE231481/peaks/")

metadata <- read.table(file = "../metadata.txt", header = T, tryLogical = F, sep = "\t")
count_matrix <- read.table("../count_matrices/backbone_0.1Mbp_count_matrix")

dgList <- DGEList(counts=count_matrix, 
                  genes=rownames(count_matrix))

countsPerMillion <- cpm(dgList)

countCheck <- countsPerMillion > 50 # 50 for 100 kbp

keep <- which(rowSums(countCheck) >= 6)
dgList <- dgList[keep,]

dgList <- calcNormFactors(dgList, method="TMM")


plotMDS(dgList, pch = ifelse(metadata$runID == 519, 21, 22), labels = metadata$sampleID)

sampleType<- rep("N", ncol(dgList))
sampleType[grep("T", colnames(dgList))] <- "T"

patientIDs <- unlist(lapply(X = colnames(dgList), FUN = function(x){
  parts <- strsplit(x, split = "_")[[1]]
  patient_id <- parts[1]
  return(patient_id)
}))

patientIDs <- factor(patientIDs)
sampleType <- factor(sampleType, levels = c("N", "T"))

designMat <- model.matrix(~ 0 + patientIDs + sampleType)
designMat

dgList <- estimateGLMCommonDisp(dgList, design=designMat)

fit <- glmFit(dgList, designMat)
lrt <- glmLRT(fit, coef="sampleTypeT")

edgeR_result <- topTags(lrt)

plotMD(lrt)
abline(h = c(-1, 0, 1), col = "blue", lty = 2)

da_output <- lrt$table
da_output$fdr <- p.adjust(da_output$PValue, method = "BH")
da_output$bool_diff_acc <- ifelse(da_output$fdr < 0.05, 1, 0)

coords <- lapply(X = rownames(da_output), 
       FUN = function(x){
        strsplit(x = x, split = "[:-]")[[1]]
         })

coords <- data.frame(do.call(rbind, coords))
colnames(coords) <- c("chr", "start", "end")
rownames(coords) <- rownames(da_output)

da_output <- cbind(coords, da_output)

write.table(x = da_output, file = "DA_output_backbone_0.1Mbp.tsv", sep = "\t", col.names = T)

