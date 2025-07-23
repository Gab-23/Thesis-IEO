plot_CNA_acc_correlation <- function(df_before, df_after, mode){
  
  if (mode == "lm") {
    formula <- as.formula(y ~ x)
  } else if (mode == "poly") {
    formula <- as.formula(y ~ poly(x, 2))
  }
  
  ampl_acc_test <- cor.test(df_before$ampl_score, df_before$logFC, method = "spearman") # perform correlation test. spearman because lack of normality for both variables
  cor_ampl_acc <- round(ampl_acc_test$estimate, 2) # take correlation value
  pval_ampl_acc <- round(ampl_acc_test$p.value, 5) # take pvalue
  
  ampl_acc <- ggplot(df_before, aes(x = ampl_score, y = logFC, colour = chr)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_smooth(method = "lm", formula = formula, se = FALSE, color = "blue") +
    annotate("text", x = max(df_before$ampl_score) * 0.7, y = max(df_before$ampl_score) * 0.9, 
             label = paste("Correlation: ", cor_ampl_acc, "\n",
                           "P-value: ", pval_ampl_acc), 
             color = "black", size = 5, fontface = "bold") +
    labs(title = "[BEFORE NORMALIZATION]",
         x = "Amplification score",
         y = "Chromatin accessibility [LogFC]") +
    theme_minimal()
  
  ampl_acc_regrout_test <- cor.test(df_after$ampl_score, df_after$logFC, method = "spearman")
  cor_ampl_acc_regrout <- round(ampl_acc_regrout_test$estimate, 2)
  pval_ampl_acc_regrout <- round(ampl_acc_regrout_test$p.value, 5)
  
  ampl_acc_regrout <- ggplot(df_after, aes(x = ampl_score, y = logFC, colour = chr)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_smooth(method = "lm", formula = formula, se = FALSE, color = "blue") +
    annotate("text", x = max(df_after$ampl_score) * 0.7, y = max(df_after$ampl_score) * 0.9, 
             label = paste("Correlation: ", cor_ampl_acc_regrout, "\n",
                           "P-value: ", pval_ampl_acc_regrout), 
             color = "black", size = 5, fontface = "bold") +
    labs(title = "[REGRESS OUT AMPLIFICATION + DELETION EFFECT]",
         x = "Amplification score",
         y = "Chromatin accessibility [LogFC]") +
    theme_minimal()
  
  del_acc_test <- cor.test(df_before$del_score, df_before$logFC, method = "spearman") # perform correlation test. spearman because lack of normality for both variables
  cor_del_acc <- round(del_acc_test$estimate, 2) # take correlation value
  pval_del_acc <- round(del_acc_test$p.value, 5) # take pvalue
  
  del_acc <- ggplot(df_before, aes(x = del_score, y = logFC, colour = chr)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_smooth(method = "lm", formula = formula, se = FALSE, color = "blue") +
    annotate("text", x = max(df_before$del_score) * 0.7, y = max(df_before$del_score) * 0.9, 
             label = paste("Correlation: ", cor_del_acc, "\n",
                           "P-value: ", pval_del_acc), 
             color = "black", size = 5, fontface = "bold") +
    labs(title = "[BEFORE NORMALIZATION]",
         x = "Deletion score",
         y = "Chromatin accessibility [LogFC]") +
    theme_minimal()
  
  del_acc_regrout_test <- cor.test(df_after$del_score, df_after$logFC, method = "spearman")
  cor_del_acc_regrout <- round(del_acc_regrout_test$estimate, 2)
  pval_del_acc_regrout <- round(del_acc_regrout_test$p.value, 5)
  
  del_acc_regrout <- ggplot(df_after, aes(x = del_score, y = logFC, colour = chr)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_smooth(method = "lm", formula = formula, se = FALSE, color = "blue") +
    annotate("text", x = max(df_after$del_score) * 0.7, y = max(df_after$del_score) * 0.9, 
             label = paste("Correlation: ", cor_del_acc_regrout, "\n",
                           "P-value: ", pval_del_acc_regrout), 
             color = "black", size = 5, fontface = "bold") +
    labs(title = "[REGRESS OUT AMPLIFICATION + DELETION EFFECT]",
         x = "Deletion score",
         y = "Chromatin accessibility [LogFC]") +
    theme_minimal()
  
  gridplot <- gridExtra::grid.arrange(ampl_acc, ampl_acc_regrout, 
                                      del_acc, del_acc_regrout,
                                      nrow = 2, ncol = 2)
  
  return(gridplot)
  
} # plot correlation between chromatin accessibility and CNA before and after out regression step

REGRESS_OUT_COHORT <- function(cancer_type, backbone_path, 
                               accessibility_table_path, 
                               ampl_del_table_path, 
                               mode = "lm"){
  
  load(ampl_del_table_path)
  ampl_del_table <- ML_dataset_with_target %>% 
    filter(Type == cancer_type) %>% 
    dplyr::select(bin, ampl_score, del_score)
  
  backbone <- read.table(backbone_path, 
                         header = T, sep = "\t")
  
  acc_table <- read.table(accessibility_table_path, 
                          header = T, sep = "\t")
  
  acc_table <- acc_table %>% dplyr::select(chr, start, end, logFC, bool_diff_acc)
  acc_table <- merge(x = backbone, y = acc_table, by = c("chr", "start", "end"))
  
  df <- merge(x = acc_table, y = ampl_del_table, by = "bin")
  
  if (mode == "lm") {
    formula <- logFC ~ ampl_score + del_score
    model_both <- summary(lm(formula = formula, data = df))  
  } else if (mode == "poly") {
    formula <- mean_log2FC ~ poly(ampl_score, degree = 3, raw = TRUE) + poly(del_score, degree = 3, raw = TRUE)
    model_both <- summary(lm(formula = formula, data = df))
  }
  
  df_both <- df
  df_both$logFC <- model_both$residuals
  
  gridplot <- plot_CNA_acc_correlation(df, df_both, mode)
  
  log2k_1 <- log2(1) 
  log2k_2 <- log2(2)
  log2k_3 <- log2(3)
  
  sign_mean_log2FC_1 <- ifelse(df_both$logFC < -log2k_1, -1, ifelse(df_both$logFC > log2k_1, 1, 0))
  sign_mean_log2FC_2 <- ifelse(df_both$logFC < -log2k_2, -1, ifelse(df_both$logFC > log2k_2, 1, 0))
  sign_mean_log2FC_3 <- ifelse(df_both$logFC < -log2k_3, -1, ifelse(df_both$logFC > log2k_3, 1, 0))
  
  df_both$ampl_score <- NULL
  df_both$del_score <- NULL
  
  df_both$sign_mean_log2FC_1 <- sign_mean_log2FC_1
  df_both$sign_mean_log2FC_2 <- sign_mean_log2FC_2
  df_both$sign_mean_log2FC_3 <- sign_mean_log2FC_3
  
  output_list <- list(plots = gridplot,
                      results = df_both)
  
  return(output_list) 
}

mode = "lm"
cancer_type <- "BRCA"
backbone_path <- "/home/ieo7429/Desktop/THESIS_GAB/tables/backbone_0.1Mbp.tsv"
accessibility_table_path <- "/home/ieo7429/Desktop/THESIS_GAB/bulk_ATAC/GSE231481/peaks/DA_output_backbone_0.1Mbp.tsv"
ampl_del_table_path <- "/home/ieo7429/Desktop/THESIS_GAB/tables/training_dataset_wo_target_w_CNA_0.1Mbp_hg19.RData"

outlist <- REGRESS_OUT_COHORT(cancer_type = cancer_type, 
                              backbone_path = backbone_path, 
                              accessibility_table_path = accessibility_table_path, 
                              ampl_del_table_path = ampl_del_table_path)

regrout_df <- outlist$results
regrout_df <- regrout_df %>% 
                dplyr::select(bin, logFC, bool_diff_acc, 
                              sign_mean_log2FC_1, sign_mean_log2FC_2, sign_mean_log2FC_3)


outpath <- "/home/ieo7429/Desktop/THESIS_GAB/bulk_ATAC/GSE231481/peaks/regrout_DA_output_backbone_0.1Mbp.tsv"
write.table(x = regrout_df, 
            file = outpath, 
            sep = "\t", col.names = T)