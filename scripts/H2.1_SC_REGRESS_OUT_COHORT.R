REGRESS_OUT_COHORT <- function(cancer_type, backbone_path, accessibility_table_path, ampl_del_table_path, mode = "lm"){
  
  backbone <- read.table(backbone_path, header = T, sep = "\t")
  accessibility_table <- read.table(accessibility_table_path, header = T, sep = "\t")
  accessibility_table <- merge(x = backbone, y = accessibility_table, by = c("chr", "start", "end"))
  accessibility_table$start <- NULL; accessibility_table$end <- NULL
  
  load(ampl_del_table_path)
  ML_dataset_with_target <- ML_dataset_with_target[ML_dataset_with_target$Type == cancer_type,]
  ampl_del_table <- ML_dataset_with_target[, c(1,3,4)]
  
  df <- merge(x = ampl_del_table, y = accessibility_table, by = "bin")
  
  if (mode == "lm") {
    formula <- mean_log2FC ~ ampl_score + del_score
    model_both <- summary(lm(formula = formula, data = df))  
  } else if (mode == "poly") {
    formula <- mean_log2FC ~ poly(ampl_score, degree = 3, raw = TRUE) + poly(del_score, degree = 3, raw = TRUE)
    model_both <- summary(lm(formula = formula, data = df))
  }
  
  df_both <- df
  df_both$mean_log2FC <- model_both$residuals
  
  gridplot <- plot_CNA_acc_correlation(df, df_both, mode)
  
  log2k_1 <- log2(1) 
  log2k_2 <- log2(2)
  log2k_3 <- log2(3)
  
  sign_mean_log2FC_1 <- ifelse(df_both$mean_log2FC < -log2k_1, -1, ifelse(df_both$mean_log2FC > log2k_1, 1, 0))
  sign_mean_log2FC_2 <- ifelse(df_both$mean_log2FC < -log2k_2, -1, ifelse(df_both$mean_log2FC > log2k_2, 1, 0))
  sign_mean_log2FC_3 <- ifelse(df_both$mean_log2FC < -log2k_3, -1, ifelse(df_both$mean_log2FC > log2k_3, 1, 0))

  df_both$weighted_log2FC <- df_both$mean_log2FC * df_both$perc_of_significant_per_bin
  df_both$weighted_coef_of_var_log2FC <- (df_both$mean_log2FC / df_both$sd_log2FC) * df_both$perc_of_significant_per_bin
  
  df_both$ampl_score <- NULL
  df_both$del_score <- NULL

  output_list <- list(plots = gridplot,
                  results = df_both)
  
  return(output_list) 
}
