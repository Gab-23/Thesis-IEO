############### REGRESS OUT CNAs ###############

# The rationale behind regressing out is to correct the values of a response variable
# using the residuals. That is because residuals are what the model can't explain.
# In my case, I fit a model to predict change in chromatin accessibility using CNAs
# Residuals are those values of log2FC that cannot be explained by CNAs

# load amplification / deletion table
load("/home/ieo7429/Desktop/THESIS_GAB/tables/datasets_different_levels_total.RData")
ampl_del_table <- dataset_tables_different_levels$tot_level
ampl_del_table <- ampl_del_table[,c(1,3,4)]

write.table(ampl_del_table, file = "/home/ieo7429/Desktop/THESIS_GAB/tables/ampl_del_bin_total.bed", 
            sep = "\t", col.names = F, row.names = F, quote = F)

# load accessibility tables
acc_table <- read.table("/home/ieo7429/Desktop/THESIS_GAB/outfiles/no_processing/collapsed_table_filtered_no_processing_BRCA_NON_BASAL_1Mbp_LR.tsv", header = T, sep = "\t")

rownames(acc_table) <- paste0(acc_table$chr, "-", acc_table$start, "-", acc_table$end)

acc_table <- acc_table[rownames(acc_table) %in% rownames(acc_table_QM),]

# create GRanges
ampl_del_granges_hg19 <- GRanges(seqnames = ampl_del_table_hg19$chr, 
                                 ranges = IRanges(start = ampl_del_table_hg19$start,
                                                  end = ampl_del_table_hg19$end))

acc_granges <- GRanges(seqnames = acc_table$chr, 
                       ranges = IRanges(start = acc_table$start,
                                        end = acc_table$end))


# compute overlaps
overlaps_hg19 <- findOverlaps(acc_granges, ampl_del_granges_hg19)

# resize the tables according to overlaps

acc_granges_filtered_hg19 <- acc_granges[queryHits(overlaps_hg19)]

ampl_del_granges_filtered_hg19 <- ampl_del_granges_hg19[subjectHits(overlaps_hg19)]

intervals_hg19 <- pintersect(acc_granges_filtered_hg19, ampl_del_granges_filtered_hg19)

intervals_width_hg19 <- width(intervals_hg19)

intervals_df_hg19 <- data.frame(query = queryHits(overlaps_hg19), # prepare a df with query region
                           subject = subjectHits(overlaps_hg19), # subject region
                           widths = intervals_width_hg19) # and width of the interval

best_hits_hg19 <- intervals_df_hg19 %>% # take the best hits in general:
  group_by(query) %>% # group by the query column
  slice_max(widths, n = 1) %>% # take only the biggest overlap
  ungroup() %>%
  group_by(subject) %>% # group by the subject column
  slice_max(widths, n = 1) %>% # take only the biggest overlap
  ungroup()

worst_hits_hg19 <- intervals_df_hg19 %>% # take the worst hits in general:
  group_by(query) %>% # group by the query column
  slice_min(widths, n = 1) %>% # take only the biggest overlap
  ungroup() %>%
  group_by(subject) %>% # group by the subject column
  slice_max(widths, n = 1) %>% # take only the biggest overlap
  ungroup()

acc_table_filtered <- acc_table[best_hits$query,] # subset the data.frames
ampl_del_table_lifted_filtered <- ampl_del_table_lifted[best_hits$subject,]


# create data.frames

df <- data.frame(chr = acc_table_filtered$chr,
                 start = acc_table_filtered$start,
                 end = acc_table_filtered$end,
                 acc = acc_table_filtered$mean_log2FC,
                 ampl = ampl_del_table_lifted_filtered$ampl_score, 
                 del = ampl_del_table_lifted_filtered$del_score)


# model accessibility variation in terms of amplifications and/or deletions
# NOTA: ifg we inspect the models we see that models from Quantile mapping data
# is worse than models on standard data. 
# That is because Quantile mapping already decorrelated variables a bit

model_ampl <- summary(lm(formula = acc ~ ampl, data = df))
model_del <- summary(lm(formula = acc ~ del, data = df))
model_both <- summary(lm(formula = acc ~ ampl + del, data = df))

# set up all needed data.frames for benchmarking
df_regrout_ampl <- df
df_regrout_ampl$acc <- model_ampl$residuals
df_regrout_del <- df
df_regrout_del$acc <- model_del$residuals

df_both <- df
df_both$acc <- model_both$residuals

############### PLOT CORRELATION BENCHMARKING + PLOT VARIATION AFTER REGRESS OUT ###############

ampl_acc_test <- cor.test(df$ampl, df$acc, method = "spearman") # perform correlation test. spearman because lack of normality for both variables
cor_ampl_acc <- round(ampl_acc_test$estimate, 2) # take correlation value
pval_ampl_acc <- round(ampl_acc_test$p.value, 5) # take pvalue

ampl_acc <- ggplot(df, aes(x = ampl, y = acc, colour = chr)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  annotate("text", x = max(df$ampl) * 0.7, y = max(df$ampl) * 0.9, 
           label = paste("Correlation: ", cor_ampl_acc, "\n",
                         "P-value: ", pval_ampl_acc), 
           color = "black", size = 5, fontface = "bold") +
  labs(title = "[BEFORE NORMALIZATION]",
       x = "Amplification score",
       y = "Chromatin accessibility [mean Log2FC]") +
  theme_minimal()

ampl_acc_regrout_test <- cor.test(df_regrout_ampl$ampl, df_regrout_ampl$acc, method = "spearman")
cor_ampl_acc_regrout <- round(ampl_acc_regrout_test$estimate, 2)
pval_ampl_acc_regrout <- round(ampl_acc_regrout_test$p.value, 5)

ampl_acc_regrout <- ggplot(df_regrout_ampl, aes(x = ampl, y = acc, colour = chr)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  annotate("text", x = max(df_regrout_ampl$ampl) * 0.7, y = max(df_regrout_ampl$ampl) * 0.9, 
           label = paste("Correlation: ", cor_ampl_acc_regrout, "\n",
                         "P-value: ", pval_ampl_acc_regrout), 
           color = "black", size = 5, fontface = "bold") +
  labs(title = "[REGRESS OUT AMPLIFICATION EFFECT]",
       x = "Amplification score",
       y = "Chromatin accessibility [mean Log2FC]") +
  theme_minimal()

ampl_acc_regrout_both_test <- cor.test(df_both$ampl, df_both$acc, method = "spearman")
cor_ampl_acc_regrout_both <- round(ampl_acc_regrout_both_test$estimate, 2)
pval_ampl_acc_regrout_both <- round(ampl_acc_regrout_both_test$p.value, 5) 

ampl_acc_regrout_both <- ggplot(df_both, aes(x = ampl, y = acc, colour = chr)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  annotate("text", x = max(df_both$ampl) * 0.7, y = max(df_both$ampl) * 0.9, 
           label = paste("Correlation: ", cor_ampl_acc_regrout_both, "\n",
                         "P-value: ", pval_ampl_acc_regrout_both), 
           color = "black", size = 5, fontface = "bold") +
  labs(title = "[REGRESS OUT AMPL & DEL EFFECT]",
       x = "Amplification score",
       y = "Chromatin accessibility [mean Log2FC]") +
  theme_minimal()

png(filename = paste0("/home/ieo7429/Desktop/THESIS_GAB/plots/correlation_plot_ampl.png"), width = 12000, 
    height = 9000, 
    res = 800)
gridExtra::grid.arrange(ampl_acc, ampl_acc_regrout, ampl_acc_regrout_both,
                        nrow = 3, ncol = 1)
dev.off()

del_acc_test <- cor.test(df$del, df$acc, method = "spearman")
cor_del_acc <- round(del_acc_test$estimate, 2)
pval_del_acc <- round(del_acc_test$p.value, 5) 

del_acc <- ggplot(df, aes(x = del, y = acc, colour = chr)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  annotate("text", x = max(df$del) * 0.7, y = max(df$del) * 0.9, 
           label = paste("Correlation: ", cor_del_acc, "\n",
                         "P-value: ", pval_del_acc), 
           color = "black", size = 5, fontface = "bold") +
  labs(title = "[BEFORE NORMALIZATION]",
       x = "Deletion score",
       y = "Chromatin accessibility [mean Log2FC]") +
  theme_minimal()

del_acc_regrout_test <- cor.test(df_regrout_del$del, df_regrout_del$acc, method = "spearman")
cor_del_acc_regrout <- round(del_acc_regrout_test$estimate, 2)
pval_del_acc_regrout <- round(del_acc_regrout_test$p.value, 5)

del_acc_regrout <- ggplot(df_regrout_del, aes(x = del, y = acc, colour = chr)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  annotate("text", x = max(df_regrout_del$del) * 0.7, y = max(df_regrout_del$del) * 0.9, 
           label = paste("Correlation: ", cor_del_acc_regrout, "\n",
                         "P-value: ", pval_del_acc_regrout), 
           color = "black", size = 5, fontface = "bold") +
  labs(title = "[REGRESS OUT DELETION EFFECT]",
       x = "Deletion score",
       y = "Chromatin accessibility [mean Log2FC]") +
  theme_minimal()

del_acc_regrout_both_test <- cor.test(df_both$del, df_both$acc, method = "spearman")
cor_del_acc_regrout_both <- round(del_acc_regrout_both_test$estimate, 2)
pval_del_acc_regrout_both <- round(del_acc_regrout_both_test$p.value, 5) 

del_acc_regrout_both <- ggplot(df_both, aes(x = del, y = acc, colour = chr)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  annotate("text", x = max(df_both$del) * 0.7, y = max(df_both$del) * 0.9, 
           label = paste("Correlation: ", cor_del_acc_regrout_both, "\n",
                         "P-value: ", pval_del_acc_regrout_both), 
           color = "black", size = 5, fontface = "bold") +
  labs(title = "[REGRESS OUT AMPL & DEL EFFECT]",
       x = "Deletion score",
       y = "Chromatin accessibility [mean Log2FC]") +
  theme_minimal()

png(filename = paste0("/home/ieo7429/Desktop/THESIS_GAB/plots/correlation_plot_del_prova.png"), width = 12000, 
    height = 9000, 
    res = 800)
gridExtra::grid.arrange(del_acc, 
                        del_acc_regrout,
                        del_acc_regrout_both,
                        nrow = 3, ncol = 1)
dev.off()


# regression out approach considering both AMPLIFICATION and DELETION is clearly the best one, 
# so we pick it and make it look like a collapsed DACRs df.

df_both$mean_log2FC <- df_both$acc
df_both$perc_of_significant_per_bin <- acc_table_filtered$perc_of_significant_per_bin
df_both$weighted_log2FC <- abs(df_both$mean_log2FC) * df_both$perc_of_significant_per_bin
df_both$acc <- NULL
df_both$ampl <- NULL
df_both$del <- NULL

df_both <- process_DACRs_coords(df_both, chrom_cum_start)

whole_genome_plot_name_regrout <- "manhattan_whole_genome_regrout_both_provas"
manhattan_plot_whole <- plot_log2FC_collapsed(df_both, TSG_ONCOGENES_DF, n_of_samples, ALPHA, CNC, test.use) # plot it

png(filename = paste0(PLOTS_OUTDIR, whole_genome_plot_name_regrout, SUFFIX, "_", test.use, ".png"), width = 21000, height = 7000, res = 800)
manhattan_plot_whole
dev.off()






