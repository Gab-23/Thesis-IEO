setwd("~/Desktop/THESIS_GAB/ML/Regression")

extract_gene_names_from_windows <- function(region_path, outname, genes, write = FALSE, what = "all"){
  
  interesting_regions <- read.table(region_path, header = T, sep = "\t")
  accepted_cols <- c("is_good_residual")
  
  for (col in accepted_cols){
    if(col %in% colnames(interesting_regions)) {
      interesting_regions[[col]] <- ifelse(interesting_regions[[col]] == "True", TRUE, FALSE)
    }
  }
  
  if (!all(startsWith(as.character(interesting_regions$chr), "chr"))){
    interesting_regions$chr <- paste0("chr",interesting_regions$chr)
  }
  
  if (what == "+") {
    interesting_regions <- interesting_regions[interesting_regions$observed > 0, ]
  } else if (what == "-") {
    interesting_regions <- interesting_regions[interesting_regions$observed < 0, ]
  }
  
  grang <- GRanges(interesting_regions)
  
  overlaps <- subsetByOverlaps(genes, grang)
  
  gene_names <- annotate::getSYMBOL(x = overlaps$gene_id, data='org.Hs.eg')
  names(gene_names) <- NULL
  gene_names <- sort(gene_names)
  
  if(write){
    
    df <- as.data.frame(gene_names)
    write.table(x = df, file = outname, 
                col.names = F, 
                row.names = F, 
                sep = "\t", 
                quote = F)
    
  }
  
  return(gene_names)
  
}

# region_paths <- c("cluster_0_new.csv", "cluster_1_new.csv", "cluster_2_new.csv", "cluster_3_new.csv", "cluster_4_new.csv")
# outnames <- c()
# for (path in region_paths) {
#   outname <- paste0("gene_names_",strsplit(path, split = "[_.]")[[1]][2],"_new") 
#  outnames <- c(outnames, outname)
# } 

region_paths_2 <- c("regions_005_0.1.csv", "regions_001_0.1.csv", "regions_0001_0.1.csv")
names(region_paths_2) <- c("gene_names_filter_005", "gene_names_filter_001")

genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

gene_lists <- list()

for (idx in seq_along(along.with = region_paths_2)) {
  
  path <- region_paths_2[idx]
  outname <- region_paths_2[idx]
  alternatives <- c("all", "+", "-")
  
  for(alt in alternatives){
    
    gene_list <- extract_gene_names_from_windows(path, paste0(outname, "_", alt), genes, FALSE, what = alt)
    gene_lists[[paste0(path, "_", alt)]] <- gene_list
  
  }
}

# this snippet takes as input a named list of gene names (a list of vectors)
# and performs GO enrichment analysis over a set of ontologies 
# (Biological Process, Cellular Component, Molecular Function)

pvalueCutoff <- 0.01 # p-value
qvalueCutoff <- 0.1  # q-value, NOTA: q-value secondo clusterProfiler non è la correzione del p-value, ma il minimo p-value adjusted per avere risultati significativi
showCategory <- 15 # top terms to show

go_res <- lapply(X = seq_along(gene_lists), FUN = function(idx){
  
  list_name <- names(gene_lists[idx]) # take sublist name
  list <- gene_lists[[idx]] # take sublist
  
  print(paste0("Processing ", list_name))
  
  gene_df <- bitr(list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # take list of names and convert to entrez ID
  entrez_ids <- gene_df$ENTREZID
  
  ontologies <- c("BP","CC","MF") # define ontologies
  
  to.return <- lapply(X = seq_along(ontologies), FUN = function(idx_ont){
    
    ont <- ontologies[idx_ont] # extract ontology
    print(paste0("Processing ", ont))
    
    go_results <- enrichGO(gene             = entrez_ids, # input genes
                           OrgDb            = org.Hs.eg.db, # database
                           keyType          = "ENTREZID",
                           ont              = ont,
                           pAdjustMethod    = "BH", # p-value correction
                           pvalueCutoff     = pvalueCutoff,
                           qvalueCutoff     = qvalueCutoff)
    
    if (nrow(go_results) > 0) { # if at least one term is enriched
      go_results <- pairwise_termsim(go_results)
      dotplot <- dotplot(go_results, showCategory = showCategory) 
      barplot <- barplot(go_results, showCategory = showCategory)
      emapplot <- emapplot(go_results)
      
      return(list(dotplot = dotplot, #return plots and output
                  barplot = barplot, 
                  emapplot = emapplot, 
                  result = go_results))
    
      } else {
      print(paste0("No enrichment results for ", list_name, " - ", ont))
      return(NULL)
    }
    
  })
  names(to.return) <- ontologies
  return(to.return)
})
names(go_res) <- names(gene_lists)

# does the same but with KEGG
kegg_res <- lapply(X = seq_along(gene_lists), FUN = function(idx){
  
  list_name <- names(gene_lists[idx])
  list <- gene_lists[[idx]]
  
  print(paste0("Processing ", list_name))
  
  gene_df <- bitr(list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  entrez_ids <- gene_df$ENTREZID
  
  kegg_results <- enrichKEGG(gene          = entrez_ids,
                             organism      = 'hsa',
                             pAdjustMethod = "BH", 
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.1)
  

  if (nrow(kegg_results) > 0) {
    
    kegg_results <- pairwise_termsim(kegg_results)
    dotplot <- dotplot(kegg_results, showCategory = 10)
    barplot <- barplot(kegg_results, showCategory = 10)

    return(list(dotplot = dotplot, 
                barplot = barplot, 
                result = kegg_results))
    
  } else {
    print(paste0("No enrichment results for KEGG - ", list_name))
    return(NULL)
  }
})
names(kegg_res) <- names(gene_lists)


cancer_gene_list <- read.table("/home/ieo7429/Desktop/THESIS_GAB/tables/ONCO_TSG_TABLE_hg19_0.1Mbp.tsv", header = T, sep = "\t", fill = TRUE)
