library(biomaRt)
library(clusterProfiler)
library(ReactomePA)
source("/path/to/INfORM_functions.R")

setwd("/path/to/gene/expression/matrices")


pso_pamr <- read.table("GE_Mic_PSO_Pamr.txt", header = TRUE, sep = "\t", quote = "")

### Biomart rownames conversion ###
listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
head(datasets)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
searchAttributes(mart = ensembl, pattern = "hgnc")

conversion_table <- biomaRt::getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = rownames(pso_pamr), mart = ensembl)

conversion_table$hgnc_symbol[which(conversion_table$hgnc_symbol=="")] <- NA
sum(is.na(conversion_table$hgnc_symbol))

#pso_pamr <- pso_pamr[-which(rownames(pso_pamr) %in% setdiff(rownames(pso_pamr), conversion_table$ensembl_gene_id)),]

pso_pamr$symbols <- conversion_table$hgnc_symbol
pso_pamr <- na.omit(pso_pamr)
rownames(pso_pamr) <- pso_pamr$symbols
pso_pamr$symbols <- NULL

###
  
generatematrices=get_ranked_consensus_matrix(gx_table = pso_pamr, iMethods = c("clr","aracne","mrnet"),
                                               iEst = c("pearson","spearman","kendall","mi.empirical","mi.mm","mi.shrink","mi.sg"),
                                               iDisc=c("none","equalfreq","equalwidth","globalequalwidth"), ncores = 20, debug_output = TRUE, updateProgress = TRUE)
  
  
  #Parse ranked matrix and get bin_mat and edge_rank
  # Get edge rank list and binary inference matrix from edge rank matrix computed by get_ranked_consensus_matrix().
  # parse_edge_rank_matrix parses the edge rank matrix created by using the internal function get_ranked_consensus_matrix_matrix() to get a ranked edge list and a binary matrix.
  
rankMat.parsed=parse_edge_rank_matrix(edge_rank_matrix = generatematrices, edge_selection_strategy = "default",
                                        mat_weights = "rank", topN = 10, debug_output = TRUE, updateProgress = TRUE)
  
conGraph <- get_iGraph(rankMat.parsed$bin_mat)
  
saveRDS(conGraph, file="PSO_pamr_conGraph.rds")



