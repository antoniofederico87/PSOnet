rm(list=ls())

install.packages("openxlsx", dependencies = TRUE)
library(openxlsx)
library(tidyverse)
library(clusterProfiler)
library(ReactomePA)
library(TopKLists)

mainDir <- "/path/to/differential/expression/files"
#subDir_PSO <- "/GEO/LS_vs_NL/LS_vs_NL_eUtopia/PSO"
#subDir_AD <- "/GEO/LS_vs_NL/LS_vs_NL_eUtopia/AD"


setwd(paste0(mainDir))


DEfiles <- list.files(path = ".", recursive = TRUE, pattern = "Differential_Expression_Table_GSE*")
DEfiles <- DEfiles[-21]  ### Remove GSE6710 dataset, only 1 DEG in PSO
#DEfiles <- DEfiles[-c(5, 6, 7, 8, 9, 11)]  ### Remove GSE27887, GSE32924, GSE36842, GSE5667_GPL96, GSE5667_GPL97, GSE59294 no DEGs found in AD


datasets_list <- list()

for (file in 1:length(DEfiles)) {
  print(paste("Loading ", DEfiles[file], "in the global environment."))
  xlsx_file <- read.xlsx(DEfiles[file], sheet = 2, colNames = TRUE, rowNames = FALSE)
  xlsx_file <- as_tibble(xlsx_file)
  xlsx_file <- xlsx_file %>% mutate(ID=sapply(xlsx_file$ID, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][1]))
  datasets_list[[file]] <- xlsx_file
}

names(datasets_list) <- sapply(DEfiles, function(x) strsplit(x, split = "/", fixed = TRUE)[[1]][1])


### Fix GSE57255 Identifiers problem 
library(biomaRt)
# define biomart object
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# query biomart
results <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id"),
                 filters = "ensembl_transcript_id", values = datasets_list[["GSE57225"]]$ID,
                 mart = mart)
results

length(datasets_list[["GSE57225"]]$ID)

datasets_list[["GSE57225"]]$goodID <- NA
datasets_list[["GSE57225"]]$goodID <- results$ensembl_gene_id[match(datasets_list[["GSE57225"]]$ID, results$ensembl_transcript_id)]
datasets_list[["GSE57225"]] <- na.omit(datasets_list[["GSE57225"]])
datasets_list[["GSE57225"]]$ProbeID <- NULL
datasets_list[["GSE57225"]]$ID <- NULL
datasets_list[["GSE57225"]] <- datasets_list[["GSE57225"]] %>% group_by(goodID) %>% summarise_all(.funs = median)
datasets_list[["GSE57225"]]$ID <- datasets_list[["GSE57225"]]$goodID
datasets_list[["GSE57225"]]$goodID <- NULL

###############

pathways_list <- list()

# The dataset GSE57225 doesn not enrich any pathway [index 20]

for (data in 21:length(datasets_list)) {
  eg = clusterProfiler::bitr(datasets_list[[data]]$ID, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  path <- enrichPathway(gene=eg$ENTREZID, organism = "human", pAdjustMethod = "BH", pvalueCutoff=0.05, qvalueCutoff = 0.1, readable=T)
  
  p1 <- dotplot(path) + ggtitle(paste0("Dotplot ", names(datasets_list[data])))
  ggsave(plot = p1, filename = paste0("Dotplot_", names(datasets_list[data]), ".png"), device = "png", width = 8.6, height = 6, units = "in")
  p2 <- cnetplot(path, circular=TRUE, colorEdge=TRUE) + ggtitle(paste0("Cnetplot ", names(datasets_list[data])))
  ggsave(plot = p2, filename = paste0("Cnetplot_", names(datasets_list[data]), ".png"), device = "png", width = 8.6, height = 6, units = "in")
  p3 <- emapplot(path) + ggtitle(paste0("Emapplot ", names(datasets_list[data])))
  ggsave(plot = p3, filename = paste0("Emapplot_", names(datasets_list[data]), ".png"), device = "png", width = 8.6, height = 6, units = "in")
  
  p4 <- heatplot(path) + ggtitle(paste0("Heatplot ", names(datasets_list[data])))
  #ggsave(plot = p4, filename = paste0("Heatplot ", names(datasets_list[data]), ".png"), device = "png")
  
  #cowplot::plot_grid(p1, p2, p3, p4, ncol=2)
  
  sigpath <- path@result[which(path@result$p.adjust<0.05),]
  dataset <- strsplit(DEfiles[data], split = "/", fixed = TRUE)[[1]][1]
  print(dataset)
  sigpath$dataset <- rep(dataset, length(sigpath$ID))
  pathways_list[[data]] <- sigpath
}

names(pathways_list) <- sapply(DEfiles, function(x) strsplit(x, split = "/", fixed = TRUE)[[1]][1])

#dotplot(path, showCategory=15)
#bp=barplot(sigpath, showCategory=10, title=paste0("barplot_mod_", mod))
#enrichMap(sigpath, layout=igraph::layout_nicely, vertex.label.cex = 1, title=paste0("enrichmap_mod_", mod))

pathways_overall <- do.call("rbind", pathways_list) ### Dataframe containing all the pathways enriched in all the datasets

length(unique(pathways_overall$ID))
path_occurrence <- sort(table(pathways_overall$Description), decreasing = TRUE)

# Build occurrence matrix pathways/datasets

pathways_occ_mat <- matrix(0, nrow = length(unique(pathways_overall$ID)), ncol = length(DEfiles), dimnames = list(unique(pathways_overall$ID), names(pathways_list)))

for (pos in 1:length(names(pathways_list))) {
  datapath <- pathways_overall[which(pathways_overall$dataset %in% names(pathways_list)[pos]), "ID"]
  if(length(datapath)==0) {
    next
  }
  for (path in 1:length(datapath)){
    pathways_occ_mat[datapath[path], names(pathways_list)[pos]] <- 1
  }
}

saveRDS(pathways_occ_mat, file = "pathways_occurrence.rds")

###### DEG Borda ranking #####

DEG_ranks_list <- list()
logFC_ranks_list <- list()
DEG_vec <- c()
logFC_vec <- c()
darioscore_list <- list()

for (DEGlist in 1:length(datasets_list)) {
  DEG_vec <- c(abs(datasets_list[[DEGlist]]$logFC*-log(datasets_list[[DEGlist]]$adj.P.Val)))
  logFC_vec <- datasets_list[[DEGlist]]$logFC
  names(DEG_vec) <- names(logFC_vec) <-  datasets_list[[DEGlist]]$ID
  logFC_ranks_list[[DEGlist]] <- data.frame(gene=names(DEG_vec), logFC=datasets_list[[DEGlist]]$logFC)
  DEG_vec <- sort(DEG_vec, decreasing = TRUE)
  logFC_vec <- sort(logFC_vec, decreasing = FALSE)
  darioscore_list[[DEGlist]] <- DEG_vec
  logFC_ranks_list[[DEGlist]] <- logFC_vec
  finalDEGrank_conv <- clusterProfiler::bitr(names(DEG_vec), fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop = FALSE)
  DEG_ranks_list[[DEGlist]] <- na.omit(finalDEGrank_conv)
}

names(DEG_ranks_list) <- names(datasets_list)


list_to_borda <- list()

for (i in 1:length(DEG_ranks_list)) {
  list_to_borda[[i]] <- DEG_ranks_list[[i]]$SYMBOL
  
}

names(list_to_borda) <- names(datasets_list)

finalDEGrank <- TopKLists::Borda(list_to_borda)
finalDEGrank <- finalDEGrank$TopK$median # Rank of "most to least DEGs across studies"

logFC_ranks_names <- list()
for (i in 1:length(logFC_ranks_list)){
  logFC_ranks_names[[i]] <- names(logFC_ranks_list[[i]])
}

borda_on_logFC <- TopKLists::Borda(logFC_ranks_names)
borda_on_logFC <- borda_on_logFC$TopK$median

write.table(finalDEGrank_conv, file = "finalrank_PSO.txt", sep = "\t", quote = FALSE)

genes_overall <- unlist(list_to_borda) ### vector containing all the genes enriched in all the datasets

genes_occurrence <- sort(table(genes_overall), decreasing = TRUE) # Rank of "most occurrent genes across studies"
genes_occurrence[1:10]

### Rank entire Limma results ###

entire_files <- list.files(path = ".", recursive = TRUE, pattern = "_unfiltered.xlsx")
#entire_files <- entire_files[-21] # PSO
#entire_files <- entire_files[-c(5, 6, 7, 8, 9, 11)] # AD

unfiltered_limma <- list()

for (i in 1:length(entire_files)){
  xlsx_file <- read.xlsx(entire_files[i], sheet = 2, colNames = TRUE, rowNames = FALSE)
  xlsx_file <- as_tibble(xlsx_file)
  xlsx_file <- xlsx_file %>% mutate(ID=sapply(xlsx_file$ID, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][1]))
  unfiltered_limma[[i]] <- xlsx_file
}

names(unfiltered_limma) <- sapply(entire_files, function(x) strsplit(x, split = "/", fixed = TRUE)[[1]][1])


### Fix GSE57255 Identifiers problem 
# library(biomaRt)
# # define biomart object
# mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# # query biomart
# results <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id"),
#                  filters = "ensembl_transcript_id", values = unfiltered_limma[["GSE57225"]]$ID,
#                  mart = mart)
# results
# 
# length(unfiltered_limma[["GSE57225"]]$ID)
# 
# unfiltered_limma[["GSE57225"]]$goodID <- NA
# unfiltered_limma[["GSE57225"]]$goodID <- results$ensembl_gene_id[match(unfiltered_limma[["GSE57225"]]$ID, results$ensembl_transcript_id)]
# unfiltered_limma[["GSE57225"]] <- na.omit(unfiltered_limma[["GSE57225"]])
# unfiltered_limma[["GSE57225"]]$ProbeID <- NULL
# unfiltered_limma[["GSE57225"]]$ID <- NULL
# unfiltered_limma[["GSE57225"]] <- unfiltered_limma[["GSE57225"]] %>% group_by(goodID) %>% summarise_all(.funs = median)
# unfiltered_limma[["GSE57225"]]$ID <- unfiltered_limma[["GSE57225"]]$goodID
# unfiltered_limma[["GSE57225"]]$goodID <- NULL
# 
###############

for (lst in 1:length(unfiltered_limma)){
  darioscore <- c(abs(unfiltered_limma[[lst]]$logFC)*-log(unfiltered_limma[[lst]]$adj.P.Val))
  unfiltered_limma[[lst]]$darioscore <- darioscore
  unfiltered_limma[[lst]] <- unfiltered_limma[[lst]] %>% arrange(desc(darioscore))
  
}

unfiltered_list_to_borda <- list()
for (lst in 1:length(unfiltered_limma)){
  id_conv <- clusterProfiler::bitr(unfiltered_limma[[lst]]$ID, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop = FALSE)
  unfiltered_list_to_borda[[lst]] <- na.omit(id_conv$SYMBOL)
  
}

unfiltered_borda <- TopKLists::Borda(unfiltered_list_to_borda)
unfiltered_borda <- unfiltered_borda$TopK$median

which(unfiltered_borda %in% names(genes_occurrence[which(genes_occurrence==10)]))


### Visualization ###

### Barplots ###

df_genes <- data.frame(genes_occurrence)

p <- ggplot(df_genes[1:40,], aes(x=reorder(genes_overall, Freq), y = Freq)) + geom_bar(stat="identity", fill="steelblue") + 
  coord_flip()
p ### genes DE over x datasets


path_vec <- c()

for(i in 1:length(unique(df_path$Freq))) {
  path_vec <- c(path_vec, length(df_path$path[which(df_path$Freq==unique(df_path$Freq)[i])]))
  
}

df_path_vec <- data.frame(path_vec)
df_path_vec$Freq <- unique(df_path$Freq)

p <- ggplot(df_path_vec, aes(x=reorder(Freq, -Freq), y = path_vec)) + geom_bar(stat="identity", fill="steelblue")
p ### This graph expresses the number of pathways enriched in x datasets



df_path <- data.frame(path_occurrence)
colnames(df_path) <- c("path", "Freq")

p <- ggplot(df_path[1:20,], aes(x=reorder(path, Freq), y = Freq)) + geom_bar(stat="identity", fill="steelblue") + 
  coord_flip()
p

num_of_degs <- sort(sapply(list_to_borda, function(x) length(x)), decreasing = TRUE)
num_of_degs <- as.data.frame(num_of_degs)
num_of_degs$dataset <- rownames(num_of_degs)
p <- ggplot(num_of_degs, aes(x=reorder(dataset, -num_of_degs), y = num_of_degs)) + geom_bar(stat="identity", fill="steelblue") + theme(axis.text.x = element_text(angle=45))
p




  