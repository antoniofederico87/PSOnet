rm(list=ls())

library(biomaRt)
library(clusterProfiler)
library(ReactomePA)
library(igraph)
library(bc3net)
library(RColorBrewer)
library(doParallel)
library(foreach)
library(ggplot2)

setwd("/path/to/master/folder/")

source("/path/to/INfORM_functions.R")
source("BIOMAP_analysis_functions.R")



conGraph <- readRDS("PSO_pamr_conGraph.rds")
conGraph_nl <- readRDS("PSO_pamr_conGraph_nonlesional.rds")
conGraph.modules <- readRDS("conGraph_modules_PSO.rds")
conGraph_nl.modules <- readRDS("conGraph_modules_PSO_nl.rds")
named_ranklist <- readRDS("inform_rank_pso.rds")
named_ranklist_nl <- readRDS("inform_rank_pso_nl.rds")
unfiltered_limma_list <- readRDS("unfiltered_list_to_borda.rds")

#Check the length of the detected modules
allmodules <- c()

for (i in 1:length(conGraph.modules)) {
  print(length(conGraph.modules[[i]]))
  allmodules <- c(allmodules, length(conGraph.modules[[i]]))
}

allmodules_nl <- c()

for (i in 1:length(conGraph_nl.modules)) {
  print(length(conGraph_nl.modules[[i]]))
  allmodules_nl <- c(allmodules_nl, length(conGraph_nl.modules[[i]]))
}


modvec <- c()
for (k in 1:length(conGraph.modules)) {
  if (length(conGraph.modules[[k]])>=10){
    print(length(conGraph.modules[[k]]))
    modvec <- c(modvec, length(conGraph.modules[[k]]))
  }
}

modvec_nl <- c()
for (k in 1:length(conGraph_nl.modules)) {
  if (length(conGraph_nl.modules[[k]])>=10){
    print(length(conGraph_nl.modules[[k]]))
    modvec_nl <- c(modvec_nl, length(conGraph_nl.modules[[k]]))
  }
}

members=igraph::membership(conGraph.modules)
members_list <- list()
for(i in 1:length(conGraph.modules)){
  members_list[[i]] <- names(which(members==i))
}
names(members_list) <- paste0("mod_", 1:length(conGraph.modules))

members[members==1]
names(members[members==1])

members_nl=igraph::membership(conGraph_nl.modules)
members_nl[members_nl==1]
names(members_nl[members_nl==1])

members_list_nl <- list()
for(i in 1:length(conGraph_nl.modules)){
  members_list_nl[[i]] <- names(which(members_nl==i))
}
names(members_list_nl) <- paste0("mod_", 1:length(conGraph_nl.modules))



mainDir <- "/path/to/master/folder/"
subDir <- "Functional_annotation_PSO"
getwd()


get_reactome_from_modules(igraph_modules = conGraph.modules, geneID = "SYMBOL", pval_cutoff = 0.05, outPath = "/home/MOKA/veera/Vittorio/Functional_annotation_PSO", layout = "overall")
get_reactome_from_modules(igraph_modules = conGraph_nl.modules, geneID = "SYMBOL", pval_cutoff = 0.05, outPath = "/home/MOKA/veera/Vittorio/Functional_annotation_PSO_nl/", layout = "overall")

get_bubbleplot_from_pathways(igraph_modules = conGraph.modules, geneID = "SYMBOL")
get_bubbleplot_from_pathways(igraph_modules = conGraph_nl.modules, geneID = "SYMBOL")

### Immunologic signatures

m_df_immuno = msigdbr(species = "Homo sapiens", category = "C7")

signature_immuno <- c()
hallmarks_immuno <- list()

for (sig in 1:length(names(table(m_df_immuno$gs_name)))) {
  
  signature_immuno <- m_df_immuno[which(m_df_immuno$gs_name==names(table(m_df_immuno$gs_name))[sig]),"human_gene_symbol"] %>% pull()
  hallmarks_immuno[[sig]] <- signature_immuno
}

names(hallmarks_immuno) <- names(table(m_df_immuno$gs_name))

occ_mat_immuno <- matrix(data = NA, nrow = length(conGraph.modules), ncol = length(hallmarks_immuno), dimnames = list(paste0("mod_", seq(1:length(conGraph.modules))), names(hallmarks_immuno)))

for(mod in 1:length(conGraph.modules)){
  for(hal in 1:length(hallmarks_immuno)) {
    
    ji <- length(intersect(conGraph.modules[[mod]], hallmarks_immuno[[hal]]))/length(union(conGraph.modules[[mod]], hallmarks_immuno[[hal]]))
    occ_mat_immuno[mod, hal] <- ji
    
  }
}

colnames(occ_mat_immuno) <- sapply(colnames(occ_mat_immuno), function(x) substr(x, start = 10, stop = nchar(x)))
heatmap.2(occ_mat_immuno, Rowv = FALSE, dendrogram = "column", symbreaks = FALSE, trace = "none", cexCol = 0.8, margins = c(15,5))

### Enrichment analysis ###
### Import gene lists ###
genes_occurrence <- readRDS("genes_occurrence_PSO.rds")
biomarkers_wp6 <- readxl::read_xlsx("Biomarkers Integrity_v2 Psoriasis.xlsx", sheet = 2, col_names = TRUE)
eg = clusterProfiler::bitr(unique(biomarkers_wp6$ENTREZ), fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
gwas_PSO_dand <- read.table("Lena_gene_lists/GWAS_PSO.Dand_et_al.txt", header = FALSE, sep="\t", quote = "")
meth_PSO_lena <- read.table("Lena_gene_lists/Methylation_PSO.Möbus_et_al.txt", header = FALSE, sep="\t", quote = "")
meth_multiple_lena <- read.table("Lena_gene_lists/Methylation_in_2_or_more_studies_PSO.Möbus_et_al.txt", header = FALSE, sep="\t", quote = "")
PSO_list <- list(names(genes_occurrence), eg$SYMBOL, gwas_PSO_dand$V1, meth_PSO_lena$V1, meth_multiple_lena$V1)
names(PSO_list) <- c("DEGs_occurrence", "wp6_biomarker", "GWAS_PSO_Dand", "meth_PSO_lena", "meth_multiple_lena")
relevant_genes <- Reduce(union, PSO_list)

relevant_on_modules <- map_genes_on_modules(geneSet = relevant_genes, igraph_modules = conGraph.modules)

relevant_score_mat <- matrix(0, nrow = length(relevant_genes), ncol = length(PSO_list), dimnames = list(relevant_genes, names(PSO_list)))

for(i in 1:length(PSO_list)){
  pos <- which(relevant_genes %in% PSO_list[[i]])
  relevant_score_mat[pos,i] <- 1
}
relevant_score_mat <- as.data.frame(relevant_score_mat)
relevant_score_mat[names(genes_occurrence), "DEGs_occurrence"] <- genes_occurrence
relevant_score_mat$sum_score <- rowSums(relevant_score_mat)

relevant_score_mat_mapped <- relevant_score_mat[which(rownames(relevant_score_mat) %in% as_ids(V(conGraph))),]
relevant_score_mat_mapped$info_perc <- relevant_score_mat_mapped$sum_score/sum(relevant_score_mat_mapped$sum_score)

relevant_score_mat_unmapped <- relevant_score_mat[-which(rownames(relevant_score_mat) %in% as_ids(V(conGraph))),]


### Lists of genes falling in each module ###
mod_score_list <- list()
for (i in 1:length(conGraph.modules)){
  mod_score_list[[i]] <- relevant_score_mat_mapped[which(rownames(relevant_score_mat_mapped) %in% conGraph.modules[[i]]),]
}

mod_score_list_nl <- list()
for (i in 1:length(conGraph_nl.modules)){
  mod_score_list_nl[[i]] <- relevant_score_mat_mapped[which(rownames(relevant_score_mat_mapped) %in% conGraph_nl.modules[[i]]),]
}

### Modules score vector ###
modules_score <- c()
for(i in 1:length(conGraph.modules))
  modules_score <- c(modules_score, sum(mod_score_list[[i]]$info_perc))

names(modules_score) <- paste0("mod_", 1:length(conGraph.modules))

modules_score_nl <- c()
for(i in 1:length(conGraph_nl.modules))
  modules_score_nl <- c(modules_score_nl, sum(mod_score_list_nl[[i]]$info_perc))

names(modules_score_nl) <- paste0("mod_", 1:length(conGraph_nl.modules))

bar_df <- data.frame(modules_score)
bar_df_nl <- data.frame(modules_score_nl)

ggplot(data=bar_df, aes(x=rownames(bar_df), y=modules_score)) +
  geom_bar(stat="identity", fill="steelblue")+
 # geom_text(aes(label=len), vjust=-0.3, size=3.5)+
  theme_minimal()

ggplot(data=bar_df_nl, aes(x=rownames(bar_df_nl), y=modules_score_nl)) +
  geom_bar(stat="identity", fill="steelblue")+
  # geom_text(aes(label=len), vjust=-0.3, size=3.5)+
  theme_minimal()


nodes_combn <- t(combn(x = rownames(relevant_score_mat_mapped), m = 2))
nodes_combn <- as.data.frame(nodes_combn)
nodes_combn$info_score_v1 <- NA
nodes_combn$info_score_v2 <- NA
nodes_combn$info_score_v1 <- relevant_score_mat_mapped$info_perc[match(nodes_combn$V1, rownames(relevant_score_mat_mapped))]
nodes_combn$info_score_v2 <- relevant_score_mat_mapped$info_perc[match(nodes_combn$V2, rownames(relevant_score_mat_mapped))]


shortpath <- shortest.paths(conGraph, v = rownames(relevant_score_mat_mapped), to = rownames(relevant_score_mat_mapped), weights = NULL)
shortpath_nl <- shortest.paths(conGraph_nl, v = rownames(relevant_score_mat_mapped), to = rownames(relevant_score_mat_mapped), weights = NULL)

outvec <- list()
#ptm <- proc.time()
outvec <- mclapply(1:length(nodes_combn$V1), function(i){
  shortpath[nodes_combn$V1[i], nodes_combn$V2[i]]
}, mc.preschedule = TRUE, mc.cores = 40)

outvec_nl <- list()
#ptm <- proc.time()
outvec_nl <- mclapply(1:length(nodes_combn$V1), function(i){
  shortpath_nl[nodes_combn$V1[i], nodes_combn$V2[i]]
}, mc.preschedule = TRUE, mc.cores = 40)


#proc.time() - ptm
nodes_combn$shortpath <- unlist(outvec)
nodes_combn$shortpath_nl <- unlist(outvec_nl)

#nodes_combn <- readRDS(file = "nodes_combn.rds")

nodes_combn$final_shortpathscore <- rowMeans(cbind(nodes_combn$info_score_v1, nodes_combn$info_score_v2))/nodes_combn$shortpath
nodes_combn$final_shortpathscore_nl <- rowMeans(cbind(nodes_combn$info_score_v1, nodes_combn$info_score_v2))/nodes_combn$shortpath_nl


### Calculate intermediate genes between genes with shortpath=1

nodes_combn_shortpath_1 <- nodes_combn[which(nodes_combn$shortpath==2),]
shortpath_1_list <- list()

# ptm <- proc.time()
# shortpath_1_list <- mclapply(1:nrow(nodes_combn_shortpath_1), function(i) {
#   get.all.shortest.paths(conGraph, from = nodes_combn_shortpath_1[i,1], to = nodes_combn_shortpath_1[i,2])
# }, mc.preschedule = TRUE, mc.cores = 40)
# proc.time() - ptm

shortpath_1_list <- get.all.shortest.paths(conGraph, from = union(nodes_combn_shortpath_1$V1, nodes_combn_shortpath_1$V2), to = union(nodes_combn_shortpath_1$V1, nodes_combn_shortpath_1$V2))
tmp <- unlist(shortpath_1_list$res)
tmp2 <- table(names(tmp))
intermediate_genes <- sort(tmp2[-which(names(tmp2) %in% union(nodes_combn_shortpath_1$V1, nodes_combn_shortpath_1$V2))], decreasing = TRUE)


nodes_combn_shortpath_1_nl <- nodes_combn[which(nodes_combn$shortpath_nl==2),]
shortpath_1_list_nl <- list()
shortpath_1_list_nl <- igraph::get.all.shortest.paths(conGraph_nl, from = union(nodes_combn_shortpath_1_nl$V1, nodes_combn_shortpath_1_nl$V2), to = union(nodes_combn_shortpath_1_nl$V1, nodes_combn_shortpath_1_nl$V2))
tmp <- unlist(shortpath_1_list_nl$res)
tmp2 <- table(names(tmp))
intermediate_genes_nl <- sort(tmp2[-which(names(tmp2) %in% union(nodes_combn_shortpath_1_nl$V1, nodes_combn_shortpath_1_nl$V2))], decreasing = TRUE)

final_bridge_genes <- setdiff(names(intermediate_genes), names(intermediate_genes_nl))
final_bridge_genes_lesional <- intermediate_genes[which(names(intermediate_genes) %in% final_bridge_genes)]

convgenes = clusterProfiler::bitr(final_bridge_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
bridge_pathways <- ReactomePA::enrichPathway(gene=convgenes$ENTREZID,pvalueCutoff=1, readable=T)
ReactomePA::dotplot(bridge_pathways, showCategory=15)
kk <- enrichKEGG(gene         = convgenes$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
dotplot(kk)


bridge_df <- data.frame(final_bridge_genes_lesional)
colnames(bridge_df) <- c("Relevant_Bridge_Genes", "Number_of_shortest_paths_visits")

ggplot(data=bridge_df[1:10,], aes(x=Relevant_Bridge_Genes, y=Number_of_shortest_paths_visits)) +
  geom_bar(stat="identity", fill="skyblue2")+
  # geom_text(aes(label=len), vjust=-0.3, size=3.5)+
  theme_minimal() + theme(axis.text.x = element_text(angle = 45))


#### final ranks relevant genes and intermediate genes ####
diff_centrality_rank <- read.csv("/path/to/centrality/ranks/rank_comparison_pso_psonl.csv")
diff_centrality_rank <- diff_centrality_rank[order(diff_centrality_rank$MEDIAN.Abs.Difference, decreasing = TRUE),]

relevant_genes_inform <- sort(named_ranklist[which(names(named_ranklist) %in% rownames(relevant_score_mat_mapped))], decreasing = FALSE)
relevant_genes_info_score <- relevant_score_mat_mapped$info_perc
names(relevant_genes_info_score) <- rownames(relevant_score_mat_mapped)
relevant_genes_info_score <- sort(relevant_genes_info_score, decreasing = TRUE)

relevant_genes_centrality_diff <- diff_centrality_rank$Gene[which(diff_centrality_rank$Gene %in% rownames(relevant_score_mat_mapped))]

final_rank_relevant <- TopKLists::Borda(input = list(names(relevant_genes_inform), names(relevant_genes_info_score), relevant_genes_centrality_diff))
final_rank_relevant <- final_rank_relevant$TopK$median


intermediate_genes_inform <- sort(named_ranklist[which(names(named_ranklist) %in% final_bridge_genes)], decreasing = FALSE)
intermediate_genes_centrality_diff <- diff_centrality_rank$Gene[which(diff_centrality_rank$Gene %in% final_bridge_genes)]
final_rank_intermediate <- TopKLists::Borda(input = list(final_bridge_genes, names(intermediate_genes_inform), intermediate_genes_centrality_diff))
final_rank_intermediate <- final_rank_intermediate$TopK$median

###
relevant_genes_inform_nl <- sort(named_ranklist_nl[which(names(named_ranklist_nl) %in% rownames(relevant_score_mat_mapped))], decreasing = FALSE)
final_rank_relevant_nl <- TopKLists::Borda(input = list(names(relevant_genes_inform_nl), names(relevant_genes_info_score), relevant_genes_centrality_diff))
final_rank_relevant_nl <- final_rank_relevant_nl$TopK$median

intermediate_genes_inform_nl <- sort(named_ranklist_nl[which(names(named_ranklist_nl) %in% final_bridge_genes)], decreasing = FALSE)


#### GSEA modules against ranked relevant genes ####
source("/path/to/signed_gsea/gsea_ks.R")
source("/path/to/signed_gsea/signed_ks_test.R")


# Relevant genes
members_list_gsea <- list()
for(i in 1:length(members_list)){
  members_list_gsea[[i]]=members_list[[i]][which(members_list[[i]] %in% final_rank_relevant)]
  
}
names(members_list_gsea) <- paste0("mod_", 1:length(members_list_gsea))

idx <- which(unlist(lapply(members_list_gsea, length))<=2)
members_list_gsea <- members_list_gsea[-idx]
score_list <- list(rank=length(final_rank_relevant):1)
names(score_list$rank) <- final_rank_relevant
res_gsea_fingerprint_relevant = run_gsea_ks_on_multiple_ranked_lists(score_list = score_list, gsets = members_list_gsea, rank.type="one.tail")
res_gsea_fingerprint_relevant$PVAL


# Bridge genes

members_list_gsea_bg <- list()
for(i in 1:length(members_list)){
  members_list_gsea_bg[[i]]=members_list[[i]][which(members_list[[i]] %in% final_bridge_genes)]
  
}
names(members_list_gsea_bg) <- paste0("mod_", 1:length(members_list_gsea_bg))

idx <- which(unlist(lapply(members_list_gsea_bg, length))<=2)
members_list_gsea_bg <- members_list_gsea_bg[-idx]
score_list <- list(rank=length(final_bridge_genes):1)
names(score_list$rank) <- final_bridge_genes
res_gsea_fingerprint_bridge = run_gsea_ks_on_multiple_ranked_lists(score_list = score_list, gsets = members_list_gsea_bg, rank.type="one.tail")
res_gsea_fingerprint_bridge$PVAL

### GSEA only on differentially central genes rank ###

members_list_gsea_cr <- list()
for(i in 1:length(members_list)){
  members_list_gsea_cr[[i]]=members_list[[i]][which(members_list[[i]] %in% diff_centrality_rank$Gene)]
  
}
names(members_list_gsea_cr) <- paste0("mod_", 1:length(members_list_gsea_cr))

idx <- which(unlist(lapply(members_list_gsea_cr, length))<=2)
members_list_gsea_cr <- members_list_gsea_cr[-idx]
score_list <- list(rank=length(diff_centrality_rank$Gene):1)
names(score_list$rank) <- diff_centrality_rank$Gene
res_gsea_fingerprint_cr = run_gsea_ks_on_multiple_ranked_lists(score_list = score_list, gsets = members_list_gsea_cr, rank.type="one.tail")
res_gsea_fingerprint_cr$PVAL


# Human Protein atlas bulkRNA-Seq

# t_reg <- read_delim("/path/to/human_protein_atlas_immune_cells/bulk_rnaseq/blood_cell_category_rna_T-reg_Cell.tsv", 
#                     "\t", escape_double = FALSE, trim_ws = TRUE)
# t_reg <- t_reg[order(t_reg$`RNA single cell type specificity score`, decreasing = TRUE),]
# gdT <- read.table("/path/to/human_protein_atlas_immune_cells/bulk_rnaseq/blood_cell_category_rna_gdT-cell_Cell.tsv", header = TRUE, sep = "\t", quote = "")
# gdT <- gdT[order(gdT$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]
# MAIT <- read.table("/path/to/human_protein_atlas_immune_cells/bulk_rnaseq/blood_cell_category_rna_MAIT.tsv", header = TRUE, sep = "\t", quote = "")
# MAIT <- MAIT[order(MAIT$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]
# naive_B_cells <- read.table("/path/to/human_protein_atlas_immune_cells/bulk_rnaseq/blood_cell_category_rna_naive_B_cells.tsv", header = TRUE, sep = "\t", quote = "")
# naive_B_cells <- naive_B_cells[order(naive_B_cells$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]
# memory_B_cells <- read.table("/path/to/human_protein_atlas_immune_cells/bulk_rnaseq/blood_cell_category_rna_memory_B_cells.tsv", header = TRUE, sep = "\t", quote = "")
# memory_B_cells <- memory_B_cells[order(memory_B_cells$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]
# memory_CD4_T_cells <- read.table("/path/to/human_protein_atlas_immune_cells/bulk_rnaseq/blood_cell_category_rna_memory-CD4_T_cells.tsv", header = TRUE, sep = "\t", quote = "")
# memory_CD4_T_cells <- memory_CD4_T_cells[order(memory_CD4_T_cells$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]
# naive_CD4_T_cells <- read.table("/path/to/human_protein_atlas_immune_cells/bulk_rnaseq/blood_cell_category_rna_naive-CD4_T_cells.tsv", header = TRUE, sep = "\t", quote = "")
# naive_CD4_T_cells <- naive_CD4_T_cells[order(naive_CD4_T_cells$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]
# memory_CD8_T_cells <- read.table("/path/to/human_protein_atlas_immune_cells/bulk_rnaseq/blood_cell_category_rna_memory_CD8_T_cell.tsv", header = TRUE, sep = "\t", quote = "")
# memory_CD8_T_cells <- memory_CD8_T_cells[order(memory_CD8_T_cells$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]
# naive_CD8_T_cells <- read.table("/path/to/human_protein_atlas_immune_cells/bulk_rnaseq/blood_cell_category_rna_naive-CD8_T_cell.tsv", header = TRUE, sep = "\t", quote = "")
# naive_CD8_T_cells <- naive_CD8_T_cells[order(naive_CD8_T_cells$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]
# NK_cells <- read.table("/path/to/human_protein_atlas_immune_cells/bulk_rnaseq/blood_cell_category_rna_NK-cell_Cell.tsv", header = TRUE, sep = "\t", quote = "")
# NK_cells <- NK_cells[order(NK_cells$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]
# classical_monocytes <- read.table("/path/to/human_protein_atlas_immune_cells/bulk_rnaseq/blood_cell_category_rna_classical_monocytes.tsv", header = TRUE, sep = "\t", quote = "")
# classical_monocytes <- classical_monocytes[order(classical_monocytes$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]
# non_classical_monocytes <- read.table("/path/to/human_protein_atlas_immune_cells/bulk_rnaseq/blood_cell_category_rna_non-classical_monocytes.tsv", header = TRUE, sep = "\t", quote = "")
# non_classical_monocytes <- non_classical_monocytes[order(non_classical_monocytes$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]
# intermediate_monocytes <- read.table("/path/to/human_protein_atlas_immune_cells/bulk_rnaseq/blood_cell_category_rna_intermediate_monocytes.tsv", header = TRUE, sep = "\t", quote = "")
# intermediate_monocytes <- intermediate_monocytes[order(intermediate_monocytes$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]
# basophil <- read.table("/path/to/human_protein_atlas_immune_cells/bulk_rnaseq/blood_cell_category_rna_basophil_Cell.tsv", header = TRUE, sep = "\t", quote = "")
# basophil <- basophil[order(basophil$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]
# eosinophil <- read.table("/path/to/human_protein_atlas_immune_cells/bulk_rnaseq/blood_cell_category_rna_eosinophil_Cell.tsv", header = TRUE, sep = "\t", quote = "")
# eosinophil <- eosinophil[order(eosinophil$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]
# neutrophil <- read.table("/path/to/human_protein_atlas_immune_cells/bulk_rnaseq/blood_cell_category_rna_neutrophil_Cell.tsv", header = TRUE, sep = "\t", quote = "")
# neutrophil <- neutrophil[order(neutrophil$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]
# myeloid <- read.table("/path/to/human_protein_atlas_immune_cells/bulk_rnaseq/blood_cell_category_rna_myeloid.tsv", header = TRUE, sep = "\t", quote = "")
# myeloid <- myeloid[order(myeloid$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]
# plasmacytoid <- read.table("/path/to/human_protein_atlas_immune_cells/bulk_rnaseq/blood_cell_category_rna_plasmacytoid.tsv", header = TRUE, sep = "\t", quote = "")
# plasmacytoid <- plasmacytoid[order(plasmacytoid$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]

# cell_list <- list(t_reg, gdT, MAIT, naive_B_cells, memory_B_cells, memory_CD4_T_cells, naive_CD4_T_cells, memory_CD8_T_cells, naive_CD8_T_cells, NK_cells, classical_monocytes,
#                   non_classical_monocytes, intermediate_monocytes, basophil, eosinophil, neutrophil, myeloid, plasmacytoid)
# names(cell_list) <- c("t_reg", "gdT", "MAIT", "naive_B_cells", "memory_B_cells", "memory_CD4_T_cells", "naive_CD4_T_cells", "memory_CD8_T_cells", "naive_CD8_T_cells", "NK_cells", "classical_monocytes",
#                       "non_classical_monocytes", "intermediate_monocytes", "basophil", "eosinophil", "neutrophil", "myeloid", "plasmacytoid")

granulocytes <- read.table("/path/to/human_protein_atlas_immune_cells/cell_type_category_rna_Granulocytes_Cell.tsv", header = TRUE, sep = "\t", quote = "")
granulocytes <- granulocytes[order(granulocytes$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]
hofbauer <- read.table("/path/to/human_protein_atlas_immune_cells/cell_type_category_rna_Hofbauer.tsv", header = TRUE, sep = "\t", quote = "")
hofbauer <- hofbauer[order(hofbauer$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]
macrophages <- read.table("/path/to/human_protein_atlas_immune_cells/cell_type_category_rna_Macrophages_Cell.tsv", header = TRUE, sep = "\t", quote = "")
macrophages <- macrophages[order(macrophages$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]
monocytes <- read.table("/path/to/human_protein_atlas_immune_cells/cell_type_category_rna_Monocytes_Cell.tsv", header = TRUE, sep = "\t", quote = "")
monocytes <- monocytes[order(monocytes$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]
t_cells <- read.table("/path/to/human_protein_atlas_immune_cells/cell_type_category_rna_T-cells_Cell.tsv", header = TRUE, sep = "\t", quote = "")
t_cells <- t_cells[order(t_cells$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]
basal_keratinocytes <- read.table("/path/to/human_protein_atlas_immune_cells/cell_type_category_rna_Basal.tsv", header = TRUE, sep = "\t", quote = "")
basal_keratinocytes <- basal_keratinocytes[order(basal_keratinocytes$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]
suprabasal_keratinocytes <- read.table("/path/to/human_protein_atlas_immune_cells/cell_type_category_rna_Suprabasal.tsv", header = TRUE, sep = "\t", quote = "")
suprabasal_keratinocytes <- suprabasal_keratinocytes[order(suprabasal_keratinocytes$X.RNA.single.cell.type.specificity.score., decreasing = TRUE),]


cell_list <- list(granulocytes, hofbauer, macrophages, monocytes, t_cells, basal_keratinocytes, suprabasal_keratinocytes)
names(cell_list) <- c("granulocytes", "hofbauer", "macrophages", "monocytes", "t_cells", "basal_keratinocytes", "suprabasal_keratinocytes")

members_list_gsea_cells <- list()

for(k in 1:length(cell_list)){
  members_list_gsea_cells[[k]] <- list()
  for(i in 1:length(members_list)){
    members_list_gsea_cells[[k]][[i]]=members_list[[i]][which(members_list[[i]] %in% cell_list[[k]]$Gene)]
  }
  names(members_list_gsea_cells[[k]]) <- paste0("mod_", paste0(1:length(members_list)))
}

names(members_list_gsea_cells) <- c("granulocytes", "hofbauer", "macrophages", "monocytes", "t_cells", "basal_keratinocytes", "suprabasal_keratinocytes")


for(k in 1:length(members_list_gsea_cells)){
  idx <- which(unlist(lapply(members_list_gsea_cells[[k]], length))<=2)
  members_list_gsea_cells[[k]] <- members_list_gsea_cells[[k]][-idx]
  score_list <- list(rank=length(cell_list[[k]]$Gene):1)
  names(score_list$rank) <- cell_list[[k]]$Gene
  res_gsea_fingerprint_cells = run_gsea_ks_on_multiple_ranked_lists(score_list = score_list, gsets = members_list_gsea_cells[[k]], rank.type="one.tail")
  print(names(members_list_gsea_cells[k]))
  print(res_gsea_fingerprint_cells$PVAL)
  
}



### calculate module dissimilarities betweeen disease and control ###
### nodes ###
mod_distances <- matrix(data = NA, nrow = length(conGraph.modules), ncol = length(conGraph_nl.modules), dimnames = list(paste0("lesional_mod_", seq(1:length(conGraph.modules))), paste0("non_lesional_mod_", seq(1:length(conGraph_nl.modules)))))

for(mod in 1:length(conGraph.modules)){
  for(hal in 1:length(conGraph_nl.modules)) {
    
    ji <- length(intersect(conGraph.modules[[mod]], conGraph_nl.modules[[hal]]))/length(union(conGraph.modules[[mod]], conGraph_nl.modules[[hal]]))
    mod_distances[mod, hal] <- ji
    
  }
}


cols <- brewer.pal(3, "YlOrRd")
my_palette <- colorRampPalette(cols)

heatmap.2(t(mod_distances[-c(10,12,13), -c(9,10)]),
          #density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins=c(6,12),     # widens margins around plot
          colsep = c(24,109),    # enable color transition at specified limits
          dendrogram="none",     # only draw a row dendrogram
          keysize =  1 ,  #alter key size
          Colv = FALSE ,      #turn off column clustering
          Rowv =  FALSE ,    # turn off row clustering
          col = my_palette,
          key = TRUE,
          key.xlab = paste("") , #add label to key 
          cexRow = (1.3) , # alter row label font size
          cexCol = (1.3) , # alter column label font size
          notecex = (1.5) , # Alter cell font size
          lmat = rbind(c(0, 3, 0), c(2, 1, 0), c(0, 4, 0) ) , 
          lhei = c(2, 4, 0.9) , # Alter dimensions of display array cell heighs
          lwid = c(0.1, 4, 0.1)) #tweak specific key paramters

### Edges ###
mod_distances_edges <- matrix(data = 0, nrow = length(conGraph.modules), ncol = length(conGraph_nl.modules), dimnames = list(paste0("lesional_mod_", seq(1:length(conGraph.modules))), paste0("non_lesional_mod_", seq(1:length(conGraph_nl.modules)))))

for(mod in 1:length(conGraph.modules)){
  if(length(conGraph.modules[[mod]])<=2){next}
  ind_subgraph_pso <- induced.subgraph(conGraph, vids = conGraph.modules[[mod]])
  edgelist <- as_ids(E(ind_subgraph_pso))
  
  for(hal in 1:length(conGraph_nl.modules)) {
    if(length(conGraph_nl.modules[[hal]])<=2){next}
    ind_subgraph_nl <- induced.subgraph(conGraph_nl, vids = conGraph_nl.modules[[hal]])
    edgelist_nl <- as_ids(E(ind_subgraph_nl))
    ji <- length(intersect(edgelist, edgelist_nl))/length(union(edgelist, edgelist_nl))
    mod_distances_edges[mod, hal] <- ji
    
  }
}

cols <- brewer.pal(3, "YlOrRd")
my_palette <- colorRampPalette(cols)

final_mat <- (mod_distances + mod_distances_edges)/2

heatmap.2(t(final_mat[-c(10,12,13), -10]),
          #density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins=c(6,12),     # widens margins around plot
          colsep = c(24,109),    # enable color transition at specified limits
          dendrogram="none",     # only draw a row dendrogram
          keysize =  1 ,  #alter key size
          Colv = FALSE ,      #turn off column clustering
          Rowv =  FALSE ,    # turn off row clustering
          col = my_palette,
          key = TRUE,
          key.xlab = paste("") , #add label to key 
          cexRow = (1.3) , # alter row label font size
          cexCol = (1.3) , # alter column label font size
          notecex = (1.5) , # Alter cell font size
          lmat = rbind(c(0, 3, 0), c(2, 1, 0), c(0, 4, 0) ) , 
          lhei = c(2, 4, 0.9) , # Alter dimensions of display array cell heighs
          lwid = c(0.1, 4, 0.1)) #tweak specific key paramters

### Fisher test among modules ###

modules_fisher_test <- list()
for (mod in 1:length(members_list)){
  modules_fisher_test[[mod]] <- bc3net::enrichment(genes = members_list[[mod]], genesets = members_list_nl, reference = as_ids(V(conGraph)), adj = "fdr")
}
names(modules_fisher_test) <- paste0("lesional_mod_", 1:length(members_list))


### Pathway analysis only on interesting modules ###
### Pathway on selected modules
lst <- list()
#for(mod in c(2, 3, 5, 6,7)){
for(mod in 1:length(table(members))){
  x <- NA
  x <- names(members[members==mod])
  if(length(x)>10) {
    convgenes = clusterProfiler::bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    print(head(convgenes))
    lst[[mod]] <- convgenes$ENTREZID
  }
}

names(lst) <- seq_along(lst)
lst[sapply(lst, is.null)] <- NULL

res <- compareCluster(lst, fun="enrichPathway")
print(clusterProfiler::dotplot(res))

### Pathway on unselected modules
lst <- list()
for(mod in c(1,4,8, 9, 11)){
  x <- NA
  x <- names(members[members==mod])
  if(length(x)>10) {
    convgenes = clusterProfiler::bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    print(head(convgenes))
    lst[[mod]] <- convgenes$ENTREZID
  }
}

names(lst) <- seq_along(lst)
lst[sapply(lst, is.null)] <- NULL

res <- compareCluster(lst, fun="enrichPathway")
print(clusterProfiler::dotplot(res))


#### Local differences ####

graph_diff <- igraph::difference(conGraph_nl, conGraph)
length(E(graph_diff))
graph_inter <- igraph::intersection(conGraph, conGraph_nl)
length(E(graph_inter))

nodes_diff <- data.frame(lesional=igraph::degree(conGraph, v = V(conGraph)), non_lesional=igraph::degree(conGraph_nl, v = V(conGraph_nl)))
rownames(nodes_diff) <- as_ids(V(conGraph))
nodes_diff$difference <- abs(nodes_diff$lesional-nodes_diff$non_lesional)

betweenness_diff <- betweenness(conGraph, v = V(conGraph), directed = FALSE, normalized = FALSE)
betweenness_diff <- data.frame(lesional=betweenness(conGraph, v = V(conGraph), directed = FALSE, normalized = FALSE), non_lesional=betweenness(conGraph_nl, v = V(conGraph_nl), directed = FALSE, normalized = FALSE))
rownames(betweenness_diff) <- as_ids(V(conGraph_nl))
betweenness_diff$difference <- abs(betweenness_diff$lesional-betweenness_diff$non_lesional)

edgelist_diff <- igraph::as_edgelist(igraph::difference(non_lesional_mod_1, lesional_mod_3), names = TRUE)


### Final modules ranks ###
final_global_rank <- TopKLists::Borda(list(names(named_ranklist), diff_centrality_rank$Gene))
final_global_rank <- final_global_rank$TopK$median

global_df <- data.frame(final_global_rank)
global_df$module <- NA
for(i in 1:length(members_list)){
  global_df$module[which(global_df$final_global_rank %in% members_list[[i]])] <- names(members_list)[i]
  
}

global_df$relevant <- 0
global_df$relevant[which(global_df$final_global_rank %in% rownames(relevant_score_mat_mapped))] <- 1
global_df$bridge <- 0
global_df$bridge[which(global_df$final_global_rank %in% names(intermediate_genes))] <- 1


### Print files  ####
# saveRDS(relevant_score_mat_mapped, file = "relevant_genes_info_score.rds")
# saveRDS(global_df, file = "relevant_and_bridge_on_modules.rds")
# relevant_ranks <- data.frame(inform=names(relevant_genes_inform), info_score=names(relevant_genes_info_score), centrality_diff=relevant_genes_centrality_diff)
# bridge_ranks <- data.frame(shortpath=names(intermediate_genes), inform=names(intermediate_genes_inform), centrality_diff=intermediate_genes_centrality_diff)
# 
# saveRDS(relevant_ranks, file = "relevant_genes_ranks.rds")
# saveRDS(bridge_ranks, file = "bridge_genes_ranks.rds")

