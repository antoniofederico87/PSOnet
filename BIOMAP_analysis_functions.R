map_genes_on_modules <- function(geneSet, igraph_modules){
  
  members=igraph::membership(igraph_modules)
  members_list <- list()
  
  for (i in 1:length(table(members))){
    members_list[[i]] <- names(members[members==i])
  }
  names(members_list) <- paste0("module_", 1:length(table(members)))
  DEGmap_on_modules <- list()
  for(i in 1:length(members_list)){
    DEGmap_on_modules[[i]] <- intersect(geneSet, members_list[[i]])
  }
  return(DEGmap_on_modules)
  
}


get_reactome_from_multiple_datasets <- function(topTable_list, geneID="SYMBOL", pval_cutoff=0.05, layout="both") {
  
  if(layout=="overall"){
    pathways_list <- list()
    sigpath <- c()
    #pdf(file = paste0(subDir1, "report_annotazione_funzionale.pdf"), paper = "a4" , height = 1600, width = 900)
    for(mod in 1:length(topTable_list)){
      x <- NA
      x <- rownames(topTable_list)
      if(length(x)>10) {
        eg = clusterProfiler::bitr(x, fromType=geneID, toType="ENTREZID", OrgDb="org.Hs.eg.db")
        print(head(eg))
        sigpath <- enrichPathway(gene=eg$ENTREZID, pvalueCutoff=pval_cutoff, readable=T)
        sigpath <- as.data.frame(sigpath)
        print(head(sigpath))
        sigpath <- sigpath@result[which(sigpath@result$p.adjust<pval_cutoff),]
        dataset <- names(topTable_list)[mod]
        print(dataset)
        sigpath$dataset <- rep(dataset, length(sigpath$ID))
        pathways_list[[data]] <- sigpath
        
      }
      
    }
  }
  return(pathways_list)
}


get_reactome_from_modules <- function(igraph_modules, geneID="SYMBOL", pval_cutoff=0.05, outPath, layout="both") {
  if (file.exists(outPath)){
    setwd(file.path(outPath))
  } else {
    dir.create(file.path(outPath))
    setwd(file.path(outPath))
  }
  
  cat("The files will be exported in ", getwd())
  members=igraph::membership(igraph_modules)
  if(layout=="overall"){
    sigpath.overall <- data.frame()
    sigpath <- c()
    #pdf(file = paste0(subDir1, "report_annotazione_funzionale.pdf"), paper = "a4" , height = 1600, width = 900)
    for(mod in 1:length(igraph_modules)){
      x <- NA
      x <- names(members[members==mod])
      if(length(x)>10) {
        eg = clusterProfiler::bitr(x, fromType=geneID, toType="ENTREZID", OrgDb="org.Hs.eg.db")
        print(head(eg))
        sigpath <- enrichPathway(gene=eg$ENTREZID, pvalueCutoff=pval_cutoff, readable=T)
        sigpath <- as.data.frame(sigpath)
        print(head(sigpath))
        if (length(sigpath$ID)>0){
          sigpath$Module <- mod
          sigpath.overall <- rbind(sigpath.overall, sigpath)
          write.csv(as.data.frame(sigpath), file = paste0("Significant_enriched_pathways_module_", mod, ".csv"), quote = FALSE, row.names = FALSE)
        }
        
      }
      
    }
    write.csv(sigpath.overall, file = "Pathway_results_overall.csv", row.names = FALSE, quote = FALSE)
    
  } else if (layout=="single") {
    sigpath.overall <- data.frame()
    sigpath <- c()
    #pdf(file = paste0(subDir1, "report_annotazione_funzionale.pdf"), paper = "a4" , height = 1600, width = 900)
    for(mod in 1:length(igraph_modules)){
      x <- NA
      x <- names(members[members==mod])
      if(length(x)>10) {
        eg = clusterProfiler::bitr(x, fromType=geneID, toType="ENTREZID", OrgDb="org.Hs.eg.db")
        print(head(eg))
        sigpath <- enrichPathway(gene=eg$ENTREZID, pvalueCutoff=pval_cutoff, readable=T)
        sigpath <- as.data.frame(sigpath)
        print(head(sigpath))
        if (length(sigpath$ID)>0){
          sigpath$Module <- mod
          #sigpath.overall <- rbind(sigpath.overall, sigpath)
          write.csv(as.data.frame(sigpath), file = paste0("Significant_enriched_pathways_module_", mod, ".csv"), quote = FALSE, row.names = FALSE)
        }
        
      }
      
    }
    
  }
  
  return(cat("Analysis completed!!!"))
  
}

get_bubbleplot_from_pathways <- function(igraph_modules, geneID="SYMBOL") {
  
  lst <- list()
  members=igraph::membership(igraph_modules)
  for(mod in 1:length(conGraph.modules)){
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
  return(res)
  
}
