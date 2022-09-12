#' Get gene modules
#'
#' Get gene modules (clusters of genes according to interactions) and analyse individual gene clusters/modules (required by analyse_ppi_network())
#' @param geneNetwork An igraph object
#' @return PPI: A list of gene clusters with details on enriched KEGG Pathway (GeneCluster.txt)
#' @return A directory named PPI/GeneClusters containing the list of genes in each gene module
#' @examples
#' geneList <- readLines(gene_file)
#' genelist <- as.data.frame(geneList)
#' mapped_genes <- suppressWarnings(string_db$map(genelist, "geneList", removeUnmappedRows = TRUE, quiet = TRUE ))
#' string_db <- STRINGdb$new( version="11.5", species=9606, score_threshold=700, input_directory=".")
#' gene_network <- suppressMessages(string_db$get_subnetwork(mapped_genes$STRING_id))
#' gene_network2 <- igraph::delete.vertices(gene_network, degree(gene_network)==0)
#' gene_network2 <- igraph::simplify(gene_network2, remove.multiple=TRUE, remove.loops=TRUE)
#' results <- get_gene_modules(gene_network2)
#' @export
#' @importFrom graphics par
get_gene_modules <- function(geneNetwork){
  summaryFile <- "GeneCluster.txt"
  print ("Get modules from clustering..")
  string_db <- STRINGdb::STRINGdb$new( version="11.5", species=9606, score_threshold=700, input_directory=".")
  g <- geneNetwork
  c <- igraph::walktrap.community(g,modularity=TRUE)
  c_keep_ids <- as.numeric(names(sizes(c)[sizes(c) >= 4]))
  c_keep_v_idxs <- which(c$membership %in% c_keep_ids)

  g_sub <- igraph::induced_subgraph(g, igraph::V(g)[c_keep_v_idxs])
  c_sub <- c
  c_sub$names <- c$names[c_keep_v_idxs]
  c_sub$membership <- c$membership[c_keep_v_idxs]
  c_sub$vcount <- length(c_sub$names)
  c_sub$modularity <- igraph::modularity(g_sub, c_sub$membership, igraph::E(g_sub)$weight)
  #sort by edges
  df_levels <- unique(igraph::membership(c_sub))
  df_edges <- list()
  df_genes <- list()
  for (x in seq_along(df_levels)){
	edges <- igraph::gsize(igraph::induced_subgraph(g_sub,unlist(c_sub[x],use.names=FALSE)))
	df_edges[[names(c_sub[x])]] <- edges
    df_genes[[names(c_sub[x])]] <- paste(unlist(c_sub[x],use.names=FALSE),collapse=",")
  }
  df_clusters1 <- data.frame(cluster_id=names(df_edges),edges=unlist(df_edges))
  df_clusters2 <- data.frame(cluster_id=names(df_edges),genes=unlist(df_genes))
  df_clusters <- merge(df_clusters1,df_clusters2,by="cluster_id")
  dfClusters <- df_clusters[order(df_clusters$edges,decreasing=TRUE),]
  if (nrow(dfClusters) > 3) {
    par(mfrow=c(2,2))
    for(i in seq(1:4)){
      suppressMessages(string_db$plot_network(unlist(stringr::str_split(dfClusters$genes[i],","),use.names=FALSE), add_link=FALSE, add_summary=TRUE))
    }
  }
  else {
    par(mfrow=c(2,2))
    for(i in seq_along(dfClusters$cluster_id)){
      suppressMessages(string_db$plot_network(unlist(stringr::str_split(dfClusters$genes[i],","),user.names=FALSE), add_link=FALSE, add_summary=TRUE))
    }
  }
  write(paste("\n** Number of clusters in main network: ",length(df_edges)," **",collapse=" "),file=summaryFile,append=TRUE)
  print ("Get pathways for individual clusters..")
  gz_filename <- "9606.protein.info.v11.5.txt.gz"
  annot_read <- readLines(gz_filename)
  annot_read2 <- strsplit(annot_read,'\t')
  annot_table <- lapply(annot_read2, function(x) {x[c(1,2)]})
  string_db <- suppressMessages(STRINGdb::STRINGdb$new( version="11.5", species=9606, score_threshold=700, input_directory="."))
  dir.create("GeneClusters")
  for(i in seq_along(dfClusters$cluster_id)){
      clustergenes <- unlist(stringr::str_split(dfClusters$genes[i],","),use.names=FALSE)
      cluster_aliases <- Filter(function(x) grepl(paste(clustergenes,collapse = "|"), x[1]), annot_table)
      cluster_aliases2 <- unique(unlist(lapply(cluster_aliases, function(x) {x[2]})))
      cluster_entrez <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = cluster_aliases2, c("SYMBOL","ENTREZID"), keytype = "SYMBOL"))
      write(paste("**** CLUSTER ",i," (",length(cluster_aliases2)," genes"," )"," ****",sep=""),file=summaryFile,append=TRUE)
      write(paste(cluster_aliases2,collapse="\n"),file=paste("GeneClusters/cluster_",i,".txt",sep=""))
      write(paste("Gene Symbols: ",paste(cluster_aliases2,collapse=", "),sep=""),file=summaryFile,append=TRUE)
      write(paste("Gene IDs: ",paste(cluster_entrez$ENTREZID,collapse=", "),sep=""),file=summaryFile,append=TRUE)
      ORAresults <- clusterProfiler::enrichKEGG(gene = cluster_entrez$ENTREZID, organism = "hsa", keyType = "ncbi-geneid", pvalueCutoff = 0.01, qvalueCutoff = 0.05, pAdjustMethod = "BH", use_internal_data = FALSE, minGSSize = 3,maxGSSize = 10000)
      sumORA <- as.data.frame(ORAresults)
      titlepic <- "**** Enriched KEGG Pathways ****"
      if (nrow(sumORA)>1){
          if (nrow(sumORA)>4){
            sumORA <- sumORA[1:5,c(1,2,3,6)]
            sumORA$p.adjust <- sprintf("%e",sumORA$p.adjust)
            sumORA2 <- knitr::kable(sumORA,row.names=TRUE,caption=titlepic,"simple")
            write(sumORA2,file=summaryFile,append=TRUE)
            write("\n",file=summaryFile,append=TRUE)
          } else {
            sumORA <- sumORA[,c(1,2,3,6)]
            sumORA$p.adjust <- sprintf("%e",sumORA$p.adjust)
            sumORA2 <- knitr::kable(sumORA,row.names=TRUE,caption=titlepic,"simple")
            write(sumORA2,file=summaryFile,append=TRUE)
            write("\n",file=summaryFile,append=TRUE)
          }
      } else {
          write(paste(titlepic, " - No results",sep=""),file=summaryFile,append=TRUE)
      }
  }
  return (dfClusters)
}
