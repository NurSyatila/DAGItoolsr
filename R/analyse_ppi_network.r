#' Analyse PPI Network
#'
#' A workflow to get predict protein-protein interaction, identify gene modules and perform functional enrichment for a given gene list
#' @param geneList A vector of Entrez gene identifiers
#' @return A directory named PPI
#' @return PPI: A summary of PPI analysis (PPISummary.txt)
#' @return PPI: A PDF file containing visualization plots from PPI analysis (STRING_PPI_Network.pdf)
#' @return PPI: A list of filtered genes from PPI where isolated nodes are discarded (FilteredGenes.txt)
#' @return PPI: Top genes from PPI with highest interaction (Top_PPI_genes.txt)
#' @return PPI: A list of gene clusters with details on enriched KEGG Pathway (GeneCluster.txt)
#' @return A directory named PPI/GeneClusters containing the list of genes in each gene module
#' @return A directory named PPI/FilteredGenes containing functional enrichment results for filtered genes
#' @return A directory named PPI/TopGenes containing functional enrichment results for top 10 or 20 genes
#' @examples
#' gene_file <- 'DatabaseGenes.txt';
#' geneList <- readLines(gene_file)
#' results <- analyse_ppi_network(geneList);
#' @export
#' @importFrom grDevices dev.off pdf rainbow
#' @importFrom utils head write.table
analyse_ppi_network <- function(geneList){
  # get the list of genes
  #geneList <- readLines(gene_file)
  print ("Perform PPI and functional enrichment analysis..")
  summaryFile <- "PPISummary.txt"
  options(download.file.method="curl",timeout = 300)
  # get PPI network
  string_db <- suppressMessages(STRINGdb::STRINGdb$new( version="11.5", species=9606, score_threshold=700, input_directory="."))
  write(paste("** Number of input genes: ",length(geneList),"gene(s) **",collapse=" "),file=summaryFile,append=TRUE)
  genelist <- as.data.frame(geneList)
  mapped_genes <- suppressWarnings(string_db$map(genelist, "geneList", removeUnmappedRows = TRUE, quiet = TRUE ))
  write(paste("** Number of mapped genes (Entrez -> STRINGdb ids): ",length(mapped_genes$STRING_id),"gene(s) **",collapse=" "),file=summaryFile,append=TRUE)
  print ("Get PPI network..")
  pdf('STRING_PPI_Network.pdf')
  suppressMessages(string_db$plot_network( mapped_genes$STRING_id, add_link=TRUE, add_summary=TRUE))

  mapped_genes <- suppressMessages(string_db$map(genelist, "geneList", removeUnmappedRows = TRUE ))
  gz_filename <- "9606.protein.info.v11.5.txt.gz"
  annot_read <- readLines(gz_filename)
  dat_annot <- as.data.frame(do.call(rbind, strsplit(annot_read, split="\t")))
  dat_annot <- dat_annot[,c(1,2)]
  #string_db <- STRINGdb$new( version="11.5", species=9606, score_threshold=700, input_directory=".")
  gene_network <- suppressMessages(string_db$get_subnetwork(mapped_genes$STRING_id))
  # remove singletons and duplicates
  gene_network_new <- igraph::delete.vertices(gene_network, igraph::degree(gene_network)==0)
  gene_network_new2 <- igraph::simplify(gene_network_new, remove.multiple=TRUE, remove.loops=TRUE)
  # get filtered genes (excluding singletons)
  filtered_genes <- as.data.frame(sort(igraph::degree(gene_network_new2),decreasing=TRUE))
  filtered_genes <- cbind(rownames(filtered_genes), data.frame(filtered_genes, row.names=NULL))
  colnames(filtered_genes) <- c("ensembl_id","degree")
  filtered_aliases <- dplyr::filter(dat_annot, grepl(paste(filtered_genes$ensembl_id,collapse = "|"),dat_annot[,1]))
  colnames(filtered_aliases) <- c("ensembl_id","gene_symbol")
  # analyse filtered genes and sort according to the degree (number of interactions with other genes)
  filtered_genes2 <- merge(filtered_genes, filtered_aliases,all = TRUE)
  filtered_entrez <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = filtered_genes2$gene_symbol, c("SYMBOL","ENTREZID"), keytype = "SYMBOL"))
  colnames(filtered_entrez) <- c("gene_symbol","gene_id")
  filtered_genes3 <- merge(filtered_genes2, filtered_entrez,all = TRUE)
  filtered_genes3 <- filtered_genes3[order(filtered_genes3$degree,decreasing=TRUE),]
  # save output in a text file
  write(paste("** Total genes from PPI network (isolated nodes discarded): ", length(filtered_genes3$ensembl_id)," **\n",sep=""),file=summaryFile,append=TRUE)
  suppressMessages(string_db$plot_network( filtered_genes3$ensembl_id, add_link=TRUE, add_summary=TRUE))
  # perform functional enrichment analysis for teh following terms: KEGG, GO (BP, MF, CC)
  print ("Functional enrichment analysis..")
  write(paste("** FilteredGenes: Functional enrichment using clusterProfiler [",length(filtered_genes3$ensembl_id),"gene(s)] (isolated nodes discarded) **\n",collapse=" "),file=summaryFile,append=TRUE)
  suppressWarnings(write.table(filtered_genes3,file="FilteredGenes.txt",quote=FALSE,append=TRUE, col.names = TRUE, row.names = FALSE, sep = "\t"))
  dir.create("FilteredGenes")
  setwd("FilteredGenes")
  ORAnames <- c("BP","MF","CC","KEGG","DO")
  for (ORAname in ORAnames){
    print (paste("**", ORAname, "**", sep=" "))
    suppressMessages(get_enrichment_results(filtered_genes3$gene_id,ORAname,"ORA"))
  }
  setwd("../")
  # get top 20 genes (with highest degree) and perform functional enrichment analysis for teh following terms: KEGG, GO (BP, MF, CC)
  if (length(filtered_genes$ensembl_id)>9){
    print ("Get top 10/20 genes..")
    if (length(filtered_genes$ensembl_id)>19){ top_genes_degrees <- head(as.data.frame(sort(igraph::degree(gene_network_new2),decreasing=TRUE)),20)}
    else {top_genes_degrees <- head(as.data.frame(sort(igraph::degree(gene_network_new2),decreasing=TRUE)),10)}
    top_genes_degrees <- cbind(rownames(top_genes_degrees), data.frame(top_genes_degrees, row.names=NULL))
    colnames(top_genes_degrees) <- c("ensembl_id","degree")
    top_aliases <- dplyr::filter(dat_annot, grepl(paste(top_genes_degrees$ensembl_id,collapse = "|"),dat_annot[,1]))
    colnames(top_aliases) <- c("ensembl_id","gene_symbol")
    top_genes_degrees2 <- merge(top_genes_degrees, top_aliases,all = TRUE)
    top_entrez <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = top_genes_degrees2$gene_symbol, c("SYMBOL","ENTREZID"), keytype = "SYMBOL"))
    colnames(top_entrez) <- c("gene_symbol","gene_id")
    top_genes_degrees3 <- merge(top_genes_degrees2, top_entrez,all = TRUE)
    top_genes_degrees3 <- top_genes_degrees3[order(top_genes_degrees3$degree,decreasing=TRUE),]
    # get the PPI plot for top 20 genes
    gene_network_top <- string_db$get_subnetwork(top_genes_degrees3$ensembl_id)
    geneAnnotation_sorted <- dat_annot[match(igraph::V(gene_network_top)$name, dat_annot$V1),]
    igraph::V(gene_network_top)$name <- geneAnnotation_sorted$V2
    igraph::V(gene_network_top)$vertex.frame.color <- "white"
    edgeweights <- igraph::E(gene_network_top)$weight * 2.0
    suppressWarnings(plot(
      gene_network_top,
      layout=igraph::layout.fruchterman.reingold,
      edge.curved=FALSE,
      vertex.label.color="black",
      asp=FALSE,
      vertex.label.cex=0.6,
      vertex.shape="circle",
      vertex.color=rainbow(igraph::betweenness(gene_network_top), start=0, end=2/6),
      edge.width=edgeweights,
      edge.arrow.mode=0,
      main="PPI Network (Top genes with the highest degree)"))

    write("** Top genes from PPI network (with highest degree) **",file=summaryFile,append=TRUE)
    top_20_genes_table <- knitr::kable(top_genes_degrees3,row.names=TRUE,caption="List of top genes from PPI network","simple")
    suppressMessages(write(top_20_genes_table,file=summaryFile,append=TRUE))
    # perform functional enrichment analysis for top 20 genes
    dir.create("TopGenes")
    setwd("TopGenes")
    write(paste(top_genes_degrees3$gene_symbol,collapse="\n"),file="../Top_PPI_genes.txt")
    suppressMessages(get_enrichment_results(top_genes_degrees3$gene_id,"KEGG","ORA"))
    setwd("../")
  } else {
    write(paste("** Top genes cannot be generated **\n",collapse=" "),file=summaryFile,append=TRUE)
  }
  print ("Get gene modules/clusters..")
  # get gene modules
  suppressMessages(get_gene_modules(gene_network_new2))
  dev.off()
}
