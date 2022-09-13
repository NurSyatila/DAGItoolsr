#' Get protein coding genes
#'
#' Get protein coding genes from a list of genes
#' @param geneList A vector of genes (identifiers/symbols)
#' @param gene_file A text-formatted file name to store the list of genes in tabular form (ENTREZID, SYMBOL, GENETYPE)
#' @param pcgene_file A text-formatted file name to store the list of protein-coding genes (Entrez identifiers)
#' @param key_type A text-formatted file name to store the list of protein-coding genes (Entrez identifiers)
#' @return Type of geneList, e.g. "ENSEMBL","ENTREZID","SYMBOL","UNIPROT" etc. Refer to keytypes(org.Hs.eg.db)
#' @examples
#' geneList <- c("1588","3586","5241","100048912","10151");
#' gene_file <- "Genes.txt"
#' pcgene_file <- "PCGenes.txt"
#' results <- get_PCGenes(geneList,gene_file,pcgene_file,"ENTREZID");
#' @export
#' @importFrom utils write.table
get_PCGenes <- function(geneList,gene_file,pcgene_file,key_type){
  Genes <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = geneList, c("ENTREZID", "SYMBOL","GENETYPE"), keytype = key_type))
  suppressWarnings(write.table(Genes,file=gene_file, quote=FALSE,sep="\t",row.names = FALSE))
  PCGenes <- dplyr::filter(Genes, grepl("protein",tolower(GENETYPE)))
  PCGenes <- unique(PCGenes$ENTREZID)
  MapGenes <- paste(PCGenes, collapse="\n")
  write(MapGenes,file=pcgene_file)
  return (PCGenes)
}
