#' Single Differential Expression Worflow
#'
#' Get differentially expressed genes from gene expression data associated to a disease, perform PPI and functional enrichment analyses for a single GEO Series
#' @param queryTerms A vector of disease name and/or keywords to be used for sample classification
#' @param GSEaccession GEO Series (GSE) accession
#' @param GSEplatform GEO Series (GSE) platform
#' @return GeneExpression: A CSV-formmated file containing results of differentially expressed protein-coding genes (FilteredDE.csv) for individual GEO Series in a separate directory
#' @return GeneExpression: A text file containing Entrez identifiers of differentially expressed protein-coding genes (DEGenes.txt) for individual GEO Series in a separate directory
#' @return GeneExpression: A PDF-formatted file containing visualization plots (Plots.pdf) for individual GEO Series in a separate directory
#' @return GeneExpression: A text file containing top 20 genes in tabulated form (Top20DEGenes.txt) for individual GEO Series in a separate directory
#' @return PPI: A summary of PPI analysis (PPISummary.txt)
#' @return PPI: A list of filtered genes from PPI where isolated nodes are discarded (FilteredGenes.txt)
#' @return PPI: Top genes from PPI with highest interaction (Top_PPI_genes.txt)
#' @return PPI: A list of gene clusters with details on enriched KEGG Pathway (GeneCluster.txt)
#' @return PPI: A PDF file containing visualization plots from PPI analysis (STRING_PPI_Network.pdf)
#' @return Functional Enrichment: A separate directory of functional enrichment results for filtered genes (FilteredGenes/) and top genes (TopGenes) containing visualization plots (.pdf) and enrichment results (.csv)
#' @examples
#' queryTerms <- c("endometriosis");
#' GSEaccession <- "GSE23339";
#' GSEplatform <- "GPL6102";
#' results <- differentialexpression_single_workflow(queryTerms,GSEaccession,GSEplatform);
#' @export
#' @importFrom stats na.omit
differentialexpression_single_workflow <- function(queryTerms, GSEaccession,GSEplatform){
  # Example: run the code for a GSE dataset
  print (GSEaccession)
  dir.create(GSEaccession)
  setwd(GSEaccession)
  results_df <- analyse_deg(queryTerms,GSEaccession,GSEplatform)
  if(nrow(results_df) > 0){
    gselist <- results_df[,3]
    names(gselist) <- results_df[,7]
    gselist2 <- na.omit(gselist)
    gselist2 <- gselist2[duplicated(names(gselist2))]
    geneList <- sort(gselist2, decreasing = TRUE)
    alldegenes <- unique(unlist(results_df[,7],recursive=TRUE))
    if (length(alldegenes)<400){
      dir.create("PPI")
      setwd("PPI")
      analyse_ppi_network(alldegenes)
      setwd("../")
    } else {
      print ("Number of genes > 400. PPI analysis cannot be performed. Please use analyse_ppi_network(). Running enrichment analysis instead.")
      GSEnames <- c("BP","MF","CC","KEGG","DO")
      for (GSEname in GSEnames){
        print (paste("**", GSEname, "**", sep=" "))
        suppressMessages(get_enrichment_results(geneList,GSEname,"GSE"))
      }
    }
  }
  setwd ("../")
}
