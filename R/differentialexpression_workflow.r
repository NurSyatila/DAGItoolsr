#' Differential Expression Worflow
#'
#' Get differentially expressed genes from gene expression data associated to a disease, perform PPI and functional enrichment analyses
#' @param queryTerms A vector of disease name and/or keywords
#' @return A directory named DifferentialExpression
#' @return GeneExpression: A CSV-formmated file containing GEO Datasets associated to a disease (GSESummary.csv)
#' @return GeneExpression: A CSV-formmated file containing GEO Datasets associated to a disease with 'diseased state' condition (DiseaseGSESummary.csv)
#' @return Individual directories for each GEO Series matching the query terms containing result files
#' @return GeneExpression: A CSV-formmated file containing results of differentially expressed protein-coding genes (FilteredDE.csv) for individual GEO Series in a separate directory
#' @return GeneExpression: A text file containing Entrez identifiers of differentially expressed protein-coding genes (DEGenes.txt) for individual GEO Series in a separate directory
#' @return GeneExpression: A PDF-formatted file containing visualization plots (Plots.pdf) for individual GEO Series in a separate directory
#' @return GeneExpression: A text file containing top 20 genes in tabulated form (Top20DEGenes.txt) for individual GEO Series in a separate directory
#' @return A directory named PPI
#' @return PPI: A summary of PPI analysis (PPISummary.txt)
#' @return PPI: A list of filtered genes from PPI where isolated nodes are discarded (FilteredGenes.txt)
#' @return PPI: Top genes from PPI with highest interaction (Top_PPI_genes.txt)
#' @return PPI: A list of gene clusters with details on enriched KEGG Pathway (GeneCluster.txt)
#' @return PPI: A PDF file containing visualization plots from PPI analysis (STRING_PPI_Network.pdf)
#' @return Functional Enrichment: A separate directory of functional enrichment results for filtered genes (FilteredGenes/) and top genes (TopGenes) containing visualization plots (.pdf) and enrichment results (.csv)
#' @examples
#' queryTerms <- c("endometriosis");
#' results <- differentialexpression_workflow(queryTerms);
#' @export
differentialexpression_workflow <- function(queryTerms){
  dir.create("DifferentialExpression")
  setwd("DifferentialExpression")
  GSESum <- get_gse_datasets(queryTerms)
  alldegenes <- list()
  if (length(GSESum[,1])>0){
    for (i in 1:length(GSESum[1,])){
      if (!is.na(GSESum[i,]$gse)==TRUE){
        GSEaccession <- GSESum[i,]$gse
        GSEplatform <- GSESum[i,]$gpl
        print (GSEaccession)
        dir.create(GSEaccession)
        setwd(GSEaccession)
        results_df <- analyse_deg(queryTerms,GSEaccession,GSEplatform)
        if(nrow(results_df) > 0){
          alldegenes <- c(alldegenes,unlist(results_df[,7],recursive=TRUE))
        } else {print ("No results found.")}
        setwd("../")
      }
    }
    if (length(alldegenes)>1){
      alldegenes <- unique(unlist(alldegenes,recursive=TRUE))
      write(alldegenes,file="DEGenes.txt")
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
          suppressMessages(get_enrichment_results(alldegenes,GSEname,"ORA"))
        }
      }
    } else {print ("No results found.")}
  } else {print('No results from GEO DataSets.(1)')}
  setwd ("../")
}
