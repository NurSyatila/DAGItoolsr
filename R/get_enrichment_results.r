#' Get enrichment results
#'
#' Perform enrichment analysis for a gene list (over-representation analysis) or a dataframe consisting gene identifiers and log value (gene set enrichment analysis)
#' @param geneList A vector / dataframe of Entrez gene identifiers
#' @param ORAterm Functional term, BP, MF, CC, KEGG or DO
#' @param enrichmenttype ORA (over-representation analysis) or GSE (gene set enrichment analysis)
#' @return A CSV-formatted file of enrichment results
#' @return A PDF file containing visualization plots from enrichment analysis
#' @examples
#' #gene_file <- 'DatabaseGenes.txt';
#' #geneList <- readLines(gene_file)
#' geneList <- c("1588","3586","5241","7124", "998", "54361", "55591", "283455", "3480","8390","83608","9687","4211");
#' results <- get_enrichment_results(geneList,"KEGG","ORA");
#' @export
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics barplot
#' @importFrom utils data write.csv
get_enrichment_results <- function(geneList,ORAterm,enrichmenttype) {
    if (ORAterm=="BP"){
        if (enrichmenttype=="GSE")
          { ORAresults <- clusterProfiler::gseGO(gene = geneList, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID", nPerm = 10000, ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "none", minGSSize = 3, maxGSSize = 10000)}
        else { ORAresults <- clusterProfiler::enrichGO(gene = geneList, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.01, qvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 3, maxGSSize = 10000)}
    }
    if (ORAterm=="MF"){
      if (enrichmenttype=="GSE")
        { ORAresults <- clusterProfiler::gseGO(gene = geneList, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID", nPerm = 10000, ont = "MF", pvalueCutoff = 0.05, pAdjustMethod = "none", minGSSize = 3, maxGSSize = 10000) }
      else { ORAresults <- clusterProfiler::enrichGO(gene = geneList, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID", ont = "MF",pvalueCutoff = 0.01, qvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 3, maxGSSize = 10000) }
    }
    if (ORAterm=="CC"){
      if (enrichmenttype=="GSE")
        { ORAresults <- clusterProfiler::gseGO(gene = geneList, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID", nPerm = 10000, ont = "CC", pvalueCutoff = 0.05, pAdjustMethod = "none", minGSSize = 3, maxGSSize = 10000) }
      else { ORAresults <- clusterProfiler::enrichGO(gene = geneList, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID", ont = "CC", pvalueCutoff = 0.01, qvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 3, maxGSSize = 10000) }
    }
    if (ORAterm=="KEGG"){
      if (enrichmenttype=="GSE")
        { ORAresults <- clusterProfiler::gseKEGG(gene = geneList, organism = "hsa", keyType = "ncbi-geneid", nPerm = 10000, pvalueCutoff = 0.05, pAdjustMethod = "none", use_internal_data = FALSE, minGSSize = 3,maxGSSize = 10000) }
      else { ORAresults <- clusterProfiler::enrichKEGG(gene = geneList, organism = "hsa", keyType = "ncbi-geneid", pvalueCutoff = 0.01, qvalueCutoff = 0.05, pAdjustMethod = "BH", use_internal_data = FALSE, minGSSize = 3,maxGSSize = 10000) }
    }
    if (ORAterm=="DO"){
      if (enrichmenttype=="GSE")
        { ORAresults <- DOSE::gseDO(gene = geneList, nPerm = 10000, pvalueCutoff = 0.05, pAdjustMethod = "none", minGSSize = 3, maxGSSize = 10000) }
      else { ORAresults <- DOSE::enrichDO(gene = geneList, pvalueCutoff = 0.01, qvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 3, maxGSSize = 10000) }
    }
    sumORA <- as.data.frame(ORAresults)
    if (nrow(sumORA)>0){
        write.csv(sumORA , file = paste(ORAterm,"_Results.csv",sep=''), row.names=TRUE)
        pdf(file=paste(ORAterm,"_Plots.pdf",sep=''))
        if (enrichmenttype=="GSE"){
          plot1 <- suppressWarnings(enrichplot::dotplot(ORAresults, showCategory=10, split=".sign") + ggplot2::facet_grid(.~.sign))
          ORAresultsx <- DOSE::setReadable(ORAresults, org.Hs.eg.db::org.Hs.eg.db, 'ENTREZID')
          plot2 <- suppressWarnings(enrichplot::cnetplot(ORAresultsx, categorySize="pvalue", foldChange=geneList,node_label="all"))
          plot3 <- suppressWarnings(enrichplot::heatplot(ORAresultsx,showCategory=5, foldChange=geneList))
          ORAresultsx2 <- enrichplot::pairwise_termsim(ORAresultsx)
          suppressWarnings(options(ggrepel.max.overlaps = Inf))
          plot4 <- suppressWarnings(enrichplot::emapplot(ORAresultsx2))
          print(plot1)
          print(plot2)
          print(plot3)
          print(plot4)
        }
        else {
          plot1 <- suppressWarnings(barplot(ORAresults, showCategory=10))
          plot2 <- suppressWarnings(enrichplot::dotplot(ORAresults, showCategory=10))
          ORAresultsx <- DOSE::setReadable(ORAresults, org.Hs.eg.db::org.Hs.eg.db, 'ENTREZID')
          plot3 <- suppressWarnings(enrichplot::cnetplot(ORAresultsx, categorySize="pvalue",node_label="all"))
          plot4 <- suppressWarnings(enrichplot::heatplot(ORAresultsx,showCategory=5))
          ORAresultsx2 <- enrichplot::pairwise_termsim(ORAresultsx)
          suppressWarnings(options(ggrepel.max.overlaps = Inf))
          plot5 <- suppressWarnings(enrichplot::emapplot(ORAresultsx2))
          print(plot1)
          print(plot2)
          print(plot3)
          print(plot4)
          print(plot5)
        }
        dev.off()
        if (ORAterm=="KEGG"){
            utils::data("bods", package = "pathview")
            utils::data("gene.idtype.bods", package = "pathview")
            print ("pathview")
            enriched_pathways <- sumORA[1:5,1]
            dir.create("Pathview")
            setwd("Pathview")
            for (pth_id in enriched_pathways){
                pathview::pathview(gene.data=geneList, pathway.id=pth_id, species = "hsa")
            }
            setwd("../")

        }
    } else {
        print(paste(ORAterm, " - No results",sep=""))
    }
}
