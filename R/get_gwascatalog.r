#' Gene search in GWAS Catalog
#'
#' Get disease-associated genes (Entrez Gene identifiers) from GWAS Catalog (via gwasrapidd)
#' @param queryTerms A vector of disease name and/or keywords
#' @return A list of protein-coding genes
#' @return A CSV-formmated file containing results from GWAS Catalog database
#' @return A text-formatted file containing a list of genes (GWASCatalog_Genes.txt)
#' @return A text-formatted file containing a list of protein-coding genes (GWASCatalog_PCGenes.txt)
#' @examples
#' queryTerms <- c("endometriosis");
#' results <- get_gwascatalog(queryTerms);
#' @export
#' @importFrom stats na.omit
#' @importFrom utils write.csv
#' @importFrom magrittr %>%
get_gwascatalog <- function(queryTerms) {
    gene_file <- "GWASCatalog_Genes.txt"
    pcgene_file <- "GWASCatalog_PCGenes.txt"
    gwas_list <- character()
    summarylist <- list()
    print ("Processing (GWAS Catalog)... ")
    for (qt in queryTerms){
        my_associations <- gwasrapidd::get_associations(efo_trait = qt)
        if (length(my_associations@associations$association_id)>0){
            dplyr::filter(my_associations@associations, pvalue < 1e-6) %>% # Filter by p-value
            tidyr::drop_na(pvalue) %>%
            dplyr::pull(association_id) -> association_ids # Extract column association_id
            my_associations2 <- my_associations[association_ids]
            new_associations <- subset(as.data.frame(my_associations2@entrez_ids), !is.na(entrez_id))
            new_associations$pvalue <- my_associations[new_associations$association_id]@associations$pvalue
            new_associations$variant_id <- paste(unique(my_associations[new_associations$association_id]@risk_alleles$variant_id),collapse=",")
            new_associations$queryterm <- qt
            if (nrow(new_associations>0)){
                gwas_list <- c(gwas_list, new_associations$entrez_id)
                summarylist[[qt]] <- new_associations
            }
        }
    }
    if (length(summarylist) > 0){
        summarytable <- unique(do.call(rbind,summarylist))
        colnames(summarytable) <- c("association_id","locus_id","gene_name","entrez_id","pvalue","variant_id","queryterm")
        write.csv(summarytable, file="GWASCatalogSummary.csv",row.names = FALSE)
        gwas_list2 <- unique(na.omit(gwas_list))
        if (length(gwas_list2)>0){
            gwas_protein_coding_genesx <- get_PCGenes(gwas_list2, gene_file,pcgene_file,"ENTREZID")
        } else {
          print ("No results from GWAS Catalog.")
          gwas_protein_coding_genesx <- list() }
    } else {
      gwas_protein_coding_genesx <- list()
      print ("No results from GWAS Catalog.")}
    return(gwas_protein_coding_genesx)
}