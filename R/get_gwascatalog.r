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
            if (length(my_associations2@associations$association_id)>0){
                gwas_geneids <- my_associations2@entrez_ids$entrez_id # Get Entrez IDs
                gwas_list <- c(gwas_list, gwas_geneids)
                for (i in 1:length(my_associations2@associations$association_id)){
                    var1 <- my_associations2[i,]@associations$association_id
                    var2 <- my_associations2[i,]@associations$pvalue
                    var3 <- my_associations2[i,]@risk_alleles$variant_id
                    var4 <- my_associations2[i,]@entrez_ids$gene_name
                    var5 <- my_associations2[i,]@entrez_ids$entrez_id
                    getvar6 <- function(){
                        tryCatch({ var6 <- gwasrapidd::get_traits(association_id = my_associations2[i,]@entrez_ids$association_id)@traits$efo_id }, error=function(e) { var6 <- "NA" })
                    }
                    getvar7 <- function(){
                        tryCatch({ var7 <- gwasrapidd::get_traits(association_id = my_associations2[i,]@entrez_ids$association_id)@traits$trait }, error=function(e) { var6 <- "NA" })
                    }
                    var6 <- getvar6()
                    var7 <- getvar7()
                    df <- data.frame(paste(var1,collapse=', '),paste(var2,collapse=', '),paste(var3,collapse=', '),paste(var4,collapse=', '),paste(var5,collapse=', '),paste(var6,collapse=', '),paste(var7,collapse=', '))
                    summarylist[[paste(var1,collapse=', ')]] <- df
                }
            }
        }
    }
    if (length(summarylist) > 0){
        summarytable <- unique(do.call(rbind,summarylist))
        colnames(summarytable) <- c("AssociationID","pValue","VariantID","GeneName","EntrezID","EFOID","Trait")
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