#' Gene search in CLinVar
#'
#' Get disease-associated genes (Entrez Gene identifiers) from ClinVar (via rentrez)
#' @param queryTerms A vector of disease name and/or keywords
#' @return A list of protein-coding genes
#' @return A CSV-formmated file containing results from ClinVar
#' @return A text-formatted file containing a list of genes (ClinVar_Genes.txt)
#' @return A text-formatted file containing a list of protein-coding genes (ClinVar_PCGenes.txt)
#' @examples
#' queryTerms <- c("endometriosis");
#' results <- get_clinvar(queryTerms);
#' @export
#' @importFrom utils head write.csv
get_clinvar <- function(queryTerms) {
   gene_file <- "ClinVar_Genes.txt"
   pcgene_file <- "ClinVar_PCGenes.txt"
   new_query_terms <- character()
   for (var in queryTerms){
       query_term <- paste('"',var,'"[DIS]',sep='')
       new_query_terms <- c(new_query_terms, query_term)
   }
   new_query_terms <- paste(new_query_terms, collapse=' OR ')
   qTerms <- paste('(((',new_query_terms,') AND 9606[TID]) NOT "clinsig conflicts"[FILT]) NOT "clinsig vus"[FILT]',sep='')
   dbSearch <- suppressMessages(rentrez::entrez_search(db='clinvar',term=qTerms,use_history=TRUE,retmax=10000))
   if (length(dbSearch$ids)>0){
       getData <- function(){
           tryCatch({
               print ("Processing (ClinVar)...")
               summarylist <- list()
               genelist <- list()
               if (length(dbSearch$ids)>250){
                    chunks <- split(dbSearch$ids, ceiling(seq_along(dbSearch$ids)/250))
                    #loop over
                    for (i in 1:length(chunks)){
                      summary <- suppressMessages(rentrez::entrez_summary(db="clinvar",id=unlist(chunks[1],recursive=TRUE)))
                      esummaries <- suppressMessages(rentrez::extract_from_esummary(summary, c("uid","accession","title","clinical_significance","genes", "trait_set")))
                      for (i in 1:length(esummaries[1,])){
                          df <- data.frame(paste(esummaries[,i]$uid,collapse=', '),paste(esummaries[,i]$accession,collapse=', '),paste(esummaries[,i]$title,collapse=', '),paste(esummaries[,i]$clinical_significance$description,collapse=', '),paste(esummaries[,i]$genes$geneid,collapse=', '),paste(esummaries[,i]$trait_set$trait_name,collapse=', '))
                          summarylist[[paste(esummaries[,i]$uid,collapse=', ')]] <- df
                      }
                    }
                    dbLinks <- suppressMessages(rentrez::entrez_link(dbfrom="clinvar", db="gene", cmd="neighbor_history", web_history=dbSearch$web_history))
                    if( !is.null(dbLinks$web_histories$clinvar_gene)) {
                      dbFetch <- suppressMessages(rentrez::entrez_fetch(db="gene", web_history=dbLinks$web_histories$clinvar_gene,rettype="uilist",retmode="text"))
                      dbList <- head(unlist(strsplit(dbFetch,"\n"),recursive=TRUE),-1)
                    } else {dbList <- list()}
               } else {
                    summary <- suppressMessages(rentrez::entrez_summary(db="clinvar",id=unlist(dbSearch$ids,recursive=TRUE)))
                    esummaries <- suppressMessages(rentrez::extract_from_esummary(summary, c("uid","accession","title","clinical_significance","genes", "trait_set")))
                    if (length(unlist(dbSearch$ids,recursive=TRUE))>1){
                      for (i in 1:length(esummaries[1,])){
                          df <- data.frame(paste(esummaries[,i]$uid,collapse=', '),paste(esummaries[,i]$accession,collapse=', '),paste(esummaries[,i]$title,collapse=', '),paste(esummaries[,i]$clinical_significance$description,collapse=', '),paste(esummaries[,i]$genes$geneid,collapse=', '),paste(esummaries[,i]$trait_set$trait_name,collapse=', '))
                        summarylist[[paste(esummaries[,i]$uid,collapse=', ')]] <- df
                      }
                    } else {
                       df <- data.frame(paste(esummaries$uid,collapse=', '),paste(esummaries$accession,collapse=', '),paste(esummaries$title,collapse=', '),paste(esummaries$clinical_significance$description,collapse=', '),paste(esummaries$genes$geneid,collapse=', '),paste(esummaries$trait_set$trait_name,collapse=', '))
                        summarylist[[paste(esummaries$uid,collapse=', ')]] <- df
                    }
                    dbLinks <- suppressMessages(rentrez::entrez_link(dbfrom="clinvar", db="gene", id=dbSearch$ids))
                    if( !is.null(dbLinks$links$clinvar_gene)) {
                      dbFetch <- suppressMessages(rentrez::entrez_fetch(db="gene", id=dbLinks$links$clinvar_gene,rettype="uilist",retmode="text"))
                      dbList <- head(unlist(strsplit(dbFetch,"\n"),recursive=TRUE),-1)
                    } else {dbList <- list()}
                }
                summarytable <- do.call(rbind,summarylist)
                colnames(summarytable) <- c("UID","Accession","Title","Clinical Significance","GeneIDs","Traits")
                write.csv(summarytable, file="ClinvarSummary.csv",row.names = FALSE)
                #Get Protein Coding Gene IDs
               if (length(dbList)>0){ MappedGenesx <- get_PCGenes(dbList, gene_file,pcgene_file,"ENTREZID") }
               else {
                 print ("No results from Clinvar.")
                 MappedGenesx <- list() }
            }, error=function(e)
                {
             print ("Error: Please check the output files.")
             MappedGenesx <- list()
           }) }
            MappedGenesx <- getData()
   } else {
       print ("No results from ClinVar.")
       MappedGenesx <- list()
    }
   return(MappedGenesx)
}