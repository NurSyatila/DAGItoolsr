#' Gene search in MedGen
#'
#' Get disease-associated genes (Entrez Gene identifiers) from MedGen (via Rentrez)
#' @param queryTerms A vector of disease name and/or keywords
#' @return A list of protein-coding genes
#' @return A CSV-formmated file containing results from MedGen
#' @return A text-formatted file containing a list of genes (MedGen_Genes.txt)
#' @return A text-formatted file containing a list of protein-coding genes (MedGen_PCGenes.txt)
#' @examples
#' queryTerms <- c("endometriosis");
#' results <- get_medgen(queryTerms);
#' @export
#' @importFrom utils head write.csv
get_medgen <- function(queryTerms) {
   gene_file <- "MedGen_Genes.txt"
   pcgene_file <- "MedGen_PCGenes.txt"
   new_query_terms <- character()
    for (var in queryTerms){
        query_term <- paste('"',var,'"[TITL]',sep='')
        new_query_terms <- c(new_query_terms, query_term)
    }
    new_query_terms <- paste(new_query_terms, collapse=' OR ')
    qTerms <- paste('(',new_query_terms,')',sep='')
    dbSearch <- suppressMessages(rentrez::entrez_search(db='medgen',term=qTerms,use_history=TRUE,retmax=10000))
   if (length(dbSearch$ids)>0){
       getData <- function(){
           tryCatch({
               print ("Processing (MedGen)...")
               summarylist <- list()
               genelist <- list()
               if (length(dbSearch$ids)>250){
                    chunks <- split(dbSearch$ids, ceiling(seq_along(dbSearch$ids)/250))
                    #loop over
                    for (i in 1:length(chunks)){
                      summary <- suppressMessages(rentrez::entrez_summary(db="medgen",id=unlist(chunks[1],recursive=TRUE)))
                      esummaries <-suppressMessages(rentrez::extract_from_esummary(summary, c("uid","semanticid","semantictype","conceptid","title","definition")))
                      for (i in 1:length(esummaries[1,])){
                          df <- data.frame(paste(esummaries[,i]$uid,collapse=', '),paste(esummaries[,i]$semanticid,collapse=', '),paste(esummaries[,i]$semantictype$value,collapse=', '),paste(esummaries[,i]$conceptid,collapse=', '),paste(esummaries[,i]$title,collapse=', '),paste(esummaries[,i]$definition$value,collapse=', '))
                         summarylist[[paste(esummaries[,i]$uid,collapse=', ')]] <- df
                      }
                    }
                    dbLinks <- suppressMessages(rentrez::entrez_link(dbfrom="medgen", db="gene", cmd="neighbor_history", web_history=dbSearch$web_history))
                    if( !is.null(dbLinks$web_histories$medgen_gene)) {
                        dbFetch <- suppressMessages(rentrez::entrez_fetch(db="gene", web_history=dbLinks$web_histories$medgen_gene,rettype="uilist",retmode="text"))
                        dbList <- head(unlist(strsplit(dbFetch,"\n"),recursive=TRUE),-1)
                    } else { dbList <- list() }
               } else {
                    summary <- suppressMessages(rentrez::entrez_summary(db="medgen",id=unlist(dbSearch$ids,recursive=TRUE)))
                    esummaries <-suppressMessages(rentrez::extract_from_esummary(summary, c("uid","semanticid","semantictype","conceptid","title","definition")))
                    if (length(unlist(dbSearch$ids,recursive=TRUE))>1){
                      for (i in 1:length(esummaries[1,])){
                          df <- data.frame(paste(esummaries[,i]$uid,collapse=', '),paste(esummaries[,i]$semanticid,collapse=', '),paste(esummaries[,i]$semantictype$value,collapse=', '),paste(esummaries[,i]$conceptid,collapse=', '),paste(esummaries[,i]$title,collapse=', '),paste(esummaries[,i]$definition$value,collapse=', '))
                         summarylist[[paste(esummaries[,i]$uid,collapse=', ')]] <- df
                      }
                    } else {
                       df <- data.frame(paste(esummaries$uid,collapse=', '),paste(esummaries$semanticid,collapse=', '),paste(esummaries$semantictype$value,collapse=', '),paste(esummaries$conceptid,collapse=', '),paste(esummaries$title,collapse=', '),paste(esummaries$definition$value,collapse=', '))
                         summarylist[[paste(esummaries$uid,collapse=', ')]] <- df
                    }
                    dbLinks <- suppressMessages(rentrez::entrez_link(dbfrom="medgen", db="gene", id=dbSearch$ids))
                    if( !is.null(dbLinks$links$medgen_gene)) {
                        dbFetch <- suppressMessages(rentrez::entrez_fetch(db="gene", id=dbLinks$links$medgen_gene,rettype="uilist",retmode="text"))
                        dbList <- head(unlist(strsplit(dbFetch,"\n"),recursive=TRUE),-1)
                    } else { dbList <- list() }
                }
                summarytable <- do.call(rbind,summarylist)
                colnames(summarytable) <- c("UID","SemanticID","SemanticType","ConceptID","Title","Definition")
                write.csv(summarytable, file="MedGenSummary.csv",row.names = FALSE)
                #Get Protein Coding Gene IDs
               if (length(dbList)>0){ MappedGenesx <- get_PCGenes(dbList, gene_file,pcgene_file,"ENTREZID") }
               else {
                 print ("No results from MedGen.")
                 MappedGenesx <- list() }
            }, error=function(e) {
             print ("Error: Please check the output files.")
             MappedGenesx <- list() })
            }
            MappedGenesx <- getData()
   } else {
     print ("No results from MedGen.")
     MappedGenesx <- list()
    }
   return(MappedGenesx)
}