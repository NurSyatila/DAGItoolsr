#' Gene search in OMIM
#'
#' Get disease-associated genes (Entrez Gene identifiers) from OMIM (via rentrez)
#' @param queryTerms A vector of disease name and/or keywords
#' @return A list of protein-coding genes
#' @return A CSV-formmated file containing results from OMIM
#' @return A text-formatted file containing a list of genes (OMIM_Genes.txt)
#' @return A text-formatted file containing a list of protein-coding genes (OMIM_PCGenes.txt)
#' @examples
#' queryTerms <- c("endometriosis");
#' results <- get_omim(queryTerms);
#' @export
#' @importFrom utils head write.csv
get_omim <- function(queryTerms) {
   gene_file <- "OMIM_Genes.txt"
   pcgene_file <- "OMIM_PCGenes.txt"
   new_queryTerms <- character()
    for (var in queryTerms){
        query_term <- paste('"',var,'"[DSDR]',sep='')
        new_queryTerms <- c(new_queryTerms, query_term)
    }
    new_queryTerms <- paste(new_queryTerms, collapse=' OR ')
    qTerms <- paste('(',new_queryTerms,')',sep='')
    dbSearch <- suppressMessages(rentrez::entrez_search(db='omim',term=qTerms,use_history=TRUE,retmax=10000))
    if (length(dbSearch$ids)>0){
       getData <- function(){
           tryCatch({
               print ("Processing (OMIM)...")
               summarylist <- list()
               genelist <- list()
               if (length(dbSearch$ids)>250){
                    chunks <- split(dbSearch$ids, ceiling(seq_along(dbSearch$ids)/250))
                    #loop over
                    for (i in 1:length(chunks)){
                      summary <- suppressMessages(rentrez::entrez_summary(db="omim",id=unlist(chunks[1],recursive=TRUE)))
                      esummaries <- suppressMessages(rentrez::extract_from_esummary(summary, c("uid","title","alttitles")))
                      for (i in 1:length(esummaries[1,])){
                          df <- data.frame(paste(esummaries[,i]$uid,collapse=', '),paste(esummaries[,i]$title,collapse=', '),paste(esummaries[,i]$alttitles,collapse=', '))
                         summarylist[[paste(esummaries[,i]$uid,collapse=', ')]] <- df
                      }
                    }
                    dbLinks <- suppressMessages(rentrez::entrez_link(dbfrom="omim", db="gene", cmd="neighbor_history", web_history=dbSearch$web_history))
                    if( !is.null(dbLinks$web_histories$omim_gene)) {
                      dbFetch <- suppressMessages(rentrez::entrez_fetch(db="gene", web_history=dbLinks$web_histories$omim_gene,rettype="uilist",retmode="text"))
                      dbList <- head(unlist(strsplit(dbFetch,"\n"),recursive=TRUE),-1)
                    } else {dbList <- list()}
               } else {
                    summary <- suppressMessages(rentrez::entrez_summary(db="omim",id=unlist(dbSearch$ids,recursive=TRUE)))
                    esummaries <- suppressMessages(rentrez::extract_from_esummary(summary, c("uid","title","alttitles")))
                    if (length(unlist(dbSearch$ids,recursive=TRUE))>1){
                      for (i in 1:length(esummaries[1,])){
                          df <- data.frame(paste(esummaries[,i]$uid,collapse=', '),paste(esummaries[,i]$title,collapse=', '),paste(esummaries[,i]$alttitles,collapse=', '))
                          summarylist[[paste(esummaries[,i]$uid,collapse=', ')]] <- df
                      }
                    } else {
                        df <- data.frame(paste(esummaries$uid,collapse=', '),paste(esummaries$title,collapse=', '),paste(esummaries$alttitles,collapse=', '))
                        summarylist[[paste(esummaries$uid,collapse=', ')]] <- df
                    }
                    dbLinks <- suppressMessages(rentrez::entrez_link(dbfrom="omim", db="gene", id=dbSearch$ids))
                    if( !is.null(dbLinks$links$omim_gene)) {
                      dbFetch <- suppressMessages(rentrez::entrez_fetch(db="gene", id=dbLinks$links$omim_gene,rettype="uilist",retmode="text"))
                      dbList <- head(unlist(strsplit(dbFetch,"\n"),recursive=TRUE),-1)
                    } else {dbList <- list()}
                }
                summarytable <- do.call(rbind,summarylist)
                colnames(summarytable) <- c("Accession","Title","AltTitles")
                write.csv(summarytable, file="OMIMSummary.csv",row.names = FALSE)
                #Get Protein Coding Gene IDs
               if (length(dbList)>0){ MappedGenesx <- get_PCGenes(dbList, gene_file,pcgene_file,"ENTREZID") }
               else {
                 print ("No results from OMIM.")
                 MappedGenesx <- list() }
            }, error=function(e) {
             print ("Error: Please check the output files.")
             MappedGenesx <- list() })
       }
            MappedGenesx <- getData()
   } else { MappedGenesx <- list() }
   return(MappedGenesx)
}