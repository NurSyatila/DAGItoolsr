#' Get GSE datasets associated to a disease
#'
#' Retrieve GSE Datasets associated to a disease based on MeSH Terms
#' @param queryTerms A vector of disease name and/or keywords
#' @return A dataframe of results (GEO Dataset with 'diseased state' condition) from GEO Dataset web search through Rentrez
#' @return A CSV-formmated file containing GEO Datasets associated to a disease (GSESummary.csv)
#' @return A CSV-formmated file containing GEO Datasets associated to a disease with 'diseased state' condition (DiseaseGSESummary.csv)
#' @examples
#' queryTerms <- c("endometriosis");
#' results <- get_gse_datasets(queryTerms);
#' @export
#' @importFrom utils write.csv
get_gse_datasets <- function(queryTerms) {
    MeSHTerms <- get_meshterms(queryTerms,"mesh_term")
    print (MeSHTerms)
    new_query_terms <- character()
    for (var in MeSHTerms){
        query_term <- paste('"',var,'"[MESH]',sep='')
        new_query_terms <- c(new_query_terms, query_term)
    }
    new_query_terms <- paste(new_query_terms, collapse=' OR ')
    qTerms <- paste('(',new_query_terms,') AND "Homo sapiens"[ORGN] AND "gds"[FILT]',sep='')
    dbSearch <- rentrez::entrez_search(db='gds',term=qTerms,use_history=TRUE,retmax=10000)
    if (length(dbSearch$ids)>0){
        getData <- function(){
           tryCatch({
               print ("Processing (GEO DataSets)...")
               summarylist <- list()
               if (length(dbSearch$ids)>250){
                    chunks <- split(dbSearch$ids, ceiling(seq_along(dbSearch$ids)/250))
                    #loop over
                    for (i in 1:length(chunks)){
                      summary <- suppressMessages(rentrez::entrez_summary(db="gds",id=unlist(chunks[1],recursive=TRUE)))
                      esummaries <- suppressMessages(rentrez::extract_from_esummary(summary, c("accession","gpl","gse","seriestitle","subsetinfo","samples")))
                      for (i in 1:length(esummaries[1,])){
                          df <- data.frame(paste(esummaries[,i]$accession,collapse=', '),paste(esummaries[,i]$gpl,collapse=', '),paste(esummaries[,i]$gse,collapse=', '),paste(esummaries[,i]$seriestitle,collapse=', '),paste(esummaries[,i]$subsetinfo,collapse=', '),paste(esummaries[,i]$samples$title,collapse='||'))
                          summarylist[[paste(esummaries[,i]$accession,collapse=', ')]] <- df
                      }
                    }
               } else {
                    summary <- suppressMessages(rentrez::entrez_summary(db="gds",id=unlist(dbSearch$ids,recursive=TRUE)))
                    esummaries <- suppressMessages(rentrez::extract_from_esummary(summary, c("accession","gpl","gse","seriestitle","subsetinfo","samples")))
                    if (length(unlist(dbSearch$ids,recursive=TRUE))>1){
                      for (i in 1:length(esummaries[1,])){
                          df <- data.frame(paste(esummaries[,i]$accession,collapse=', '),paste(esummaries[,i]$gpl,collapse=', '),paste(esummaries[,i]$gse,collapse=', '),paste(esummaries[,i]$seriestitle,collapse=', '),paste(esummaries[,i]$subsetinfo,collapse=', '),paste(esummaries[,i]$samples$title,collapse='||'))
                          summarylist[[paste(esummaries[,i]$accession,collapse=', ')]] <- df
                      }
                    } else {
                       df <- data.frame(paste(esummaries$accession,collapse=', '),paste(esummaries[,i]$gpl,collapse=', '),paste(esummaries$gse,collapse=', '),paste(esummaries$seriestitle,collapse=', '),paste(esummaries$subsetinfo,collapse=', '),paste(esummaries$samples$title,collapse='||'))
                          summarylist[[paste(esummaries$accession,collapse=', ')]] <- df
                    }
                }
                summarytable <- do.call(rbind,summarylist)
                colnames(summarytable) <- c("accession","gpl","gse","seriestitle","subsetinfo","samples")
                summarytable$gpl <- paste("GPL",summarytable$gpl,sep="")
                summarytable$gse <- paste("GSE",summarytable$gse,sep="")
                write.csv(summarytable, file="GSESummary.csv",row.names = FALSE)
                Filteredsummarytable <- dplyr::filter(summarytable, grepl("disease",tolower(summarytable$subsetinfo)))
                Filteredsummarytable <- Filteredsummarytable[!duplicated(Filteredsummarytable[,"gse"]),]
                write.csv(Filteredsummarytable, file="DiseaseGSESummary.csv",row.names = FALSE)
                return (Filteredsummarytable)
            }, error=function(e)
                {
                    print("Please check the output files.")
                    Filteredsummarytable <- data.frame(matrix(ncol = 6, nrow = 0))
                    return (Filteredsummarytable)
                }
                )
            }
            Filteredsummarytable <- getData()
   } else {
       print ("No results from GEO Datasets.")
       Filteredsummarytable <- data.frame(matrix(ncol = 6, nrow = 0))
       return (Filteredsummarytable)
    }
} #require get_meshterms()
