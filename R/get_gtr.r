#' Gene search in GTR
#'
#' Get disease-associated genes (Entrez Gene identifiers) from GTR (via rentrez)
#' @param queryTerms A vector of disease name and/or keywords
#' @return A list of protein-coding genes
#' @return A CSV-formmated file containing results from GTR database
#' @return A text-formatted file containing a list of genes (GTR_Genes.txt)
#' @return A text-formatted file containing a list of protein-coding genes (GTR_PCGenes.txt)
#' @examples
#' queryTerms <- c("endometriosis");
#' results <- get_gtr(queryTerms);
#' @export
#' @importFrom utils write.csv
get_gtr <- function(queryTerms) {
   gene_file <- "GTR_Genes.txt"
   pcgene_file <- "GTR_PCGenes.txt"
   new_query_terms <- character()
   for (var in queryTerms){
       query_term <- paste('"',var,'"[DISNAME]',sep='')
       new_query_terms <- c(new_query_terms, query_term)
   }
   qTerms <- paste('(',paste(new_query_terms,collapse=' OR '),')',sep="")
    dbSearch <- suppressMessages(rentrez::entrez_search(db='gtr',term=qTerms,use_history=TRUE,retmax=10000))
   if (length(dbSearch$ids)>0){
       getData <- function(){
           tryCatch({
               print ("Processing (GTR)...")
               cuilist <- list()
               if (length(dbSearch$ids)>250){
                    chunks <- split(dbSearch$ids, ceiling(seq_along(dbSearch$ids)/250))
                    #loop over
                    for (i in 1:length(chunks)){
                      summary <- suppressMessages(rentrez::entrez_summary(db="gtr",id=unlist(chunks[1],recursive=TRUE)))
                      esummaries <- suppressMessages(rentrez::extract_from_esummary(summary, "conditionlist"))
                      for (i in 1:length(esummaries[1,])){
                          dfcond <- data.frame(cbind(esummaries[,i]$name, esummaries[,i]$cui))
                          dfcondcui <- dplyr::filter(dfcond, grepl(paste(queryTerms,collapse="|"),tolower(dfcond[,1])))
                          cuilist <- c(cuilist,dfcondcui$X2)
                      }
                    }
               } else {
                    summary <- suppressMessages(rentrez::entrez_summary(db="gtr",id=unlist(dbSearch$ids,recursive=TRUE)))
                    esummaries <- suppressMessages(rentrez::extract_from_esummary(summary, "conditionlist"))
                    if (length(unlist(dbSearch$ids,recursive=TRUE))>1){
                      for (i in 1:length(esummaries[1,])){
                          dfcond <- data.frame(cbind(esummaries[,i]$name, esummaries[,i]$cui))
                          dfcondcui <- dplyr::filter(dfcond, grepl(paste(queryTerms,collapse="|"),tolower(dfcond[,1])))
                          cuilist <- c(cuilist,dfcondcui$X2)
                      }
                    } else {
                        dfcond <- data.frame(cbind(esummaries$name, esummaries$cui))
                        dfcondcui <- dplyr::filter(dfcond, grepl(paste(queryTerms,collapse="|"),tolower(dfcond[,1])))
                        cuilist <- c(cuilist,dfcondcui[,2])
                    }
               }
               cuilist <- unique(unlist(cuilist,recursive=TRUE))
               if (length(cuilist)>0){
                  genelist <- list()
                  summarylist <- list()
                  for (cuid in cuilist){
                    url <- gsub(' ','',paste('https://www.ncbi.nlm.nih.gov/gtr/conditions/',cuid,collapse=''))
                    page_html <- gsub('\n','',toString(rvest::read_html(url)))
                    gtrmatches <- stringr::str_match(page_html,'<ul class="associate_genes gtr-reset-list">.*</ul>')[[1]]
                    gtrmatches2 <- strsplit(strsplit(gtrmatches,"</ul>")[[1]][1], '<li>')[[1]][-1]
                    for (g in gtrmatches2){
                      geneid <- gsub('/gtr/genes/','',stringr::str_match_all(g,'/gtr/genes/\\d+')[1])
                      genesymbol <- gsub("<.*?>", "", strsplit(g,"</a>")[[1]][1])
                      gtemp <- gsub("<.*?>", "", strsplit(g,"</a>")[[1]][3])
                      genename <- strsplit(gtemp,"Summary: ")[[1]][2]
                      synonym <- gsub("Also known as: ","",strsplit(gtemp,"Summary: ")[[1]][1])
                      namex <- paste(cuid,"_",geneid,sep="")
                      df <- data.frame(cuid, geneid, genesymbol, genename, synonym)
                      summarylist[[namex]] <- df
                      genelist <- c(genelist,geneid)
                    }
                  }
                  summarytable <- do.call(rbind,summarylist)
                  colnames(summarytable) <- c("CUID","GeneID","GeneSymbol","GeneName","Synonym")
                  write.csv(summarytable, file="GTRSummary.csv",row.names = FALSE)
                  genelist <- unique(unlist(genelist,recursive=TRUE))
                  if (length(genelist)>0){
                    MappedGenesx <- get_PCGenes(genelist, gene_file,pcgene_file,"ENTREZID")
                  } else {
                    print ("No results from GTR.")
                    MappedGenesx <- list() }
               } else {
                 print ("No results from GTR.")
                 MappedGenesx <- list()
               }
            }, error=function(e)
                {
                    print("Error: Please check the output files.")
                    MappedGenesx <- list()
                }) }
       MappedGenesx <- getData()
   } else {
       print ("No results from GTR.")
       MappedGenesx <- list()
    }
   return(MappedGenesx)
}