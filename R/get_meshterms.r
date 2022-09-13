#' Get MeSH terms or identifiers
#'
#' Get MeSH terms or identifiers associated to a disease (via rentrez)
#' @param queryTerms A vector of disease name and/or keywords
#' @param meshtype type of MeSH output, "mesh_id" or "mesh_term"
#' @return A list of MeSH terms or identifers
#' @examples
#' queryTerms <- c("endometriosis");
#' results <- get_meshterms(queryTerms,"mesh_id");
#' results <- get_meshterms(queryTerms,"mesh_term");
#' @export
get_meshterms <- function(queryTerms,meshtype) {
    new_queryTerms <- character()
    for (var in queryTerms){
        query_term <- paste('"',var,'"[MESH]',sep='')
        new_queryTerms <- c(new_queryTerms, query_term)
    }
    new_queryTerms <- paste(new_queryTerms, collapse=' OR ')
    qTerms <- paste('(("main heading"[TYPE]) AND (',new_queryTerms,')',sep='')
    dbSearch <- suppressMessages(rentrez::entrez_search(db='mesh',term=qTerms,use_history=TRUE,retmax=10000))
    if (length(dbSearch$ids)>0){
       getData <- function(){
           tryCatch({
               print ("Processing (MESH)...")
               MESHitems <- list()
               if (length(dbSearch$ids)>250){
                    chunks <- split(dbSearch$ids, ceiling(seq_along(dbSearch$ids)/250))
                    #loop over
                    for (i in 1:length(chunks)){
                      summary <- suppressMessages(rentrez::entrez_summary(db="mesh",id=unlist(chunks[1],recursive=TRUE)))
                      if (meshtype=="mesh_id"){ MESHitems <- c(MESHitems, gsub('^68','D',summary$uid))}
                      if (meshtype=="mesh_term"){ MESHitems <- c(MESHitems, summary$ds_meshterms)}
                    }
               } else {
                    summary <- suppressMessages(rentrez::entrez_summary(db="mesh",id=unlist(dbSearch$ids,recursive=TRUE)))
                    if (meshtype=="mesh_id"){ MESHitems <- c(MESHitems, gsub('^68','D',summary$uid))}
                    if (meshtype=="mesh_term"){ MESHitems <- c(MESHitems, summary$ds_meshterms)}
                }
                MESHitems <- unique(unlist(MESHitems,recursive=TRUE))
            }, error=function(e) {
             print ("Error: Please check the output files.")
             MESHitems <- list() }) }
            MESHitems <- getData()
   } else {
       print ("MESH terms cannot be identified.")
       MESHitems <- list()
    }
    return (MESHitems)
}