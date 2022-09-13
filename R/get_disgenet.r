#' Gene search in DisGeNET
#'
#' Get disease-associated genes (Entrez Gene identifiers) from DisGeNET (via DisGeNET2r)
#' @param queryTerms A vector of disease name and/or keywords
#' @param userEmail Registered email to DisGeNET account
#' @param userPassword Password to DisGeNET account
#' @return A list of protein-coding genes
#' @return A CSV-formmated file containing results from DisGeNET
#' @return A text-formatted file containing a list of genes (DisGeNET_Genes.txt)
#' @return A text-formatted file containing a list of protein-coding genes (DisGeNET_PCGenes.txt)
#' @examples
#' queryTerms <- c("endometriosis");
#' userEmail <- "you@gmail.com";
#' userPassword <- "yourpassword";
#' results <- get_disgenet(queryTerms,userEmail,userPassword);
#' @export
#' @importFrom utils write.csv
get_disgenet <- function(queryTerms,userEmail,userPassword) {
    gene_file <- "DisGeNET_Genes.txt"
    pcgene_file <- "DisGeNET_PCGenes.txt"
    MeSHTerms <- get_meshterms(queryTerms,"mesh_id")
    print (MeSHTerms)
    if (length(MeSHTerms)>0){
        get_disgenetx <- function(){
            tryCatch({
                print ('Processing (DisGeNET)..')
                genelist <- list()
                disgenet_api_key <- disgenet2r::get_disgenet_api_key(email = userEmail,password = userPassword )
                Sys.setenv(DISGENET_API_KEY= disgenet_api_key)
                results <- disgenet2r::disease2gene( disease = MeSHTerms, vocabulary = "MESH", database = "CURATED", score = c( 0.4,1 ))
                summarytable <- disgenet2r::extract(results)
                write.csv(summarytable, file="DisGeNETSummary.csv",row.names = FALSE)
                dbList <- as.character(unique(summarytable$geneid))
                if (length(dbList)>0){
                   MappedGenesx <- get_PCGenes(dbList, gene_file,pcgene_file,"ENTREZID")
                } else {
                  print ("No results from DisGeNET.")
                  MappedGenesx <- list() }
            }, error=function(e)
              {
                print ("Error: Unable to get data from DisGeNET.")
                MappedGenesx <- list()
              }
            )
        }
        MappedGenesx <- get_disgenetx()
    } else { MappedGenesx <- list() }
    return(MappedGenesx)
}