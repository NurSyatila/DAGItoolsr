#' Gene search in PheGenI
#'
#' Get disease-associated genes (Entrez Gene identifiers) from PheGenI
#' @param queryTerms A vector of disease name and/or keywords
#' @return A list of protein-coding genes
#' @return A CSV-formmated file containing results from PheGenI
#' @return A text-formatted file containing a list of genes (PheGenI_Genes.txt)
#' @return A text-formatted file containing a list of protein-coding genes (PheGenI_PCGenes.txt)
#' @examples
#' queryTerms <- c("endometriosis");
#' results <- get_phegeni(queryTerms);
#' @export
#' @importFrom utils download.file write.csv
get_phegeni <- function(queryTerms) {
    gene_file <- "PheGenI_Genes.txt"
    pcgene_file <- "PheGenI_PCGenes.txt"
    options(download.file.method="curl",timeout = 300)
    url_download <- 'https://www.ncbi.nlm.nih.gov/projects/gap/eqtl/EpiViewBE.cgi?type=dl.tab'
    tryCatch(suppressMessages(download.file(url_download, destfile='phegeni.tab', quiet = TRUE)),
        error = function(e) print('phegeni.tab: did not work out'))
    #suppressMessages(download.file(url_download, destfile='phegeni.tab'))
    if (file.exists('phegeni.tab') && file.size('phegeni.tab') > 0) {
        print ("Processing (PheGenI)... ")
        phegeni <- readLines(file('phegeni.tab'))
        phegeni_tab <- strsplit(phegeni,'\t')
        phegeni_tab_short <- lapply(phegeni_tab, function(x) {x[c(1:13)]})
        patterns <- tolower(paste(queryTerms, collapse = "|"))
        summarylist <- Filter(function(x) grepl(patterns, tolower(x[2])), phegeni_tab_short)
        summarytable <- do.call(rbind,summarylist)
        if (length(summarytable[,1])>0){
            colnames(summarytable) <- c("ID","Trait", "SNP rs", "Context", "Gene", "Gene ID", "Gene 2", "Gene ID 2", "Chromosome", "Location", "P-Value", "Source", "PubMed")
            write.csv(summarytable, file="PheGenISummary.csv",row.names = FALSE)
            genelist <- unique(unlist(c(summarytable[,6],summarytable[,8]),recursive=TRUE))
            if (length(genelist)>0){
              MappedGenesx <- get_PCGenes(genelist, gene_file,pcgene_file,"ENTREZID")
            } else {
              print ("No results from PheGenI.")
              MappedGenesx <- list() }
        } else {
          print ("No results from PheGenI.")
          MappedGenesx <- list() }
        close(file('phegeni.tab'))
    } else {
        print ("Error: Unable to download file (ncbi/gap/phegeni).")
        MappedGenesx <- list()
    }
    return(MappedGenesx)
}
