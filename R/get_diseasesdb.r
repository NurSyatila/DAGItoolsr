#' Gene search in DISEASES database
#'
#' Get disease-associated genes (Entrez Gene identifiers) from DISEASES database
#' @param queryTerms A vector of disease name and/or keywords
#' @return A list of protein-coding genes
#' @return A CSV-formmated file containing results from DISEASES database
#' @return A text-formatted file containing a list of genes (DISEASES_Genes.txt)
#' @return A text-formatted file containing a list of protein-coding genes (DISEASES_PCGenes.txt)
#' @examples
#' queryTerms <- c("endometriosis");
#' results <- get_diseasesdb(queryTerms);
#' @export
#' @importFrom stats na.omit
#' @importFrom utils download.file head read.csv write.csv
get_diseasesdb <- function(queryTerms) {
    gene_file <- "DISEASES_Genes.txt"
    pcgene_file <- "DISEASES_PCGenes.txt"
    options(download.file.method="curl",timeout = 300)
    #options(download.file.method="curl")
    patterns <- paste(queryTerms, collapse = "|")
    print ("Processing DISEASES (Text Mining category)...")
    url_text_mining <- 'https://download.jensenlab.org/human_disease_textmining_filtered.tsv'
    tryCatch(suppressMessages(download.file(url_text_mining, destfile='DISEASES_human_disease_textmining_filtered.tsv', quiet = TRUE)), error = function(e) print('DISEASE textmining file: did not work out'))
    #suppressMessages(download.file(url_text_mining, destfile='DISEASES_human_disease_textmining_filtered.tsv'))
    if (file.exists('DISEASES_human_disease_textmining_filtered.tsv') && file.size('DISEASES_human_disease_textmining_filtered.tsv') > 0) {
        text_table <- read.csv(file = 'DISEASES_human_disease_textmining_filtered.tsv', sep = '\t',header=FALSE)
        text_results1 <- text_table[grepl(patterns,tolower(text_table$V4)), ]
        text_results2 <- text_results1[as.numeric(text_results1$V6)>3.0, ]
        if(nrow(text_results2) > 0){
            text_results2['Type']='Text mining'
            colnames(text_results2) <- c('EnsemblID','GeneSymbol','DO','Disease','ZScore','ConfidentScore','URL','Type')
            summarytable1 <- text_results2[,c(ncol(text_results2),1:(ncol(text_results2)-1))]
            suppressMessages(write.csv(summarytable1, file="DISEASESTextSummary.csv",row.names = FALSE))
            text_gene_list <- na.omit(summarytable1$GeneSymbol)
        } else { text_gene_list <- list() }
    } else {
        print ("Error: Unable to download file (DISEASES_human_disease_textmining_filtered.tsv).")
        text_gene_list <- list()
    }
    print ("Processing DISEASES (Knowledge category)...")
    url_knowledge <- 'https://download.jensenlab.org/human_disease_knowledge_filtered.tsv'
    tryCatch(suppressMessages(download.file(url_knowledge, destfile='DISEASES_human_disease_knowledge_filtered.tsv', quiet = TRUE)),
        error = function(e) print('DISEASE knowledge file: did not work out'))
    #suppressMessages(download.file(url_knowledge, destfile='DISEASES_human_disease_knowledge_filtered.tsv'))
    if (file.exists('DISEASES_human_disease_knowledge_filtered.tsv') && file.size('DISEASES_human_disease_knowledge_filtered.tsv') > 0) {
        knowledge_table <- read.csv(file = 'DISEASES_human_disease_knowledge_filtered.tsv', sep = '\t',header=FALSE)
        knowledge_results1 <- knowledge_table[grepl(patterns,tolower(knowledge_table$V4)), ]
        knowledge_results2 <- knowledge_results1[as.numeric(knowledge_results1$V7)>3, ]
        if(nrow(knowledge_results2) > 0){
            knowledge_results2['Type']='Knowledge'
            colnames(knowledge_results2) <- c('EnsemblID','GeneSymbol','DO','Disease','SourceDB','EvidenceType','ConfidenceScore','Type')
            summarytable2 <- knowledge_results2[,c(ncol(knowledge_results2),1:(ncol(knowledge_results2)-1))]
            suppressMessages(write.csv(summarytable2, file="DISEASESKnowledgeSummary.csv",,row.names = FALSE))
            knowledge_gene_list <- na.omit(summarytable2$GeneSymbol)
        } else { knowledge_gene_list <- list() }
    } else {
        print ("Error: Unable to download file (DISEASES_human_disease_knowledge_filtered.tsv).")
        knowledge_gene_list <- list()
    }
    print ("Processing DISEASES (Experiments category)...")
    url_experiments <- 'https://download.jensenlab.org/human_disease_experiments_filtered.tsv'
    tryCatch(suppressMessages(download.file(url_experiments, destfile='DISEASES_human_disease_experiments_filtered.tsv', quiet = TRUE)),
        error = function(e) print('DISEASE experiment file: did not work out'))
    #suppressMessages(download.file(url_experiments, destfile='DISEASES_human_disease_experiments_filtered.tsv'))
    if (file.exists('DISEASES_human_disease_experiments_filtered.tsv') && file.size('DISEASES_human_disease_experiments_filtered.tsv') > 0) {
        exps_table <- read.csv(file = 'DISEASES_human_disease_experiments_filtered.tsv', sep = '\t',header=FALSE)
        exps_results1 <- exps_table[grepl(patterns,tolower(exps_table$V4)), ]
        exps_results2 <- exps_results1[as.numeric(exps_results1$V7)>2, ]
        if(nrow(exps_results2) > 0){
            exps_results2['Type']='Experiments'
            colnames(exps_results2) <- c('EnsemblID','GeneSymbol','DO','Disease','SourceDB','RankScore','ConfidenceScore','Type')
            summarytable3 <- exps_results2[,c(ncol(exps_results2),1:(ncol(exps_results2)-1))]
            suppressMessages(write.csv(summarytable3, file="DISEASESExpSummary.csv",,row.names = FALSE))
            exps_gene_list <- na.omit(summarytable3$GeneSymbol)
        } else { exps_gene_list <- list() }
    } else {
        print ("Error: Unable to download file (DISEASES_human_disease_experiments_filtered.tsv).")
        exps_gene_list <- list()
    }
    gene_diseases_list <- unique(unlist(c(text_gene_list,knowledge_gene_list,exps_gene_list),recursive=TRUE))
    if (length(gene_diseases_list)>0){
        PCGenesx <- get_PCGenes(gene_diseases_list, gene_file,pcgene_file,"SYMBOL")
    } else {
      print ("No results from DISEASES")
      PCGenesx <- list() }
    return (PCGenesx)
}