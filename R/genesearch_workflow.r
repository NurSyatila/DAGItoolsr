#' Database Gene Search Worflow
#'
#' Get disease-associated genes (Entrez Gene identifiers) from available public resources
#' @param queryTerms A vector of disease name and/or keywords
#' @param disgenet An option to include results from DisGeNET. Note that This option will require email and password to DisGeNET account
#' @param email Registered email to DisGeNET account
#' @param password Password to DisGeNET account
#' @return A directory named GeneSearch
#' @return Gene Search: CSV-formmated files containing results from individual database search
#' @return Gene Search: Text-formatted files containing a list of genes from individual database search ([Database]_Genes.txt)
#' @return Gene Search: Text-formatted files containing a list of protein-coding genes from individual database search ([Database]_PCGenes.txt)
#' @return Gene Search: A text-formatted file containing a list of protein-coding genes from all database searches (DatabaseGenes.txt)
#' @return PPI: A summary of PPI analysis (PPISummary.txt)
#' @return PPI: A list of filtered genes from PPI where isolated nodes are discarded (FilteredGenes.txt)
#' @return PPI: Top genes from PPI with highest interaction (Top_PPI_genes.txt)
#' @return PPI: A list of gene clusters with details on enriched KEGG Pathway (GeneCluster.txt)
#' @return PPI: A PDF file containing visualization plots from PPI analysis (STRING_PPI_Network.pdf)
#' @return Functional Enrichment: A separate directory of functional enrichment results
#' @examples
#' queryTerms <- c("endometriosis");
#' results <- genesearch_workflow(queryTerms,disgenet=FALSE);
#' results <- genesearch_workflow(queryTerms,disgenet=TRUE,"you@gmail.com","yourpassword");
#' @export
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics barplot text
genesearch_workflow <- function(queryTerms,disgenet=FALSE,email=NULL,password=NULL){
  if(disgenet==FALSE){
    print ("DisGeNET will not be included.")
    dir.create("GeneSearch")
    setwd("GeneSearch")
    clinvar <- get_clinvar(queryTerms)
    medgen <- get_medgen(queryTerms)
    omim <- get_omim(queryTerms)
    gtr <- get_gtr(queryTerms)
    phegeni <- get_phegeni(queryTerms)
    gwascatalog <- get_gwascatalog(queryTerms)
    diseasesdb <- get_diseasesdb(queryTerms)
    allgenes <- c(clinvar, medgen, omim, gtr, phegeni, gwascatalog, diseasesdb)
    allgenes <- unique(unlist(allgenes,recursive=TRUE))
    if (length(allgenes)>0){
      write(allgenes,file="DatabaseGenes.txt")
      # get the distribution of PCG from all databases
      pdf('GeneSearch.pdf',paper="a4r")
      colors <- c("darkblue","red","green", "orange","cyan","magenta","black")
      dtb <- c("ClinVar", "MedGen","OMIM", "GTR", "PheGenI", "GWAS", "DISEASES")
      nogenes <- c(length(clinvar),length(medgen),length(omim),length(gtr),length(phegeni),length(gwascatalog),length(diseasesdb))
      barp <- barplot(nogenes, main = "Protein-coding genes from database search", ylab = "No. of genes", border = "black", col = colors, names=dtb, las=3 )
      text(barp, nogenes + 6.0, labels = nogenes)
      dev.off()
      if (length(allgenes)<400){
        dir.create("PPI")
        setwd("PPI")
        analyse_ppi_network(allgenes)
        setwd("../")
      } else {
        print ("Number of genes > 400. PPI analysis cannot be performed. Please use analyse_ppi_network(). Running enrichment analysis instead.")
        GSEnames <- c("BP","MF","CC","KEGG","DO")
        for (GSEname in GSEnames){
            print (paste("**", GSEname, "**", sep=" "))
            suppressMessages(get_enrichment_results(allgenes,GSEname,"ORA"))
        }
      }
    } else {print ("No results found.")}
    setwd("../")
  }
  if (disgenet==TRUE){
    if (!is.null(email) && !is.null(password)){
      print ("DisGeNET data will be included.")
      dir.create("GeneSearch")
      setwd("GeneSearch")
      clinvar <- get_clinvar(queryTerms)
      medgen <- get_medgen(queryTerms)
      omim <- get_omim(queryTerms)
      gtr <- get_gtr(queryTerms)
      phegeni <- get_phegeni(queryTerms)
      gwascatalog <- get_gwascatalog(queryTerms)
      # note that disgenet requires email and password to the registered account in DisgeNET
      #email <- 'nsag@ukm.edu.my'
      #password <- 'ns65ns65'
      diseasesdb <- get_diseasesdb(queryTerms)
      disgenetdb <- get_disgenet(queryTerms, email, password)
      allgenes <- c(clinvar, medgen, omim, gtr, phegeni, gwascatalog, disgenetdb, diseasesdb)
      allgenes <- unique(unlist(allgenes,recursive=TRUE))
      if (length(allgenes)>0){
        write(allgenes,file="DatabaseGenes.txt")
        # get the distribution of PCG from all databases
        pdf('GeneSearch.pdf',paper="a4r")
        colors <- c("darkblue","red","green", "orange","cyan","magenta","black","yellow")
        dtb <- c("ClinVar", "MedGen","OMIM", "GTR", "PheGenI", "GWAS", "DisGeNET", "DISEASES")
        nogenes <- c(length(clinvar),length(medgen),length(omim),length(gtr),length(phegeni),length(gwascatalog),length(disgenetdb),length(diseasesdb))
        barp <- barplot(nogenes, main = "Protein-coding genes from database search", ylab = "No. of genes", border = "black", col = colors, names=dtb, las=3 )
        text(barp, nogenes + 6.0, labels = nogenes)
        dev.off()
        if (length(allgenes)<400){
          dir.create("PPI")
          setwd("PPI")
          analyse_ppi_network(allgenes)
          setwd("../")
        } else {
          print ("Number of genes > 400. PPI analysis cannot be performed. Please use analyse_ppi_network(). Running enrichment analysis instead.")
          GSEnames <- c("BP","MF","CC","KEGG","DO")
          for (GSEname in GSEnames){
              print (paste("**", GSEname, "**", sep=" "))
              suppressMessages(get_enrichment_results(allgenes,GSEname,"ORA"))
          }
        }
      } else {print ("No results found.")}
      setwd("../")
    }
    else {print ("Please provide email and password to DisGeNET account.")}
  }
}
