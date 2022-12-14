% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_disgenet.r
\name{get_disgenet}
\alias{get_disgenet}
\title{Gene search in DisGeNET}
\usage{
get_disgenet(queryTerms, userEmail, userPassword)
}
\arguments{
\item{queryTerms}{A vector of disease name and/or keywords}

\item{userEmail}{Registered email to DisGeNET account}

\item{userPassword}{Password to DisGeNET account}
}
\value{
A list of protein-coding genes

A CSV-formmated file containing results from DisGeNET

A text-formatted file containing a list of genes (DisGeNET_Genes.txt)

A text-formatted file containing a list of protein-coding genes (DisGeNET_PCGenes.txt)
}
\description{
Get disease-associated genes (Entrez Gene identifiers) from DisGeNET (via DisGeNET2r)
}
\examples{
queryTerms <- c("endometriosis");
userEmail <- "you@gmail.com";
userPassword <- "yourpassword";
results <- get_disgenet(queryTerms,userEmail,userPassword);
}
