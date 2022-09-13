# DAGItoolsr

##Disease-associated gene identification and analysis tools

Implements methods to identify and analyse disease-associated genes. Several options are available: (1) gene search from databases, and (2) gene search from differential expression data. It also provides options to predict protein-protein interaction (via STRINGdb) and perform functional enrichment (via clusterProfiler).

To install: 
> install.packages("remotes")<br>
> install_github("NurSyatila/DAGItoolsr")

To use DAGItoolsr package:
> library(DAGItoolsr)

To run a gene search workflow for a disease, e.g. endometriosis:
> queryTerms <- c("endometriosis", "endometrioma") # include disease name and keywords for search
> results <- genesearch_workflow(queryTerms)

A directory named 'GeneSearch' will be created which contains results from database gene search, PPI and/or functional enrichment analyses for a list of genes matching the query disease name / keywords.

