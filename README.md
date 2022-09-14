# DAGItoolsr

##Disease-associated gene identification and analysis tools

Implements methods to identify and analyse disease-associated genes. Several options are available: (1) gene search from databases, and (2) gene search from differential expression data. It also provides options to predict protein-protein interaction (via STRINGdb) and perform functional enrichment (via clusterProfiler).

To install: 
> install.packages("remotes")<br>
> install_github("NurSyatila/DAGItoolsr")

To use DAGItoolsr package:
> library(DAGItoolsr)

## Gene Search Workflow

To run a gene search workflow for a disease, e.g. endometriosis:
> queryTerms <- c("endometriosis", "endometrioma") # include disease name and keywords for search<br>
> results <- genesearch_workflow(queryTerms)

A directory named 'GeneSearch' will be created which contains results from database gene search, PPI and/or functional enrichment analyses for a list of genes matching the query disease name / keywords.

## Differential Expression Workflow (Single GEO Series)

To run a differential expression workflow for a disease, e.g. endometriosis:
> queryTerms <- c("endometriosis", "endometrioma") # include disease name and keywords for search<br>
> GSEaccession <- "GSE23339"
> GSEplatform <- "GPL6102"
> results <- differentialexpression_single_workflow(queryTerms,GSEaccession,GSEplatform)

A directory named 'GSE23339' will be created which contains results from gene expresion analysis, PPI and/or functional enrichment analyses for a list of differentially expressed protein-coding genes. The queryTerms parameter is used to help assigning samples into groups (normal vs. diseased state).

## Differential Expression Workflow

To run a gene search workflow for a disease, e.g. endometriosis:
> queryTerms <- c("endometriosis", "endometrioma") # include disease name and keywords for search<br>
> results <- differentialexpression_workflow(queryTerms)

A directory named 'DifferentialExpression' will be created which contains results from GEO series search, gene expresion analysis, PPI and/or functional enrichment analyses for individual GEO series matching the query disease name / keywords.
