#' Assign sample groupings for a GEO Series)
#'
#' Assign sample groupings for a GEO Series through automatic classification (required by analyse_deg function)
#' @param queryTerms A vector of disease name and/or keywords
#' @param GSEaccession GSE Accession
#' @param GSEplatform GSE Platform
#' @return A vector of GEOquery results and sample groupings e.g. c(0,0,0,0,1,1,1,1)
#' @examples
#' queryTerms <- c("endometriosis");
#' GSEaccession <- "GSE23339";
#' GSEplatform <- "GPL6102";
#' results <- get_samplegroupings(queryTerms,GSEaccession,GSEplatform);
#' @export
#' @importFrom Biobase pData
get_samplegroupings <- function(queryTerms,GSEaccession,GSEplatform) {
    parse_gsefile <- function(GSEaccession){
        tryCatch({
            gset <- suppressMessages(GEOquery::getGEO(GSEaccession, GSEMatrix =TRUE, AnnotGPL=TRUE,destdir="."))
            if (length(gset) > 1) idx <- grep(GSEplatform, attr(gset, "names")) else idx <- 1
            gset <- gset[[idx]]
            if (!is.null(pData(gset)$characteristics_ch1)==TRUE){samples <- pData(gset)$characteristics_ch1}
            else {samples <- pData(gset)$description}
            if (!is.null(samples)==TRUE){
              groups <- unique(samples)
              if (length(groups)==1){
                print(paste("group A: ",groups,collapse=""))
                print ("** Error: Only one sample types detected. The GSE data set will not be used. **")
                sampleGrouping <- ""
                return (list(gset,sampleGrouping))}
              if (length(groups)==2){
                # divide groups into two
                sample_groupings <- gsub(groups[1],"0",samples)
                sampleGrouping <- gsub(groups[2],"1",sample_groupings)
                print(paste("group A: ",groups[1],collapse=""))
                print(paste("group B: ",groups[2],collapse=""))
                return (list(gset,sampleGrouping))
              }
              if (length(groups)>2){
                # divide groups based on presence of healthy condition
                conds_normal <- c("healthy","normal","control")
                patternsh <- tolower(paste(conds_normal, collapse = "|"))
                new_groups <- grepl(patternsh, tolower(groups))
                group1 <- group2 <- character()
                for(i in seq(1:length(new_groups))){
                  if (new_groups[i]==TRUE){group1 <- c(group1,groups[i])}
                  if (new_groups[i]==FALSE){group2 <- c(group2,groups[i])}
                }
                if (length(group1) !=0 && length(group2) !=0){
                  sample_groupings <- gsub(paste(group1,collapse="|"),"0",samples)
                  sampleGrouping <- gsub(paste(group2,collapse="|"),"1",sample_groupings)
                  print(paste("group A: ",paste(group1,collapse=" | "),collapse=""))
                  print(paste("group B: ",paste(group2,collapse=" | "),collapse=""))
                  return (list(gset,sampleGrouping))
                } else {
                  # divide groups based on presence of query term (disease name)
                  patterns <- tolower(paste(queryTerms, collapse = "|"))
                  new_groups <- grepl(patterns, tolower(groups))
                  group1 <- group2 <- character()
                  for(i in seq(1:length(new_groups))){
                    if (new_groups[i]==TRUE){group1 <- c(group1,groups[i])}
                    if (new_groups[i]==FALSE){group2 <- c(group2,groups[i])}
                  }
                  if (length(group1) !=0 && length(group2) !=0){
                    sample_groupings <- gsub(paste(group1,collapse="|"),"0",samples)
                    sampleGrouping <- gsub(paste(group2,collapse="|"),"1",sample_groupings)
                    print(paste("group A: ",paste(group1,collapse=" | "),collapse=""))
                    print(paste("group B: ",paste(group2,collapse=" | "),collapse=""))
                    return (list(gset,sampleGrouping))
                  } else {
                    # divide groups based on presence of disease term
                    new_groups <- grepl("disease", tolower(groups))
                    group1 <- group2 <- character()
                    for(i in seq(1:length(new_groups))){
                      if (new_groups[i]==TRUE){group1 <- c(group1,groups[i])}
                      if (new_groups[i]==FALSE){group2 <- c(group2,groups[i])}
                    }
                    if (length(group1) !=0 && length(group2) !=0){
                      sample_groupings <- gsub(paste(group1,collapse="|"),"0",samples)
                      sampleGrouping <- gsub(paste(group2,collapse="|"),"1",sample_groupings)
                      print(paste("group A: ",paste(group1,collapse=" | "),collapse=""))
                      print(paste("group B: ",paste(group2,collapse=" | "),collapse=""))
                      return (list(gset,sampleGrouping))
                    } else {
                      #divide groups based on distinct words
                      ul <- unique(unlist(strsplit(samples, " ")))
                      keyterms <- list()
                      for (i in ul){
                        if (all(grepl(i, samples))=="FALSE"){
                          groupk1 <- groupk2 <- character()
                          character()
                          for (seq in samples){
                            if (grepl(i,seq)==TRUE){groupk1 <- c(groupk1,seq)}
                            if (grepl(i,seq)==FALSE){groupk2 <- c(groupk2,seq)}
                          }
                          df <- data.frame(i,paste(groupk1,collapse='||'),paste(groupk2,collapse='||'),abs(length(groupk1)-length(groupk2)))
                          keyterms[[i]] <- df
                        }
                      }
                      summarytable <- do.call(rbind,keyterms)
                      group1 <- unlist(strsplit(summarytable[order(summarytable[,4]),][,2][1],"[||]"),recursive=TRUE)
                      group1 <- group1[group1 != ""]
                      group2 <- unlist(strsplit(summarytable[order(summarytable[,4]),][,3][1],"[||]"),recursive=TRUE)
                      group2 <- group2[group2 != ""]
                      if (length(group1) !=0 && length(group2) !=0){
                        sample_groupings <- gsub(paste(group1,collapse="|"),"0",samples)
                        sampleGrouping <- gsub(paste(group2,collapse="|"),"1",sample_groupings)
                        print(paste("group A: ",paste(group1,collapse=" | "),collapse=""))
                        print(paste("group B: ",paste(group2,collapse=" | "),collapse=""))
                        return (list(gset,sampleGrouping))
                      } else {
                        print ("** Error: More/less than two sample types detected. The GSE data set will not be used. **")
                        sampleGrouping <- ""
                        return (list(gset,sampleGrouping))
                      }
                    }
                  }}
              }
            } else {
              print("Error (2): Sample groupings cannot be determined.")
              sampleGrouping <- ""
              return (list(gset,sampleGrouping))
            }
        }, error=function(e)
                {
                    print("Error (1): Unable to download / parse the GSE file.")
                    gset <- ""
                    sampleGrouping <- ""
                    return (list(gset,sampleGrouping))
                }
        )
    }
    deg_list <- parse_gsefile(GSEaccession)
    return (deg_list)
}
