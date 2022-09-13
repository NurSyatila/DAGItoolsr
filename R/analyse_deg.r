#' Analyse Gene Expression Data (GEO Series)
#'
#' Perform differential expression analysis and retrieve differentially expressed protein-coding genes from GEO Series
#' @param queryTerms A vector of disease name and/or keywords
#' @param GSEaccession GSE Accession
#' @param GSEplatform GSE Platform
#' @return A dataframe of differentially expressed protein-coding genes
#' @return A CSV-formmated file containing results of differentially expressed protein-coding genes (FilteredDE.csv)
#' @return A text file containing Entrez identifiers of differentially expressed protein-coding genes (DEGenes.txt)
#' @return A PDF-formatted file containing visualization plots (Plots.pdf)
#' @return A text file containing top 20 genes in tabulated form (Top20DEGenes.txt)
#' @examples
#' queryTerms <- c("endometriosis");
#' GSEaccession <- "GSE23339";
#' GSEplatform <- "GPL6102";
#' results <- analyse_deg(queryTerms,GSEaccession,GSEplatform);
#' @export
#' @importFrom Biobase fData
#' @importFrom grDevices dev.off palette pdf
#' @importFrom graphics abline boxplot hist legend par
#' @importFrom stats filter model.matrix na.omit quantile
#' @importFrom utils write.csv write.table
#' @importFrom Biobase exprs annotation
analyse_deg <- function(queryTerms,GSEaccession,GSEplatform){
    output1 <- get_samplegroupings(queryTerms,GSEaccession,GSEplatform)
    gset <- output1[[1]]
    sampleGrouping <- output1[[2]]
    if (length(sampleGrouping) > 1){
        gene_annot <- suppressMessages(fData(gset))
        gene_ids <- unlist(strsplit(gene_annot$"Gene ID","///"),recursive=TRUE)
        pc_genes <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = gene_ids, c("ENTREZID", "GENETYPE"), keytype = "ENTREZID"))
        pc_genes <- dplyr::filter(pc_genes, grepl("protein",tolower(pc_genes$GENETYPE)))
        pc_genes <- unique(pc_genes$ENTREZID)
        gene_annot2 <- gene_annot[gene_annot$"Gene ID" %in% pc_genes,]
        gset2 <- gset[rownames(gset) %in% gene_annot2$"ID",]
        # log2 transformation
        ex <- exprs(gset2)
        qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
        LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
        if (LogC ) { 
            ex[which(ex <= 0)] <- NaN
            ex <- log2(ex) 
        }
        # assign samples to groups and set up design matrix
        gs <- factor(sampleGrouping)
        groups <- make.names(c("groupA","groupB"))
        levels(gs) <- groups
        gset2$group <- gs
        design <- model.matrix(~group + 0, gset2)
        colnames(design) <- levels(gs)
        fit <- limma::lmFit(gset2, design)  # fit linear model
        # set up contrasts of interest and recalculate model coefficients
        cts <- paste(groups[1], groups[2], sep="-")
        cont.matrix <- limma::makeContrasts(contrasts=cts, levels=design)
        fit2 <- limma::contrasts.fit(fit, cont.matrix)
        fit2 <- limma::eBayes(fit2, 0.01)
        full_results <- limma::topTable(fit2, adjust="BH", sort.by="B", number=Inf)
        # save DEG results
        new_results <- full_results[full_results$adj.P.Val < 0.05 & abs(full_results$logFC) > 2.0,]
        if (length(new_results)>0){
            new_results <- subset(new_results, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.ID","Gene.symbol","Gene.title"))
            new_results <- subset(new_results, !is.na(Gene.symbol))
            new_results <- subset(new_results, Gene.symbol != "")
            deg_list <- unique(na.omit(new_results$Gene.ID))
            if (length(deg_list)>0){
                write.csv(new_results,file="FilteredDE.csv")
                write(deg_list,"DEGenes.txt")
                # summarize test results as "up", "down" or "not expressed"
                dT<-limma::decideTests(fit2, lfc = 2.0,p.value=0.05,method="separate",adjust.method="BH")
                pdf("Plots.pdf")
                # Venn diagram of results
                limma::vennDiagram(dT, circle.col=palette())
                # create Q-Q plot for t-statistic
                t.good <- which(!is.na(fit2$F)) # filter out bad probes
                limma::qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")
                #get volcano plot
                vplot <- full_results %>%
                  dplyr::mutate(Significant = adj.P.Val < 0.05, abs(logFC) > 2.0 ) %>%
                  dplyr::mutate(Rank = 1:dplyr::n(), Label = ifelse(Rank < 20, Gene.symbol,"")) %>%
                  ggplot2::ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + ggplot2::geom_point() + ggrepel::geom_text_repel(col="black")
                print(vplot)
                # Build histogram of P-values for all genes. Normal test
                # assumption is that most genes are not differentially expressed.
                tT2 <- limma::topTable(fit2, adjust="BH", sort.by="p", number=Inf)
                hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj", ylab = "Number of genes", main = "P-adj value distribution")
                # volcano plot (log P-value vs log fold change)
                ct <- 1
                limma::volcanoplot(fit2, coef=1, main=colnames(fit2)[ct], pch=20, highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))
                # MD plot (log fold change vs mean log expression)
                # highlight statistically significant (p-adj < 0.05) probes
                limma::plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
                abline(h=0)
                ################################################################
                # General expression data analysis
                ex <- exprs(gset2)
                # box-and-whisker plot
                ord <- order(gs)  # order samples by group
                palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02", "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
                par(mar=c(7,4,2,1))
                title <- paste (GSEaccession, "/", annotation(gset2), sep ="")
                boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
                legend("topleft", groups, fill=palette(), bty="n")
                # expression value distribution
                par(mar=c(4,4,2,1))
                title <- paste (GSEaccession, "/", annotation(gset2), " value distribution", sep ="")
                limma::plotDensities(ex, group=gs, main=title, legend ="topright")
                # mean-variance trend, helps to see if precision weights are needed
                limma::plotSA(fit2, main=paste("Mean variance trend, ",GSEaccession,sep=""))
                # get top genes
                if (length(deg_list)>20){
                    #get top 20 genes
                    topN <- 20
                    ##
                    ids_of_interest <- dplyr::mutate(new_results, Rank = 1:dplyr::n()) %>%
                    dplyr::filter(Rank < topN) %>%
                    dplyr::pull(ID)
                    gene_names <- dplyr::mutate(new_results, Rank = 1:dplyr::n()) %>%
                    dplyr::filter(Rank < topN) %>%
                    dplyr::pull(Gene.symbol)
                    ## Get the rows corresponding to ids_of_interest and all columns
                    tryCatch({
                        gene_matrix <- exprs(gset2)[ids_of_interest,]
                        pheatmap::pheatmap(gene_matrix, labels_row = gene_names, scale="row")
                        top20genesx <- new_results[new_results$ID %in% ids_of_interest,]
                        top20genes <- subset(top20genesx, select=c("ID","logFC","Gene.ID","Gene.symbol","Gene.title"))
                        top20genes$DE <- ""
                        top20genes$DE[top20genes$logFC > 0] <- "Up"
                        top20genes$DE[top20genes$logFC < 0] <- "Down"
                        suppressMessages(write.table(top20genes, file="Top20DEGenes.txt", quote=FALSE,sep="\t",row.names = FALSE))
                    }, error=function(e)
                    {
                        print("Error: Unable to create heatmap.")
                    })

                }
                dev.off()
            }
            else {
                print("** No protein-coding genes (DEG) identified **")
                new_results <- data.frame(matrix(ncol = 3, nrow = 0))
            }
        } else {
            print ("** No differentially expressed genes found **")
            new_results <- data.frame(matrix(ncol = 3, nrow = 0))
        }
    } else {
        print ("Sample grouping cannot be detected.")
        new_results <- data.frame(matrix(ncol = 3, nrow = 0))
    }
    return (new_results)
}
