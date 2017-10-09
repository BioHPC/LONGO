#' LONGO GO analysis
#'
#' This is an internal function used to do the GO analysis.
#'
#' @param analyzeData  dataframe to be analyzed
#' @param biomaRt.df
#' @param speciesEnsembl
#' @param repNames column of the data used as a control for the
#' statistical analysis
#' @import edgeR
#' @import biomaRt
#' @import topGO
#' @return returns the GO
 ##imports topGO?
 #send in rawdata from alldata.df$rawdata as analyzeData
 GOanalysis <- function(analyzeData,
                        biomaRt.go.df,
                        speciesEnsembl,
                        repNames) {
# add num sig nodes
# add option for type of analysis
# add option for stat test
# add p value option
# add options for doing just long genes..


     library(topGO,quietly = TRUE)

        #query biomart for go terms

##    n <- proc.time()
    if(is.null(biomaRt.go.df)){
        genes.unique <- unique(analyzeData$symbol)
        biomaRt.go.df <- biomaRt::getBM(
            attributes=c("external_gene_name","go_id"),
            filters="external_gene_name",
            values=genes.unique,
            mart=speciesEnsembl
        )
        geneToGoTerm <- split(biomaRt.go.df$go_id, biomaRt.go.df$external_gene_name)
        #print(paste("null:",head(geneToGoTerm)))
        print(typeof(geneToGoTerm))
    }
    else{
        geneToGoTerm <- biomaRt.go.df
        #print(paste("not null:",head(geneToGoTerm)))
        print(typeof(geneToGoTerm))
    }

print("yooooo head(analyzeData)=")
print(head(analyzeData))

     analyzeData <- data.table::as.data.table(analyzeData)
     temp2 <- analyzeData[, lapply(.SD, mean, na.rm= TRUE),
                          by=c(colnames(analyzeData)[(ncol(analyzeData)-1)]),
                          .SDcols=c(2:(ncol(analyzeData)-2))]
     temp2 <-  as.data.frame(temp2)
     symbolnames <- as.data.frame(temp2$symbol)
     analyzeData <- temp2[,c(2:ncol(temp2))]#c(3:ncol(temp2))
     rownames(analyzeData) <- symbolnames$`temp2$symbol`

    # get fold change

     print("yooooo head(analyzeData)again2=")
     print(head(analyzeData))
    final.go.df <- as.data.frame(analyzeData)
    # data should be in cols 2:(ncol()-2)
#    for(i in 2:(ncol(analyzeData)-2)){
#      final.go.df[i] <- log10(analyzeData[i]/analyzeData[controlColumn])
#      print(colnames(analyzeData)[i])
#    }

#head(final.go.df[3])
#    final.go.df[3] <- as.numeric(final.go.df[3])
    #remove first identifier
    # final.go.df <- final.go.df[-1]
    # #remove length
    # final.go.df <- final.go.df[-ncol(final.go.df)]
    # symbolnames <- final.go.df$symbol
    # rownames(final.go.df) <- symbolnames
    # # reorder it
    # final.go.df <- final.go.df[-ncol(final.go.df)]


    for(i in 1:ncol(final.go.df)) {final.go.df[,i] <- as.numeric(final.go.df[,i])}
    y <- edgeR::DGEList(counts=final.go.df[1:ncol(final.go.df)], group=unlist(repNames))
#    print(y$samples)
    #keep <- rowSums(edgeR::cpm(y)>1) >= 2
    #y <- y[keep, , keep.lib.sizes=FALSE]
    y <- edgeR::calcNormFactors(y)
#    print(y$samples)
    y <- edgeR::estimateDisp(y)
    et <- edgeR::exactTest(y)
###   edgeR::topTags(et)

    testing.p <- et$table$PValue
    #testing.p <- as.data.frame(et$table$PValue)
    #testing.p$`et$table$PValue`<-as.numeric(testing.p$`et$table$PValue`)
    names(testing.p) <- (symbolnames$`temp2$symbol`)

#    testing.p <- as.list(et$table$PValue)
#    names(testing.p) <- symbolnames$`temp2$symbol`
#    testing.p<-as.numeric(testing.p)
    #types of domains...
    #BP biological process
    #MF molecular function
    #CC cellular component

     topDiffGenes <- function(allScore) {
         return(allScore < 0.01)
     }

     testGOdata <- new("topGOdata",
                       description="Testing if it'll work", ontology="BP",
                       allGenes = testing.p, geneSel = topDiffGenes,
                       nodeSize = 10,
                       annot = annFUN.gene2GO, gene2GO = geneToGoTerm)
##### #    resultFisher <- runTest(testGOdata, algorithm="classic", statistic="fisher")
 #    resultFisher

     resultKS <- runTest(testGOdata, algorithm="classic", statistic="ks")
#     resultKS

#######     resultKS.elim <- runTest(testGOdata, algorithm="elim", statistic="ks")
#     resultKS.elim

     # pValue.classic <- score(resultKS)
     # pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
     # gstat <- termStat(testGOdata, names(pValue.classic))
     # gSize <- gstat$Annotated / max(gstat$Annotated) * 4
     # colMap <- function(x) {
     #     .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
     #     return(.col[match(1:length(x), order(x))])
     # }
     # gCol <- colMap(gstat$Significant)
     # plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",
     #      pch = 19, cex = gSize, col = gCol)
    showSigOfNodes(testGOdata, score(resultKS), firstSigNodes = 15, useInfo = 'all')
    recordedPlot = recordPlot()
    dev.off()

   return(list(testGOdata,recordedPlot,biomaRt.go.df))
#    return()
}
