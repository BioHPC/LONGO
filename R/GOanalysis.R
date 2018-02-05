#' LONGO GO analysis
#'
#' This is an internal function used to do the GO analysis.
#'
#' @param analyzeData  dataframe to be analyzed
#' @param biomaRt.go.df dataframe of gene id to GO term
#' @param speciesEnsembl Ensembl used to call the bioMart if biomaRt.df is null
#' @param repNames Order/naming of the replications
#' @param GODomain Option for what GO domain to use for the analysis
#' @param GOScoring Option for what method should be used for the GO scoring
#' @param GoGraphing Option for what method should be used for the GO graphig
#' @param sigCount Option for the number of significant terms are needed in
#' the nodes to count
#' @param sigCutOff Option for the p value to mark significance
#' @param textScale Option to alter text size in GO graph
#' @import edgeR
#' @importFrom biomaRt getBM
#' @import topGO
#' @importFrom methods new
#' @return returns the GO data
#send in rawdata from alldata.df$rawdata as analyzeData
GOanalysis <- function(analyzeData,
                        biomaRt.go.df,
                        speciesEnsembl,
                        repNames,
                        GODomain,
                        GOScoring,
                        GoGraphing,
                        sigCount=15,
                        sigCutOff=0.01,
                        textScale=0.15) {
    library(topGO, quietly=TRUE)
    #query biomart for go terms
    if(is.null(biomaRt.go.df)){
        genes.unique <- unique(analyzeData$symbol)
        biomaRt.go.df <- biomaRt::getBM(
            attributes=c("external_gene_name","go_id"),
            filters="external_gene_name",
            values=genes.unique,
            mart=speciesEnsembl
        )
        geneToGoTerm <- split(biomaRt.go.df$go_id,
                            biomaRt.go.df$external_gene_name)
        biomaRt.go.df <- geneToGoTerm

    }
    else{
        geneToGoTerm <- biomaRt.go.df
    }

    analyzeData <- data.table::as.data.table(analyzeData)
    temp2 <- analyzeData[, lapply(.SD, mean, na.rm= TRUE),
                        by=c(colnames(analyzeData)[(ncol(analyzeData)-1)]),
                        .SDcols=c(2:(ncol(analyzeData)-2))]
    temp2 <-  as.data.frame(temp2)
    symbolnames <- as.data.frame(temp2$symbol)
    analyzeData <- temp2[,c(2:ncol(temp2))]#c(3:ncol(temp2))
    rownames(analyzeData) <- symbolnames$`temp2$symbol`


    final.go.df <- as.data.frame(analyzeData)


    for(i in 1:ncol(final.go.df)) {
        final.go.df[,i] <- as.numeric(final.go.df[,i])
    }

    y <- edgeR::DGEList(counts=final.go.df[1:ncol(final.go.df)],
                        group=(unlist(repNames))
                    )

    y <- edgeR::calcNormFactors(y)

    y <- edgeR::estimateDisp(y)
    et <- edgeR::exactTest(y)

    testing.p <- et$table$PValue
    names(testing.p) <- (symbolnames$`temp2$symbol`)

    topDiffGenes <- function(allScore) {
        return(allScore < sigCutOff)
    }

    testGOdata <- new("topGOdata",
                    description="Testing if it'll work", ontology=GODomain,
                    allGenes = testing.p, geneSel = topDiffGenes,
                    nodeSize = 10,
                    annot = annFUN.gene2GO, gene2GO = geneToGoTerm
                )

    resultKS <- topGO::runTest(testGOdata, algorithm=GoGraphing,
                                statistic=GOScoring
                            )

    par(cex = 1/as.numeric(as.character(textScale)))
    topGO::showSigOfNodes(testGOdata, score(resultKS),
                            firstSigNodes=sigCount, useInfo='def'
                        )
    recordedPlot = recordPlot()
    dev.off()

    return(list(testGOdata,recordedPlot,biomaRt.go.df))
}
