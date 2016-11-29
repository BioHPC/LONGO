#' Use LONGO through shiny interface
#'
#' This function allows the user to input data files and alter the input
#' variables to make sure the formatting is correct.
#' They can then run the LONGO package which will output the results and plots
#' in the browser and allow the user to download results as needed.
#'
#' @return Returns nothing
#' @importFrom biomaRt listDatasets
#' @importFrom biomaRt useMart
#' @importFrom biomaRt listAttributes
#' @importFrom biomaRt listFilters
#' @import shiny
#' @importFrom DT renderDataTable
#' @importFrom DT dataTableOutput
#' @importFrom DT datatable
#' @importFrom grDevices dev.off
#' @importFrom grDevices png
#' @importFrom graphics legend
#' @importFrom graphics matplot
#' @importFrom graphics par
#' @importFrom stats cor
#' @importFrom stats median
#' @importFrom utils read.csv
#' @importFrom utils write.csv
#' @importFrom utils write.table
#' @import graphics
#' @examples
#' if(interactive()) {LONGO()}
#' @export
LONGO <- function() {

    ensembl <-biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
        host="www.ensembl.org")
    datasets <- biomaRt::listDatasets(ensembl)
    datasets <- datasets[1]
    datasets <- datasets[order(datasets[1]),]
    datasets <- as.data.frame(datasets)
    ui <- shiny::shinyUI({
        navbarPage("LONGO",
            shiny::tabPanel(title="Data input",
                shiny::sidebarLayout(
                    shiny::sidebarPanel(
                        shiny::selectInput(inputId="species",
                            label="Choose species gene ensembl:",
                            choices=datasets[1],
                            selected="hsapiens_gene_ensembl"
                        ),
                        shiny::uiOutput("identifier"),
                        shiny::fileInput(inputId="datafile",
                            label="Choose CSV file:",
                            accept=c("text/csv",
                                "text/comma-separated-values,text/plain",
                                ".csv", ".tsv")
                            ),
                        shiny::checkboxInput(inputId="header",
                            label="Header",
                            value=TRUE
                        ),
                        shiny::radioButtons(inputId="sep",
                            label="Separator",
                            choices=c(
                                Tab="\t",
                                Comma=",",
                                Semicolon=";"),
                            selected="\t"
                        ),
                        shiny::checkboxInput(inputId="normalize",
                            label="Normalize",
                            value=TRUE
                        ),
                        shiny::checkboxInput(inputId="filter",
                            label="Filter for RNAseq",
                            value=FALSE
                        ),
                        shiny::actionButton(inputId="action",
                            label="Submit"
                        )
                    ),
                    shiny::mainPanel(
                        shiny::uiOutput(outputId="statusmessage"),
                        DT::dataTableOutput(outputId="data_preview")
                    )
                )
            ),
            shiny::tabPanel(title="Data Table",
                shiny::mainPanel(
                    DT::dataTableOutput("data_annotated"),
                    shiny::downloadButton(outputId="downloadRawData",
                        label="Download"
                    )
                )
            ),
            shiny::tabPanel(title="LONGO Output",
                shiny::sidebarLayout(
                    shiny::sidebarPanel(
                        shiny::radioButtons(inputId="meanmedian",
                            label="Sliding Window Method:",
                            choices=c("median", "mean"), selected="median"
                        ),
                        shiny::radioButtons(inputId="highestmean",
                            label="Handle Multi Probes to Single Gene",
                            choices=c("mean", "highest"), selected="mean"
                        ),
                        shiny::radioButtons(inputId="scale",
                            label="X-axis scale:",
                            choices=c("linear", "log"), selected="linear"
                        ),
                        shiny::radioButtons(inputId="legend",
                            label="Legend Position",
                            choices=c("topleft", "topright", "bottomleft",
                                "bottomright"),
                                selected="topright"
                        ),
                        shiny::sliderInput(inputId="bin_size",
                            label="Bin Size:",
                            min=100,
                            max=1000,
                            value=200,
                            step=100
                        ),
                        shiny::sliderInput(inputId="step_size",
                            label="Step Size",
                            min=20,
                            max=200,
                            value=40,
                            step=10
                        ),
                        shiny::radioButtons(inputId="control",
                            label="Control Column",
                            choices=c(1, 2)
                        )
                    ),
                    shiny::mainPanel(
                        shiny::plotOutput(outputId="plot1",
                            height=800,
                            width=1200
                        ),
                        shiny::downloadButton(outputId="downloadFinalData",
                            label="Download Data"
                        ),
                        shiny::plotOutput(outputId="plot2",
                            height=800,
                            width=1200
                        ),
                        shiny::downloadButton(outputId="downloadPData",
                            label="Download P Values"
                        ),
                        shiny::plotOutput(outputId="plot3",
                            height=800,
                            width=1200
                        ),
                        shiny::downloadButton(outputId="downloadDivergData",
                            label="Download Divergence Values"
                        ),
                        DT::dataTableOutput(outputId="data_final_table")
                    )
                )
            ),
            shiny::tabPanel(title="Long Gene Quotient",
                shiny::mainPanel(
                    shiny::plotOutput(outputId="plot4"),
                    shiny::downloadButton(outputId="downloadLongGeneQValues",
                        label="Download Long Gene Quotient Values"
                    )
                )
            )
        )
    })

    server <- shiny::shinyServer(function(input, output, session) {
        output$identifier <- shiny::renderUI(
        if (is.null(input$species)) {
            return(NULL)
        }
        else {
            alldata.df$species_ensembl <-
            biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org",
                dataset=input$species)
            ensemblAttributes <- biomaRt::listAttributes(
                alldata.df$species_ensembl)
            ensemblFilters <- biomaRt::listFilters(
                alldata.df$species_ensembl)
            ensemblIdentifiers <- merge(ensemblFilters[1],ensemblAttributes[1])
            shiny::selectInput(inputId="attribute",
                label="Choose gene identifier:",
                choices=unique(ensemblIdentifiers),
                selected="external_gene_name"
            )
        })

        options(shiny.maxRequestSize=30 * 1024 ^ 2)

        alldata.df <- shiny::reactiveValues(
            filedata=NULL,
            rawdata=NULL,
            finaldata=NULL,
            P_data=NULL,
            diverg_data=NULL,
            LQ_data=NULL,
            status="Please select options and start. Analysis can take up
                to a minute to complete. Please be patient",
            labels=NULL,
            species_ensembl=NULL
        )

        output$columns <- shiny::renderText(c(1, 2, 3))

        shiny::observeEvent(list(input$datafile, input$sep, input$header), {
            if (is.null(input$datafile)) {
                return()
            }
            alldata.df$status <-"Please select options and start. Analysis can
                take up to a minute to complete. Please be patient"
            file1 <- input$datafile
            alldata.df$filedata <-
            read.csv(
                file=file1$datapath,
                header=input$header,
                sep=input$sep,
                comment.char="!",
                na.strings=c("NA", " ", "")
            )
        })

        output$statusmessage <- shiny::renderUI({
            shiny::tags$h2(alldata.df$status)
        })

        shiny::observeEvent(input$action, {
                if (is.null(input$datafile)) {
                    alldata.df$status <- "Please load a file and
                        then click submit"
                    return()
                }
            temp2 <- callbiomaRt(alldata.df$filedata, input$attribute,
                alldata.df$species_ensembl)
            if(is.null(temp2)){
                alldata.df$status <-"The gene identifier or the species
                selected was incorrect for the uploaded data. Please select
                the correct species and gene identifier for your data."
                return()
            }
            alldata.df$rawdata <- dict(alldata.df$filedata, temp2)
            temp1 <- analyze(alldata.df$rawdata, input$highestmean,
                input$bin_size, input$step_size, input$meanmedian,
                input$filter, input$normalize, 2
            )
            if(is.null(temp1)){
                alldata.df$status <-"There were less than 200 genes identified
                with the species and gene identifier. Please make sure you have
                the correct inputs for your data and resubmit."
                return()
            }

            alldata.df$finaldata <- as.data.frame(temp1[1])
            alldata.df$P_data <- as.data.frame(temp1[2])
            alldata.df$diverg_data <- as.data.frame(temp1[3])
            alldata.df$LQ_data <- as.data.frame(temp1[4])
            updateRadioButtons(session=session,
                inputId="control",
                choices =
                colnames(alldata.df$filedata)[2:(ncol(alldata.df$filedata))]
            )
            alldata.df$status <- "Analysis completed"
        })

        shiny::observeEvent(
            list(input$bin_size, input$step_size, input$meanmedian,
                input$highestmean, input$control
        ),{
            if (is.null(alldata.df$rawdata)) {
                return()
            }
            temp1 <- analyze(alldata.df$rawdata, input$highestmean,
                input$bin_size, input$step_size, input$meanmedian,
                input$filter, input$normalize,
                    (which(
                        colnames(alldata.df$rawdata) == input$control
                    ))
            )
            alldata.df$finaldata <- as.data.frame(temp1[1])
            alldata.df$P_data <- as.data.frame(temp1[2])
            alldata.df$diverg_data <- as.data.frame(temp1[3])
            alldata.df$LQ_data <- as.data.frame(temp1[4])
        })

        output$data_preview <- DT::renderDataTable({
            DT::datatable(alldata.df$filedata)
        })

        output$ui.action <- shiny::renderUI({
            if (is.null(input$datafile)) {
                return (NULL)
            }
            shiny::actionButton(inputId="action", label="Submit")
        })

        output$data_annotated <- DT::renderDataTable(alldata.df$rawdata)

        output$plot1 <- shiny::renderPlot({
            data.df.analyzed <- alldata.df$finaldata
            if (input$scale == "linear") {
                x_vals <- (data.df.analyzed$kb_length)
                x_lab <- "Gene length (kb)"
            }
            else{
                # log
                x_vals <- log(data.df.analyzed$kb_length)
                x_lab <- "Log(Gene length (kb))"
            }
            ymax <-(max(data.df.analyzed[, 2:(ncol(data.df.analyzed))]) * 1.2)
            yminim <- min(data.df.analyzed[, 2:(ncol(data.df.analyzed))])
            # png("LONGO_out.png", width=6, height=6, units="in", res=300)
            matplot(x=x_vals, y=data.df.analyzed[, 2], type="l",
                col=1, xlab=x_lab, ylim=c(yminim, ymax),
                ylab="Gene expression (a.u.)", main="LONGO Plot"
            )
            for (i in 3:ncol(data.df.analyzed)) {
                par(new=TRUE)
                matplot(x=x_vals, y=data.df.analyzed[, i],
                    type="l", col=i - 1, xlab="",
                    ylab="", ylim=c(yminim, ymax),
                    axes=FALSE
                )
            }

            labels <- colnames(data.df.analyzed)
            legend(input$legend, legend=c(labels[2:length(labels)]),
                col=2:ncol(data.df.analyzed) - 1, lty=1,
                cex=1, ncol=2
            )
        })

        output$plot2 <- shiny::renderPlot({
            data.df.analyzed <- alldata.df$P_data

            if (input$scale == "linear") {
                x_vals <- (data.df.analyzed$kb_length)
                x_lab <- "Gene length (kb)"
            }
            else{
                # log
                x_vals <- log(data.df.analyzed$kb_length)
                x_lab <- "Log(Gene Length)"
            }
            ymax <-(max(data.df.analyzed[, 2:(ncol(data.df.analyzed))]) * 1.2)
            yminim <- min(data.df.analyzed[, 2:(ncol(data.df.analyzed))])
            matplot(x=x_vals, y=data.df.analyzed[, 2], type="l",
                col=1, xlab=x_lab, ylim=c(yminim, ymax),
                ylab="P value", main="LONGO P Value Plot"
            )
            for (i in 3:ncol(data.df.analyzed)) {
                par(new=TRUE)
                matplot(x=x_vals, y=data.df.analyzed[, i],
                    type="l", col=i - 1, xlab="",
                    ylab="", ylim=c(yminim, ymax),
                    axes=FALSE
            )}
            labels <- colnames(data.df.analyzed)
            legend(input$legend, legend=c(labels[2:length(labels)]),
                col=2:ncol(data.df.analyzed) - 1, lty=1,
                    cex=1, ncol=2
            )
        })

        output$plot3 <- shiny::renderPlot({
            data.df.analyzed <- alldata.df$diverg_data
            if (input$scale == "linear") {
                x_vals <- (data.df.analyzed$kb_length)
                x_lab <- "Gene length (kb)"
            }
            else{
                # log
                x_vals <- log(data.df.analyzed$kb_length)
                x_lab <- "Log(Gene Length)"
            }
            ymax <-max(data.df.analyzed[, 2:(ncol(data.df.analyzed))]) * 1.2
            yminim <- min(data.df.analyzed[, 2:(ncol(data.df.analyzed))])
            matplot(x=x_vals, y=data.df.analyzed[, 2], type="l",
                col=1, xlab=x_lab, ylim=c(yminim, ymax),
                ylab="Partial JS distance",
                main="LONGO Divergence Plot"
            )
            for (i in 3:ncol(data.df.analyzed)) {
                par(new=TRUE)
                    matplot(x=x_vals, y=data.df.analyzed[, i],
                        type="l", col=i - 1, xlab="", ylab="",
                        ylim=c(yminim, ymax), axes=FALSE
                    )
            }
            labels <- colnames(data.df.analyzed)
            legend(input$legend, legend=c(labels[2:length(labels)]),
                col=2:ncol(data.df.analyzed) - 1, lty=1,
                cex=1, ncol=2
            )
        })

        output$plot4 <- shiny::renderPlot({
            plot_data.df <- alldata.df$LQ_data
            bp <- barplot( as.matrix(plot_data.df), axes=1,
                ylim=c(-1, 1), axisnames=FALSE
            )
            abline(h=0.25, col="red")
            text(bp, par("usr")[3], labels=colnames(plot_data.df),
                srt=45, adj=c(1.1, 1.1), xpd=TRUE, cex=.9
            )
            mtext(text=input$datafile, outer=TRUE, cex=1.5)
        })

        output$data_final_table <-
            DT::renderDataTable(DT::datatable(alldata.df$finaldata))

        output$downloadRawData <- shiny::downloadHandler(
            filename=function() {
                paste("LONGO_out_raw_", input$datafile, sep="")
            },
            content=function(file) {
                write.csv(x=alldata.df$rawdata, file=file,
                    row.names=FALSE)
            }
        )

        output$downloadFinalData <- shiny::downloadHandler(
            filename=function() {
                paste("LONGO_out_final_", input$datafile, sep="")
            },
            content=function(file) {
                write.csv(x=alldata.df$finaldata, file=file,
                    row.names=FALSE)
            }
        )

        output$downloadPData <- shiny::downloadHandler(
            filename=function() {
                paste("LONGO_out_P_Values_", input$datafile, sep="")
            },
            content=function(file) {
                write.csv(x=alldata.df$P_data, file=file,
                    row.names=FALSE)
            }
        )

        output$downloadDivergData <- shiny::downloadHandler(
            filename=function() {
                paste("LONGO_out_Divergence_Distance_", input$datafile, sep="")
            },
            content=function(file) {
                write.csv(x=alldata.df$diverg_data, file=file,
                    row.names=FALSE)
            }
        )

        output$downloadLongGeneQValues <- shiny::downloadHandler(
            filename=function() {
                paste("LONGO_out_Long_Gene_Quotient_Values_", input$datafile, sep="")
            },
            content=function(file) {
                write.csv(x=alldata.df$LQ_data, file=file,
                          row.names=FALSE)
            }
        )

        session$onSessionEnded(stopApp)
    })

    shiny::shinyApp(ui, server)
}
