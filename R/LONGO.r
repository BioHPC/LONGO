#' Use LONGO through shiny interface
#'
#' This function allows the user to input data files and alter the input variables to make sure the formatting is correct.
#' They can then run thte LONGO package which will output the results and plots in the browser and allow the user to download results as needed.
#'
#' @importFrom biomaRt listDatasets
#' @importFrom biomaRt useMart
#' @importFrom biomaRt listAttributes
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
#' @export
LONGO <- function() {
  ensembl <-
    useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")
  datasets <- biomaRt::listDatasets(ensembl)
  datasets <- datasets[1]
  datasets <- datasets[order(datasets[1]), ]
  datasets <- as.data.frame(datasets)
  status <- 0
  ui <- shiny::shinyUI(navbarPage(
    "LONGO",
    shiny::tabPanel(title = "Data input",
                    shiny::sidebarLayout(
                      shiny::sidebarPanel(
                        shiny::selectInput(
                          inputId = "species",
                          label = "Choose species gene ensembl:",
                          choices = datasets[1],
                          selected = "hsapiens_gene_ensembl"
                        ),
                        shiny::uiOutput("identifier"),
                        shiny::fileInput(
                          inputId = "datafile",
                          label = "Choose CSV file:",
                          accept = c(
                            "text/csv",
                            "text/comma-separated-values,text/plain",
                            ".csv",
                            ".tsv"
                          )
                        ),
                        shiny::checkboxInput(
                          inputId = "header",
                          label = "Header",
                          value = TRUE
                        ),
                        shiny::radioButtons(
                          inputId = "sep",
                          label = "Separator",
                          choices = c(
                            Tab = "\t",
                            Comma = ",",
                            Semicolon = ";"
                          ),
                          selected = "\t"
                        ),
                        shiny::checkboxInput(
                          inputId = "normalized",
                          label = "Normalized",
                          value = TRUE
                        ),
                        shiny::checkboxInput(
                          inputId = "filtered",
                          label = "Filtered",
                          value = TRUE
                        ),
                        shiny::actionButton(inputId = "action", label = "Submit")
                        # uiOutput("ui.action"),
                        # tags$hr(),
                        # uiOutput("ui.action.text")#,
                        #   #textOutput(outputId="timer")
                      ),
                      shiny::mainPanel(
                        shiny::uiOutput(outputId = "statusmessage"),
                        DT::dataTableOutput(outputId = "data_preview")
                      )
                    )),
    shiny::tabPanel(title = "data table",
                    shiny::mainPanel(
                      DT::dataTableOutput("data_annotated"),
                      shiny::downloadButton(outputId = "downloadRawData", label = "Download")
                    )),
    shiny::tabPanel(title = "plot output",
                    shiny::sidebarLayout(
                      shiny::sidebarPanel(
                        shiny::radioButtons(
                          inputId = "meanmedian",
                          label = "Sliding Window Method:",
                          choices = c("median", "mean")
                        ),
                        shiny::radioButtons(
                          inputId = "highestmean",
                          label = "Handle Multi Probes to Single Gene",
                          choices = c("highest", "mean")
                        ),
                        shiny::radioButtons(
                          inputId = "scale",
                          label = "X-axis scale:",
                          choices = c("linear", "log")
                        ),
                        shiny::radioButtons(
                          inputId = "legend",
                          label = "Legend Position",
                          choices = c("topleft", "topright", "bottomleft", "bottomright"),
                          selected = "topright"
                        ),
                        shiny::sliderInput(
                          inputId = "bin_size",
                          label = "Bin Size:",
                          min = 100,
                          max = 1000,
                          value = 200,
                          step = 100
                        ),
                        shiny::sliderInput(
                          inputId = "step_size",
                          label = "Step Size",
                          min = 20,
                          max = 200,
                          value = 40,
                          step = 10
                        ),
                        shiny::radioButtons(inputId = "control", label = "Control Column",choices = c(1,2)),
                        shiny::downloadButton(outputId = "downloadFinalData", label = "Download Data")
                      ),
                      shiny::mainPanel(
                        shiny::plotOutput(
                          outputId = "plot1",
                          height = 800,
                          width = 1200#,
                          # click = "plot_click"
                        ),
                        shiny::plotOutput(
                          outputId = "plot2",
                          height = 800,
                          width = 1200
                        ),
                        shiny::downloadButton(outputId = "downloadPData", label = "Download P Values"),
                        shiny::plotOutput(
                          outputId = "plot3",
                          height = 800,
                          width = 1200
                        ),
                        shiny::downloadButton(outputId = "downloadJSData", label = "Download JS Values"),
                     #   shiny::plotOutput(
                     #    outputId = "plot4",
                     #    height = 800,
                     #    width = 1200
                     #   ),
                        DT::dataTableOutput(outputId = "data_final_table")
                      )
                    )),
    shiny::tabPanel(title = "statistic plots",
                    shiny::mainPanel(# plotOutput(outputId = "plot2"),
                      # plotOutput(outputId = "plot3"),
                      shiny::plotOutput(outputId = "plot4")))
  ))

  server <- shiny::shinyServer(function(input, output, session) {
    output$identifier <- shiny::renderUI(if (is.null(input$species)) {
      return(NULL)
    }
    else {
      alldata.df$species_ensembl <-
        biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                         host = "www.ensembl.org",
                         dataset = input$species)
      ensembl_attributes <- biomaRt::listAttributes(alldata.df$species_ensembl)
      ensembl_attributes <- ensembl_attributes[1]
      shiny::selectInput(
        inputId = "attribute",
        label = "Choose gene identifier:",
        choices = unique(ensembl_attributes[1]),
        selected = "external_gene_name"
      )
    })

    options(shiny.maxRequestSize = 30 * 1024 ^ 2)


    alldata.df <- shiny::reactiveValues(
      filedata = NULL,
      rawdata = NULL,
      finaldata = NULL,
      P_data = NULL,
      JS_data = NULL,
      status = "Please select options and start. Analysis can take up to a minute to complete. Please be patient",
      labels = NULL,
      species_ensembl = NULL
    )

    output$columns <- shiny::renderText(c(1,2,3))

    shiny::observeEvent(list(input$datafile, input$sep, input$header), {
      if (is.null(input$datafile)) {
        return()
      }
      alldata.df$status <-
        "Please select options and start. Analysis can take up to a minute to complete. Please be patient"
      file1 <- input$datafile
      alldata.df$filedata <-
        read.csv(
          file = file1$datapath,
          header = input$header,
          sep = input$sep,
          comment.char = "!",
          na.strings = c("NA", " ", "")
        )
    })

    output$statusmessage <- shiny::renderUI({
      shiny::tags$h2(alldata.df$status)

    })

    shiny::observeEvent(input$action, {
      if(is.null(input$datafile)){
        alldata.df$status <- "Please load a file and then click submit"
        return()
      }
      message("calling biomart")
      temp2 <-
        callbiomaRt(alldata.df$filedata, input$attribute, alldata.df$species_ensembl)
      message("calling dict")
      alldata.df$rawdata <- dict(alldata.df$filedata, temp2)
      temp1 <- analyze(
        alldata.df$rawdata,
        input$highestmean,
        input$bin_size,
        input$step_size,
        input$meanmedian,
        input$filtered,
        input$normalized,
        2
      )

      alldata.df$finaldata <- as.data.frame(temp1[1])
      alldata.df$P_data <- as.data.frame(temp1[2])
      alldata.df$JS_data <- as.data.frame(temp1[3])
      updateRadioButtons(session = session,inputId = "control", choices = colnames(alldata.df$filedata)[2:(ncol(alldata.df$filedata))])
      alldata.df$status <- "Analysis completed"
    })

    shiny::observeEvent(list(
      input$bin_size,
      input$step_size,
      input$meanmedian,
      input$highestmean,
      input$control
    ),
    {
      if (is.null(alldata.df$rawdata)) {
        return()
      }
      temp1 <- analyze(
        alldata.df$rawdata,
        input$highestmean,
        input$bin_size,
        input$step_size,
        input$meanmedian,
        input$filtered,
        input$normalized,
        (which(colnames(alldata.df$rawdata)==input$control))
      )
      alldata.df$finaldata <- as.data.frame(temp1[1])
      alldata.df$P_data <- as.data.frame(temp1[2])
      alldata.df$JS_data <- as.data.frame(temp1[3])
    })

    output$data_preview <- DT::renderDataTable({
      DT::datatable(alldata.df$filedata)
    })

    output$ui.action <- shiny::renderUI({
      if (is.null(input$datafile)) {
        return (NULL)
      }
      shiny::actionButton(inputId = "action", label = "Submit")
    })

    output$data_annotated <- DT::renderDataTable(alldata.df$rawdata)

    output$plot1 <- shiny::renderPlot({
      data.df.analyzed <- alldata.df$finaldata
      if (input$scale == "linear") {
        x_vals <- (data.df.analyzed$kb_length)
        x_lab <- "Gene length in KB"
      }
      else{
        # log
        x_vals <- log(data.df.analyzed$kb_length)
        x_lab <- "Log(Gene Length)"
      }
      ymax <-
        (max(data.df.analyzed[, 2:(ncol(data.df.analyzed))]) * 1.2)
      yminim <- min(data.df.analyzed[, 2:(ncol(data.df.analyzed))])
      # png("LONGO_out.png", width=6, height=6, units="in", res=300)
      matplot(
        x = x_vals,
        y = data.df.analyzed[, 2],
        type = "l",
        col = 1,
        xlab = x_lab,
        ylim = c(yminim, ymax),
        ylab = "expression Level",
        main = "LONGO Plot"
      )
      for (i in 3:ncol(data.df.analyzed)) {
        par(new = TRUE)
        matplot(
          x = x_vals,
          y = data.df.analyzed[, i],
          type = "l",
          col = i - 1,
          xlab = "",
          ylab = "",
          ylim = c(yminim, ymax),
          axes = FALSE
        )
      }
      labels <- colnames(data.df.analyzed)
      legend(
        input$legend,
        legend = c(labels[2:length(labels)]),
        col = 2:ncol(data.df.analyzed) - 1,
        lty = 1,
        cex = 1,
        ncol = 2
      )
    })
    output$plot2 <- shiny::renderPlot({
      data.df.analyzed <- alldata.df$P_data

      if (input$scale == "linear") {
        x_vals <- (data.df.analyzed$kb_length)
        x_lab <- "Gene length in KB"
      }
      else{
        # log
        x_vals <- log(data.df.analyzed$kb_length)
        x_lab <- "Log(Gene Length)"
      }
      ymax <-
        (max(data.df.analyzed[, 2:(ncol(data.df.analyzed))]) * 1.2)
      yminim <- min(data.df.analyzed[, 2:(ncol(data.df.analyzed))])
      matplot(
        x = x_vals,
        y = data.df.analyzed[, 2],
        type = "l",
        col = 1,
        xlab = x_lab,
        ylim = c(yminim, ymax),
        ylab = "P value",
        main = "LONGO P Value Plot"
      )
      for (i in 3:ncol(data.df.analyzed)) {
        par(new = TRUE)
        matplot(
          x = x_vals,
          y = data.df.analyzed[, i],
          type = "l",
          col = i - 1,
          xlab = "",
          ylab = "",
          ylim = c(yminim, ymax),
          axes = FALSE
        )
      }
      labels <- colnames(data.df.analyzed)
      legend(
        input$legend,
        legend = c(labels[2:length(labels)]),
        col = 2:ncol(data.df.analyzed) - 1,
        lty = 1,
        cex = 1,
        ncol = 2
      )
    })

    output$plot3 <- shiny::renderPlot({
      data.df.analyzed <- alldata.df$JS_data
      if (input$scale == "linear") {
        x_vals <- (data.df.analyzed$kb_length)
        x_lab <- "Gene length in KB"
      }
      else{
        # log
        x_vals <- log(data.df.analyzed$kb_length)
        x_lab <- "Log(Gene Length)"
      }
      ymax <-
        max(data.df.analyzed[, 2:(ncol(data.df.analyzed))]) * 1.2
      yminim <- min(data.df.analyzed[, 2:(ncol(data.df.analyzed))])
      matplot(
        x = x_vals,
        y = data.df.analyzed[, 2],
        type = "l",
        col = 1,
        xlab = x_lab,
        ylim = c(yminim, ymax),
        ylab = "JS-distance",
        main = "LONGO JS Plot"
      )
      for (i in 3:ncol(data.df.analyzed)) {
        par(new = TRUE)
        matplot(
          x = x_vals,
          y = data.df.analyzed[, i],
          type = "l",
          col = i - 1,
          xlab = "",
          ylab = "",
          ylim = c(yminim, ymax),
          axes = FALSE
        )
      }
      labels <- colnames(data.df.analyzed)
      legend(
        input$legend,
        legend = c(labels[2:length(labels)]),
        col = 2:ncol(data.df.analyzed) - 1,
        lty = 1,
        cex = 1,
        ncol = 2
      )
    })
    output$plot4 <- shiny::renderPlot({
      data.df.analyzed <- alldata.df$finaldata
      last_point <- 10
      total_points <- (nrow(data.df.analyzed) - last_point)
      control_column <- 2
      labels <- colnames(data.df.analyzed)
      temp.df <- data.frame(data.df.analyzed[1:total_points, 1])
      for (i in 2:(ncol(data.df.analyzed))) {
        temp.df[,i] <- rep(0, total_points)
      }
      for (i in total_points:1) {
        for (j in 2:(ncol(data.df.analyzed))) {
          temp.df[i, j] <-
            cor(data.df.analyzed[(nrow(data.df.analyzed):i), control_column], data.df.analyzed[(nrow(data.df.analyzed):i), j])
        }
      }
      if (input$scale == "linear") {
        x_vals <- (temp.df[,1])
        x_lab <- "Gene length in KB"
      }
      else{
        # log
        x_vals <- log(temp.df[,1])
        x_lab <- "Log(Gene Length)"
      }
      ymax <-
        (max(temp.df[, 2:(ncol(temp.df))]) * 1.5)
      yminim <- min(temp.df[, 2:(ncol(temp.df))])
      matplot(
        x = x_vals,
        y = temp.df[, 2],
        type = "l",
        col = 1,
        xlab = x_lab,
        ylim = c(yminim, ymax),
        ylab = "Correlation",
        main = "LONGO Correlation Plot"
      )
      for (i in 3:ncol(temp.df)) {
        par(new = TRUE)
        matplot(
          x = x_vals,
          y = temp.df[, i],
          type = "l",
          col = i - 1,
          xlab = "",
          ylab = "",
          ylim = c(yminim, ymax),
          axes = FALSE
        )
      }
      legend(
        input$legend,
        legend = c(labels[2:length(labels)]),
        col = 2:ncol(temp.df) - 1,
        lty = 1,
        cex = 1,
        ncol = 2
      )
    })


  #  output$plot4 <- shiny::renderPlot({
  #  })

    output$data_final_table <-
      DT::renderDataTable(DT::datatable(alldata.df$finaldata))

    #  allow user to download simplified data with a click. not automate it
    output$downloadRawData <- shiny::downloadHandler(
      filename = function() {
        paste("LONGO_out_raw_", input$datafile, sep = "")
      },
      content = function(file) {
        write.csv(x = alldata.df$rawdata,
                  file = file,
                  row.names = FALSE)
      }
    )

    #   allow user to download simplified data with a click. not automate it
    output$downloadFinalData <- shiny::downloadHandler(
      filename = function() {
        paste("LONGO_out_final_", input$datafile, sep = "")
      },
      content = function(file) {
        write.csv(x = alldata.df$finaldata,
                  file = file,
                  row.names = FALSE)
      }
    )

    output$downloadPData <- shiny::downloadHandler(
      filename = function() {
        paste("LONGO_out_P_Values_", input$datafile, sep = "")
      },
      content = function(file) {
        write.csv(x = alldata.df$P_data,
                  file = file,
                  row.names = FALSE)
      }
    )

    output$downloadJSData <- shiny::downloadHandler(
      filename = function() {
        paste("LONGO_out_JS_Values_", input$datafile, sep = "")
      },
      content = function(file) {
        write.csv(x = alldata.df$JS_data,
                  file = file,
                  row.names = FALSE)
      }
    )

    session$onSessionEnded(stopApp)
  })

  shiny::shinyApp(ui, server)
}

