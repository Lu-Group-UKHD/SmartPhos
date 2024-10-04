# This script contains the server part for the Phosphoproteomics Shiny app.
# The server script acts based on the inputs from the UI script.

library(shiny)
library(shinythemes)
library(shinyjs)
library(shinyBS)
library(shinyWidgets)
library(DT)
library(plotly)
library(matrixStats)
library(piano)
library(dplyr)


# Increasing the limit of upload file size from 5 MB to 50 MB
options(shiny.maxRequestSize=500*1024^2)

#temporary directory for unzip raw data folder
outDir <- "rawFolder"

shinyServer(function(input, output, session) {

  # a reactive val to store the multiassayexperiment object
  mae <- reactiveVal()

  # a reactive val to store the mae after normalization adjustment
  # this is used for updating PP/FP ratio plot
  maeAdj <- reactiveVal()

  # a reactive variable to stored the processed and unfiltered assay
  processedDataUF <- reactiveVal()

  # a reactive variable to stored the processed assay
  processedData <- reactiveVal()

  # a reactive value to store the saved list
  saveList <- reactiveValues(file = list.files("save/"))

  # a reactive values variable to store the value of all inputs
  inputsValue <- reactiveValues()

  # a box showing saved results
  output$saveListBox <- renderUI({
    selectInput("seleFile", label = NULL, saveList$file,
                size = 5, selectize = FALSE)
  })

  # save calculated results
  observeEvent(input$save, {
    saveObj <- mae()
    fileName <- paste0(input$text, "_", format(Sys.Date(), "%Y%m%d"), ".Rds")
    saveRDS(saveObj, file = paste0("save/", fileName))
    saveList$file <- unique(c(saveList$file, fileName))
  })

  # remove saved file
  observeEvent(input$remove, {
    fileName <- input$seleFile
    file.remove(paste0("save/", fileName))
    saveList$file <- saveList$file[saveList$file != fileName]
  })

  # load saved file
  observeEvent(input$load, {
    mae(readRDS(paste0("save/", input$seleFile)))
    maeAdj(NULL)
  })

  # download the multiAssayExperiment object
  output$download <- downloadHandler(
    filename = function()
        { paste0(input$text, "_", format(Sys.Date(), "%Y%m%d"), ".Rds") },
    content = function(file) {
      saveRDS(mae(), file = file)
    })

  # zip file
  observeEvent(input$uploadZip, {
    inputsValue$upload <- input$upload
    # removing the already existing directory before unzipping
    unlink(outDir, recursive = TRUE)
    dir.create(outDir, showWarnings = FALSE)
    unzip(input$uploadZip$datapath, exdir = outDir, junkpaths = TRUE)

    # error check
    fileTable <- NULL
    tryCatch({
      fileTable <- as.data.frame(read.delim(file.path(outDir, "fileTable.txt")))
    },
    error = function(e) {
      showModal(modalDialog(
        title = "Processing zip files failed...",
        "fileTable.txt file not found.",
        easyClose = TRUE,
        footer = NULL
      ))})
    if(!is.null(fileTable))
    {
      tryCatch({
        stopifnot(c("sampleType", "id", "fileName") %in% colnames(fileTable))
        # rendering column annotation option for column annotations
        output$colAnnoBoxPreprocess <- renderUI({
          # excluding type since it's already represented by two assays
          selectInput("colAnnoPreprocess", "Select additional column annotations:",
                      colnames(fileTable)[colnames(fileTable) != "type"],
                      selected = colnames(fileTable)[colnames(fileTable) != "type"], multiple = TRUE)
        })
      },
      error = function(e) {
        showModal(modalDialog(
          title = "Missing columns in fileTable.txt file...",
          "Please make sure fileTable.txt file contains columns with column names: fileName, sampleType and id.",
          easyClose = TRUE,
          footer = NULL
        ))})
    }
  })

  # process raw files
  observeEvent(input$processUpload, {
    withProgress(message = 'Processing files', {
      fileTable <- as.data.frame(read.delim(file.path(outDir, "fileTable.txt")))
      #fileTable$fileName <- sub("^", "./rawFolder/", fileTable$fileName)
      #outDir is a variable
      fileTable$fileName <- file.path(outDir, fileTable$fileName)

      tryCatch({
        inputsValue$tool <- input$tool
        # for data from Spectronaut
        if (input$tool == "Spectronaut") {
          testData <- SmartPhos::readExperimentDIA(fileTable, annotation_col = input$colAnnoPreprocess)
        }
        # for data from MaxQuant
        else if (input$tool == "MaxQuant") {
          testData <- SmartPhos::readExperiment(fileTable, annotation_col = input$colAnnoPreprocess)
        }
        mae(testData)
        maeAdj(NULL) #reset maeAdj
      },
      error = function(e) {
        showModal(modalDialog(
          title = "Processing of the uploaded data failed...",
          "Contact the SmartPhos package developers as there might me some issue with the SmartPhos package.",
          easyClose = TRUE,
          footer = NULL
        ))})
    })
  })

  # read the uploaded object
  observeEvent(input$uploadObject, {
    inputsValue$upload <- input$upload
    file <- input$uploadObject
    ext <- tools::file_ext(file$datapath)
    req(file)
    # making sure an rds object is being uploaded
    validate(need(ext %in% c("rds", "RDS", "Rds"), "Please upload an rds object"))
    mae(readRDS(file$datapath))
    maeAdj(NULL) #reset maeAdj
  })

  # rendering select options in the UI based on the uploaded object
  output$option <- renderUI({
    if (!is.null(mae())) {
      selectInput("assay", "Select one assay", names(experiments(mae())))
    }
  })

  # rendering options for the phospho-enriched or non-enriched sample type
  output$seleAssayBox <- renderUI({
    if (!is.null(mae())) {
      inputsValue$assay <- input$assay
      if (input$assay == "Proteome") {
        radioButtons("getPP", "Select the sample type",
                     c("Phospho-enriched" = TRUE, "Non-enriched" = FALSE),
                     selected = FALSE, inline = TRUE)
      }
      else if (input$assay == "Phosphoproteome") {
        radioButtons("getFP", "Select the sample type",
                     c("Phospho-enriched" = FALSE, "Non-enriched" = TRUE),
                     selected = FALSE, inline = TRUE)
      }
    }
  })

  # loaded object for boxplot and table
  loadedData <- reactive({
    se <- mae()[[input$assay]]
    colData(se) <- colData(mae()[, colnames(se)])
    if (input$assay == "Phosphoproteome") {
      inputsValue$getFP <- input$getFP
      if (input$getFP) {
        ppe <- se[,se$sampleType %in% c("FullProteome", "FP")]
        colData(ppe) <- colData(se)[colnames(ppe),]
      }
      else {
        ppe <- se[,se$sampleType %in% c("Phospho", "PP")]
        colData(ppe) <- colData(se)[colnames(ppe),]
      }
      ppe
    }
    else if (input$assay == "Proteome") {
      inputsValue$getPP <- input$getPP
      if (input$getPP) {
        fpe <- se[,se$sampleType %in% c("Phospho", "PP")]
        colData(fpe) <- colData(se)[colnames(fpe),]
      }
      else {
        fpe <- se[,se$sampleType %in% c("FullProteome", "FP")]
        colData(fpe) <- colData(se)[colnames(fpe),]
      }
      fpe
    }
  })


  #-----------------------------------------------------------------------------
  # launching MatrixQCvis from the phosphoproteomics app
  # observeEvent(input$launch_app, {
  #   mae <- readRDS(input$upload$datapath)
  #   se <- mae[[input$assay]]
  #   colData(se) <- colData(mae)
  #   job_env <-  new.env()
  #   job_env$se <- se
  #   rstudioapi::jobRunScript(path = "script.R", importEnv = TRUE)
  # })
  observeEvent(input$launch_app, {
    showModal(modalDialog(
      title = "Launching MatrixQCvis not possible...",
      "This functionality is yet to be implemented in the Shiny App.",
      easyClose = TRUE,
      footer = NULL
    ))
  })
  #-----------------------------------------------------------------------------

  # rendering options for the phospho-enriched or non-enriched sample type

  # whether to normalize phospho data with proteomic data
  output$ifNormByProteinBox <- renderUI({
      if (!is.null(mae())) {
          if (input$assay == "Phosphoproteome" && "Proteome" %in% names(assays(mae()))) {
              checkboxInput("ifNormByProtein","Normalize phosphorylation by protein expression", value = FALSE)
          }
      }
  })

  output$seleNormCorrect <- renderUI({
    if (!is.null(mae())) {
      if (input$assay == "Phosphoproteome" && input$getFP == FALSE) {
        checkboxInput("ifNormCorrect","Perform normalization correction (may take long time, no further normalization needed)", value = FALSE)
      }
    }
  })

  output$ifAlreadyNormBox <- renderUI({
      if(input$ifNormCorrect) {
          updateRadioButtons(session = getDefaultReactiveDomain(), inputId = "normalize", selected = FALSE)
          radioButtons("ifAlreadyNorm", "Have the data already been normalized by Spectronaut/MaxQuant ?", c("Yes","No"), selected = "No", inline = TRUE)
      }
  })

  # render select option for selecting column(s) for batch effects removal
  output$seleColBatch <- renderUI({
    if(input$batch) {
      selectizeInput("colBatch", "Choose up to 2 columns:",
                     colnames(colData(mae())),
                     multiple = TRUE,
                     options = list(maxItems = 2))
    }
  })

  observeEvent(input$processSelection, {
    withProgress(message = 'Processing files', {
      # normalization correction if selected
      # first step, whether to perform phospho normalization correction,
      # if phospho data is selected.
      maeData <- mae()

      if (input$assay == "Phosphoproteome") {
        if (!is.null(input$ifNormCorrect) && input$ifNormCorrect) {
          if (("Proteome" %in% names(mae())) && ("FullProteome" %in% unique(colData(mae())$sampleType))) {
            if (!is.null(colData(mae())$sampleName)) {
              inputsValue$ifNormCorrect <- input$ifNormCorrect
              maeData <- runPhosphoAdjustment(mae(),
                                              normalization = ifelse(input$ifAlreadyNorm == "Yes", FALSE, TRUE), # depends on whether the data were already normalized
                                              minOverlap = 3, # at least three overlapped feature between sample pair
                                              completeness = 0.5 # use feature present in at least 50% of the samples
              )
              assays(maeData[["Phosphoproteome"]])[["Intensity"]] <- assays(maeData[["Phosphoproteome"]])[["Intensity_adjusted"]]
              assays(maeData[["Phosphoproteome"]])[["Intensity_adjusted"]] <- NULL
              #maeData <- maeData[, !is.na(maeData$adjustFactorPP)]
              maeAdj(maeData)
            } else {
              showModal(modalDialog(
                title = "Normalization correction not possible!",
                "Please make sure that the uploaded R object or the fileTable.txt file contains the sampleName column with appropriate entry. Add sampleName column and upload the object or zip file again or continue without normalization correction.",
                easyClose = TRUE,
                footer = NULL
              ))}
          } else {
            showModal(modalDialog(
              title = "Normalization correction not possible!",
              "Please make sure that the uploaded R object or the fileTable.txt file contains the proteomics for the unenriched samples. Add the appropriate samples or continue without normalization correction.",
              easyClose = TRUE,
              footer = NULL
            ))}
          }
        # whether normalize the phospho intensity by protein expression,
        # if normalization adjustment is performed, this will be performed after normalization adjustment
        if (!is.null(input$ifNormByProtein) && input$ifNormByProtein) {
          if (("Proteome" %in% names(mae())) && ("FullProteome" %in% unique(colData(mae())$sampleType))) {
            if (!is.null(colData(mae())$sampleName)) {
              inputsValue$ifNormByProtein <- input$ifNormByProtein
              maeData <- SmartPhos::normByFullProteome(maeData)
            } else {
              showModal(modalDialog(
                title = "Normalize phosphorylation by protein expression!",
                "Please make sure that the uploaded R object or the fileTable.txt file contains the sampleName column with appropriate entry. Add sampleName column and upload the object or zip file again or continue without this option.",
                easyClose = TRUE,
                footer = NULL
              ))}
          } else {
            showModal(modalDialog(
              title = "Normalize phosphorylation by protein expression!",
              "Please make sure that the uploaded R object or the fileTable.txt file contains the proteomics for the unenriched samples. Add the appropriate samples or continue without this option.",
              easyClose = TRUE,
              footer = NULL
            ))}
        }
      }
      # summarizedAssayExperiment object of the selected assay
      se <- maeData[[input$assay]]
      colData(se) <- colData(maeData[, colnames(se)])
      # make sure NUll is passed to batch argument when user have not
      # selected Batch correction
      if (input$batch) batchCol <- input$colBatch else batchCol <- NULL
      inputsValue$batch <- input$batch
      inputsValue$colBatch <- batchCol

      if (input$assay == "Proteome") {
        # function from the utils.R script for the preprocessing of the data
        fp <- preprocessProteome(se, filterList = NULL,
                                 transform = input$transform,
                                 normalize = input$normalize,
                                 getPP = input$getPP,
                                 missCut = input$missFilter,
                                 removeOutlier = strsplit(input$outliers, ",\\s*")[[1]],
                                 impute = input$impute,
                                 batch = batchCol,
                                 scaleFactorTab = NULL)
        processedDataUF(fp)
      }
      else {
        pp <- preprocessPhos(se, filterList = NULL,
                             transform = input$transform,
                             normalize = ifelse(!is.null(input$ifNormCorrect) && input$ifNormCorrect ==FALSE, input$normalize, FALSE), #if normalization correction has been performed, normalization should not be performed again
                             getFP = input$getFP,
                             missCut = input$missFilter,
                             removeOutlier = strsplit(input$outliers, ",\\s*")[[1]],
                             impute = input$impute,
                             batch = batchCol,
                             scaleFactorTab = NULL)
        processedDataUF(pp)
      }

      # log
      inputsValue$transform <- input$transform
      inputsValue$normalize <- input$normalize
      inputsValue$missFilter <- input$missFilter
      inputsValue$outliers <- input$outliers
      inputsValue$impute <- input$impute
    })
  })

  # render download link if the summarized experiment object is created
  output$downloadSE <- renderUI({
    if (!is.null(processedDataUF())) {
      downloadButton('downloadSEobj', 'Download the proceesed assay',
                     style="color: #FFFFFF; background-color: #3498DB; border-color: #2E86C1;")
    }
  })

  # download summarized experiment object
  output$downloadSEobj <- downloadHandler(
    filename = function()
        { paste0("summarizedExp", "_", format(Sys.Date(), "%Y%m%d"), ".Rds") },
    content = function(file) {
      saveRDS(processedDataUF(), file = file)
    })

  # text output to show the number of samples and features
  output$dataInfo <- renderUI({
    if (!is.null(processedDataUF())) {
      HTML(sprintf("<b>Number of samples: %s<br/>Number of features: %s<b><br/>",
                   ncol(processedDataUF()), nrow(processedDataUF())))
    }
    else {
      HTML(sprintf("<b>Number of samples: %s<br/>Number of features: %s<b><br/>",
                   ncol(loadedData()), nrow(loadedData())))
    }
  })

  # data table
  output$metaData <- DT::renderDataTable({
    if (!is.null(processedDataUF())) {
      colDataTable <- mutate_if(data.frame(colData(processedDataUF())),
                                is.character, as.factor)
      datatable(colDataTable, filter = "top", rownames = FALSE,
                caption = "Subset the data for further analysis by selecting the annotations in the filtering boxes",
                selection = "none", style = "bootstrap")
    }
    else {
      colDataTable <- mutate_if(data.frame(colData(loadedData())), is.character, as.factor)
      datatable(colDataTable, filter = "top", rownames = FALSE,
                selection = "none", style = "bootstrap")
    }
  })

  # Reactive object for getting filtered samples
  processedData <- reactive({
    if(!is.null(processedDataUF())) {
      selectedRows <- as.character(data.frame(colData(processedDataUF()))[input$metaData_rows_all,] %>% pull(sample))
      dataNew <- processedDataUF()[, selectedRows]
      dataNew
    }
  })

  output$missingPlot <- renderPlot({
    plotMissing(loadedData())
  })

  # to toggle between show and hide functionality
  toggle("missingPlot")
  observeEvent(input$missingV, {
   toggle("missingPlot")
  })

  output$colorBoxUI <- renderUI({
    selectInput("colorBox", "Select color by:",
                c("none", colnames(colData(mae()))),
                selected = "none")
  })

  # plot log ratio
  output$boxPlotLogRatio <- renderPlot({
    if (is.null(maeAdj())) {
        plotLogRatio(mae(), normalization = FALSE)
    } else {
        plotLogRatio(maeAdj(), normalization = FALSE)
    }
  })

  # Plot boxplot
  output$boxPlot <- renderPlot({
    inputsValue$colorBox <- input$colorBox
    if (!is.null(processedDataUF())) {
      g <- plotIntensity(processedDataUF(), input$colorBox)
    }
    else {
      g <- plotIntensity(loadedData(), input$colorBox)
    }
    g
  })

  ####################################### PCA ##################################

  output$messagePCA <- renderText({
    if (is.null(processedData())) {
      "Make sure to press the process button inside the preprocessing options panel before proceeding."
    }
  })

  runPCA <- observeEvent(input$RunPCA, {
    if ("imputed" %in% assayNames(processedData())) {
      output$errMsgPCA <- renderText("")
      withProgress(message = "Running principal component analysis, please wait...", {
        pca <- stats::prcomp(t(assays(processedData())[["imputed"]]),
                             center = TRUE, scale.=TRUE)

        # download PCA values as tsv file
        output$downloadPCATable <- downloadHandler(
          filename = function() { paste('PCAvalues', '.tsv', sep='') },
          content = function(file) {
            write.table(pcaDf, file = file, quote = FALSE, sep = '\t',
                        col.names = NA)
        })
        # rendering select option for x-axis of PCA plot
        output$xaxisPCAui <- renderUI({
          selectInput("xaxisPCA", "Select the x-axis PC", colnames(pca[["x"]]),
                      selected = "PC1")
        })
        # rendering select option for y-axis of PCA plot
        output$yaxisPCAui <- renderUI({
          selectInput("yaxisPCA", "Select the y-axis PC", colnames(pca[["x"]]),
                      selected = "PC2")
        })
        # rendering color option for PCA plot
        output$colorPCAui <- renderUI({
          selectInput("colorPCA", "Select color by:",
                      c("none",colnames(colData(processedData()))),
                      selected = "none")
        })
        # rendering shape option for PCA plot
        output$shapePCAui <- renderUI({
          selectInput("shapePCA", "Select shape by:",
                      c("none",colnames(colData(processedData()))),
                      selected = "none")
        })

        plotpc <- reactive({
          inputsValue$xaxisPCA <- input$xaxisPCA
          inputsValue$yaxisPCA <- input$yaxisPCA
          inputsValue$colorPCA <- input$colorPCA
          inputsValue$shapePCA <- input$shapePCA

          SmartPhos::plotPCA(pca, processedData(), xaxis = input$xaxisPCA,
                  yaxis = input$yaxisPCA, color = input$colorPCA,
                  shape = input$shapePCA)

        })
        output$pcplot <- renderPlotly({
          ggplotly(plotpc())
        })
      })}

    else {
      showModal(modalDialog(
        title = "PCA not possible...",
        "Please make sure imputation is not selected none in the preprocessing tab.",
        easyClose = TRUE,
        footer = NULL
      ))
    }

    # download the PCA plot as PDF file
    output$downPCA <- downloadHandler(
      filename = function() { paste0("pca", '.pdf', sep = '') },
      content = function(file) {
        ggsave(file, plot = plotpc(),
               device = "pdf",
               width = input$figWidthPCA,
               height = input$figHeightPCA,
               limitsize = FALSE)
      }
    )
  })

  ################################# heat map ###################################

  # rendering column annotation option for heatmap
  output$colAnnoBoxHM <- renderUI({
    if (!is.null(processedData())) {
      selectInput("colAnnoHM", "Select additional column feature:",
                  colnames(colData(processedData())),
                  selected = NULL, multiple = TRUE)
    }
  })

  # render options for top variants option
  output$topGenes <- renderUI({
    numericInput("numGenes","Number of genes", value = 100)
  })

  output$colCluster <- renderUI({
    numericInput("numClustCol", "Number of column clusters", value = 1)
  })

  output$rowCluster <- renderUI({
    numericInput("numClustRow", "Number of row clusters", value = 1)
  })

  plotMap <- eventReactive(input$doPlot, {

    if (input$chooseType == "Top variant") {
      inputsValue$numGenes <- input$numGenes
      setName <- sprintf("Top %s most variant genes", input$numGenes)
      p <- plotHeatmap(type = input$chooseType,
                       se = processedData(),
                       top = input$numGenes,
                       annotationCol = input$colAnnoHM,
                       cutCol = input$numClustCol,
                       cutRow = input$numClustRow,
                       title = setName)
    }
    else if (input$chooseType == "Differentially expressed") {
      if(!is.null(filterDE())) {
        setName <- "Differentially expressed genes"
        p <- plotHeatmap(type = input$chooseType,
                         se = processedDataSub(),
                         data = filterDE(),
                         annotationCol = input$colAnnoHM,
                         cutCol = input$numClustCol,
                         cutRow = input$numClustRow,
                         clustCol = FALSE,
                         clustRow = FALSE,
                         title = setName)
      }
      else {
        showModal(modalDialog(
          title = "Plotting heatmap not possible...",
          "Make sure you have performed Differential expression analysis before.",
          easyClose = TRUE,
          footer = NULL
        ))
      }
    }
    else if (input$chooseType == "Selected time series cluster") {
      if(!is.null(selectedCluster())) {
        inputsValue$seleClusterHM <- input$seleCluster
        setName <- input$seleCluster
        p <- plotHeatmap(type = input$chooseType,
                         se = processedDataSubClust(),
                         data = selectedCluster(),
                         annotationCol = input$colAnnoHM,
                         cutCol = input$numClustCol,
                         cutRow = input$numClustRow,
                         clustCol = FALSE,
                         clustRow = FALSE,
                         title = setName)
      }
      else {
        showModal(modalDialog(
          title = "Plotting heatmap not possible...",
          "Make sure you have performed Time series clustering before.",
          easyClose = TRUE,
          footer = NULL
        ))
      }
    }
    inputsValue$chooseType <- input$chooseType
    inputsValue$colAnnoHM <- input$colAnnoHM
    inputsValue$numClustCol <- input$numClustCol
    inputsValue$numClustRow <- input$numClustRow
    p
  })

  # plot the heatmap
  output$plotHM <- renderPlot({
    if (!is.null(plotMap())) {
      output$errMsg1 <- renderText("")
      withProgress(message = "Plotting heatmap, please wait...", value = NULL, {
        grid.draw(plotMap())
      })
    } else {
      output$errMsg1 <- renderText("Please perform differential expression analysis first or load a previous result!")
    }
  })

  # download the Heatmap plot as PDF file
  output$downHM <- downloadHandler(
    filename = function() { paste0("heatmap", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotMap(),
             device = "pdf",
             width = input$figWidthHM,
             height = input$figHeightHM,
             limitsize = FALSE)
  })

  ############################ differential expression #########################


  # select sample by sample name
  output$ID1Box <- renderUI({
    selectInput("seleID1", "Select sample ID(s) for reference group",
                processedData()$sample, multiple = TRUE)
  })

  output$ID2Box <- renderUI({
    selectInput("seleID2", "Select sample ID(s) for target group",
                processedData()$sample, multiple = TRUE)
  })

  output$seleMetaColBoxDiff <- renderUI({
      useCol <- colnames(colData(processedData()))
      useCol <- useCol[!useCol %in% c("sample","sampleType","adjustFactorPP","sampleName")]
      selectInput("seleMetaColDiff", "Select metadata column for testing",
                  useCol, multiple = FALSE)
  })


  # the selection box to select treatment for the reference group
  output$treat1Box <- renderUI({
    allTreat <- unique(processedData()[[input$seleMetaColDiff]])
    selectInput("seleTreat1", "Select level(s) for the reference group",
                allTreat, multiple = TRUE)
  })

  # select time point for reference
  output$time1Box <- renderUI({
    if (!is.null(processedData()$timepoint)) {
      if (!is.null(input$seleTreat1))
        allTime <- unique(processedData()[,processedData()[[input$seleMetaColDiff]] %in% input$seleTreat1]$timepoint) else
          allTime <- unique(processedData()$timepoint)
        selectInput("seleTime1", "Select time point(s) for the reference group",
                    allTime, multiple = TRUE)
    }
  })

  # the selection box to select treatment for the target groups
  output$treat2Box <- renderUI({
    allTreat <- unique(processedData()[[input$seleMetaColDiff]])
    selectInput("seleTreat2", "Select levels(s) for the target group",
                allTreat, multiple = TRUE)
  })

  # the select the comparison time point
  output$time2Box <- renderUI({
    if (!is.null(processedData()$timepoint)) {
      if (!is.null(input$seleTreat2))
        allTime <- unique(processedData()[,processedData()[[input$seleMetaColDiff]] %in% input$seleTreat2]$timepoint) else
          allTime <- unique(processedData()$timepoint)
        selectInput("seleTime2", "Select time point(s) for the target group",
                    allTime, multiple = TRUE)
    }
  })

  output$seleMethodBox <- renderUI({
    allowChoice <- c("limma", "ProDA")
    radioButtons("deMethod", "Select DE method", allowChoice, inline = TRUE)
  })

  # a reactive value to store the se object after subsetting
  processedDataSub <- reactiveVal()
  # a reactive value to store the differential expression results
  tableDE <- reactiveVal()
  # a reactive value to monitor if plot histogram or boxplot
  ifHistogram <- reactiveValues(value = TRUE)
  # a reactive value to save the clicked row or volcano plot point value
  lastClicked <- reactiveValues()
  # a reactive value to sync the data table with the volcano plot selection
  colorRows <- reactiveValues(
    row_priority = c(),
    row_color = c()
  )
  # a reactive value to highlight the point with row selection
  ptHiglight <- reactiveValues(
    log2FC = NULL,
    pValue = NULL
  )

  # reactive event for calculating differential expression
  observeEvent(input$runDE, {
    withProgress(message = "Running DE analysis, please wait...", value = NULL, {
      tryCatch({
        if (input$seleID) {
          inputsValue$seleID1 <- input$seleID1
          inputsValue$seleID2 <- input$seleID2
          de <- performDifferentialExp(se = processedData(), assay = "Intensity",
                                       method = input$deMethod,
                                       reference = input$seleID1,
                                       target = input$seleID2)
        }
        else {
          inputsValue$seleMetaColDiff <- input$seleMetaColDiff
          inputsValue$Treat1 <- input$seleTreat1
          inputsValue$Treat2 <- input$seleTreat2
          if (!is.null(processedData()$timepoint)) {
            inputsValue$Time1 <- input$seleTime1
            inputsValue$Time2 <- input$seleTime2
            de <- performDifferentialExp(se = processedData(),
                                         assay = "Intensity",
                                         method = input$deMethod,
                                         condition = input$seleMetaColDiff,
                                         reference = input$seleTreat1,
                                         target = input$seleTreat2,
                                         refTime = input$seleTime1,
                                         targetTime = input$seleTime2)
          }
          else {
            de <- performDifferentialExp(se = processedData(),
                                         assay = "Intensity",
                                         method = input$deMethod,
                                         condition = input$seleMetaColDiff,
                                         reference = input$seleTreat1,
                                         target = input$seleTreat2)
          }
        }
        inputsValue$deMethod <- input$deMethod

        processedDataSub(de$seSub)
        tableDE(de$resDE)
        ifHistogram$value <- TRUE
      },
      error = function(e) {
        showModal(modalDialog(
          title = "Running DE analysis failed...",
          "Make sure you have selected the correct experiment options.",
          easyClose = TRUE,
          footer = NULL
        ))}
      )
    })
  })

  # filter options for pvalue, log2FC
  filterDE <- reactive({
    if (!is.null(tableDE())) {
      DEtab <- tableDE()
      if(input$ifAdjusted) {
        inputsValue$ifAdjusted <- input$ifAdjusted
        DEtab <- filter(DEtab, abs(log2FC) >= input$fcFilter,
                        padj <= as.numeric(input$pFilter))
      } else {
        DEtab <- filter(DEtab, abs(log2FC) >= input$fcFilter,
                        pvalue <= as.numeric(input$pFilter))
      }
      inputsValue$fcFilter <- input$fcFilter
      inputsValue$pFilter <- input$pFilter
      DEtab
    }
  })

  # table for differentially expressed genes
  output$DEtab <- DT::renderDataTable({
    if (!is.null(filterDE())) {
      if (!is.null(colorRows$row_color)) {
        datatable(filterDE(), selection = 'single', rownames = FALSE,
                  caption = "Differentially expressed genes") %>%
          formatStyle("ID",
                      target = "row",
                      backgroundColor = styleEqual(colorRows$row_priority,
                                                   colorRows$row_color,
                                                   default = 'white')) %>%
          formatRound(c('pvalue',"padj"),digits=3) %>%
          formatRound(c('log2FC','stat'),digits=2)
      }
      else {
        datatable(filterDE(), selection = 'single', rownames = FALSE,
                  caption = "Differentially expressed genes") %>%
          formatRound(c('pvalue',"padj"),digits=3) %>%
          formatRound(c('log2FC','stat'),digits=2)
      }
    }
  })

  output$downloadTableUI <- renderUI({
    if (!is.null(filterDE())) {
      downloadLink('downloadTable', 'Download current table')
    }
  })

  # a button to download DE gene table as tsv
  output$downloadTable <- downloadHandler(
    filename = function() { paste('DE_Table', '.tsv', sep = '') },
    content = function(file) {
      write.table(filterDE(), file = file, quote=FALSE, sep = '\t', col.names = NA)
    }
  )

  # volcano plot
  plotV <- reactive({
    v <- plotVolcano(tableDE = tableDE(),
                     pFilter = input$pFilter,
                     fcFilter = input$fcFilter)
    v
  })

  output$plotVolcano <- renderPlotly({
    p <- ggplotly(plotV(), source = "volcano")
    p %>%
      event_register("plotly_click")
    if (!is.null(ptHiglight$log2FC)) {
      # highlight the point selected either by clicking on the data table row or
      # by clicking on a point on the volcano plot
      p <- add_trace(p, x = ptHiglight$log2FC, y = ptHiglight$pValue,
                     type = "scatter", mode = 'markers',
                     marker = list(size = 10, symbol = "star"),
                     showlegend = FALSE)
    }
    p
  })

  # observe event when a point in the volcano plot is clicked
  observeEvent(event_data("plotly_click", source = "volcano"),{
    # a point in the volcano plot is clicked, turn the ifHistogram to false
    ifHistogram$value <- FALSE
    d <- event_data("plotly_click", source = "volcano")
    lastInfo <- d$customdata
    lastClicked$geneID <- filterDE()[filterDE()$ID == lastInfo,]$ID
    if (input$assay == "Phosphoproteome") {
      # site contains the information of the gene, residue and site position
      lastClicked$geneSymbol <- filterDE()[filterDE()$ID == lastInfo,]$site
    }
    else {
      lastClicked$geneSymbol <- filterDE()[filterDE()$ID == lastInfo,]$Gene
    }
    colorRows$row_priority <- filterDE()$ID
    colorRows$row_priority <- c(lastClicked$geneID, colorRows$row_priority[colorRows$row_priority != lastClicked$geneID])
    colorRows$row_color <- lapply(colorRows$row_priority, function(x) ifelse(x %in% lastClicked$geneID, "lightgreen", "white"))
    ptHiglight$log2FC <- filterDE()[filterDE()$ID == lastInfo,]$log2FC
    ptHiglight$pValue <- -log10(filterDE()[filterDE()$ID == lastInfo,]$pvalue)
  })
  # observe event when a row in the DE table is clicked
  observeEvent(input$DEtab_row_last_clicked,{
    # a row is clicked, turn the ifHistogram to false
    ifHistogram$value <- FALSE
    lastInfo <- input$DEtab_row_last_clicked
    lastClicked$geneID <- filterDE()[lastInfo,]$ID
    if (input$assay == "Phosphoproteome") {
      lastClicked$geneSymbol <- filterDE()[lastInfo,]$site
    }
    else {
      lastClicked$geneSymbol <- filterDE()[lastInfo,]$Gene
    }
    ptHiglight$log2FC <- filterDE()[lastInfo,]$log2FC
    ptHiglight$pValue <- -log10(filterDE()[lastInfo,]$pvalue)
    colorRows$row_priority <- c(lastClicked$geneID, colorRows$row_priority[colorRows$row_priority != lastClicked$geneID])
    colorRows$row_color <- lapply(colorRows$row_priority, function(x) ifelse(x %in% lastClicked$geneID, "lightgreen", "white"))
  })

  # a ui to hold the plot on the first panel
  output$ui.plot <- renderUI({
    if (!is.null(tableDE())) {
      #the size of the ui is depend on whether it's a histogram or a boxplot
      if(ifHistogram$value == TRUE) {
        fig.width <- 600
        fig.height <- 400
      } else {
        fig.width <- 600
        fig.height <- 500
      }
      plotlyOutput("plot1", width = paste0(fig.width,"px"),
                   height = paste0(fig.height,"px"))
    }
  })

  # Histogram plot or boxplot on the first panel
  output$plot1 <- renderPlotly({
    if (!is.null(tableDE())){
      if (ifHistogram$value == TRUE) {
        # plot the histogram of p values
        p <- ggplot(tableDE(), aes(x = pvalue)) +
          geom_histogram(fill = "grey", col = "blue", alpha=0.7) +
          ggtitle("P value histogram") + theme_bw() +
          theme(plot.title = element_text(hjust = 0.5))
        ggplotly(p) %>% config(displayModeBar = FALSE)
      }
      # Box-plot for comparison
      else {
        inputsValue$geneID_DE <- lastClicked$geneID
        inputsValue$geneSymbol_DE <- lastClicked$geneSymbol

        p <- plotBox(se = processedDataSub(), id = lastClicked$geneID,
                     symbol = lastClicked$geneSymbol)
      }
    }
  })
  ############################ time series clustering ##########################

  #### Widgets
  output$seleMetaColBoxTime <- renderUI({
    useCol <- colnames(colData(processedData()))
    useCol <- useCol[!useCol %in% c("sample","sampleType","adjustFactorPP", "timepoint")]
    selectInput("seleMetaColTime", "Select metadata column for testing",
                useCol, multiple = FALSE)
  })

  # the selection box to select the condition
  output$clusterTreatBox <- renderUI({
    allTreat <- unique(processedData()[[input$seleMetaColTime]])
    selectInput("seleTreat_cluster", "Select a condition", allTreat)
  })

  # the selection box to select the reference condition
  output$clusterTreatRefBox <- renderUI({
    inputsValue$seleMetaColTime <- input$seleMetaColTime
    allTreat <- unique(processedData()[[input$seleMetaColTime]])
    allTreat <- allTreat[!allTreat %in% input$seleTreat_cluster]
    selectInput("seleTreat_clusterRef","Select a reference condition", allTreat)
  })

  # link to download the cluster table
  output$downloadClusterTabUI <- renderUI({
    if( !is.null (clusterTabVal())) {
      downloadLink("downloadClusterTab","Download cluster table")
    }
  })

  output$downloadClusterTab <- downloadHandler(
    filename = function() { paste('clusterTable', '.tsv', sep='') },
    content = function(file) {
      allClusterFeature <- clusterTabVal() %>%
          distinct(feature, .keep_all = TRUE) %>% .$feature
      if (!is.null(rowData(processedData())$site)) {
        allClusterSite <- rowData(processedData())[allClusterFeature, "site"]
        allClusterSequence <- rowData(processedData())[allClusterFeature, "Sequence"]
      } else {
        allClusterSite <- rep(NA, length(allClusterFeature))
        allClusterSequence <- rep(NA,length(allClusterFeature))
      }
      clustTab <- clusterTabVal() %>% distinct(feature,.keep_all = TRUE) %>%
        mutate(Gene = rowData(processedData())[allClusterFeature, "Gene"],
               site = allClusterSite,
               Sequence = allClusterSequence) %>%
        select(any_of(c("feature", "Gene", "cluster", "prob", "cNum", "Sequence", "site"))) %>%
        dplyr::rename(id = feature, probability = prob,
                      clusterSize = cNum, PhosphoSite = site)
      write.table(clustTab, file = file, quote = FALSE, sep = '\t', col.names = NA)
    }
  )

  # reactive value to save all the timepoints
  allTime <- reactiveVal()

  # selecting time range
  output$timerangeBox <- renderUI({
    # list the time points available to the selected condition
    processedDataSub <- processedData()[, processedData()[[input$seleMetaColTime]] == input$seleTreat_cluster]
    allTimepoint <- unique(processedDataSub$timepoint)
    # remove time points with only 1 sample
    remove1sampleT <- c()
    for (time in allTimepoint) {
      idPerTime <- processedDataSub[, processedDataSub$timepoint == time]$sample
      if (length(idPerTime) < 2)
        remove1sampleT <- append(remove1sampleT, time)
    }
    allTimepoint <- allTimepoint[!allTimepoint %in% remove1sampleT]
    # if a reference treatment is selected, then only list timepoints shared between the two conditions
    if (input$clusterFor == "logFC" | input$clusterFor == "two-condition expression") {
      processedDataRef <- processedData()[, processedData()[[input$seleMetaColTime]] == input$seleTreat_clusterRef]
      timepointRef <- unique(processedDataRef$timepoint)
      # remove time points with only 1 sample
      remove1sampleT <- c()
      for (time in timepointRef) {
        idPerTime <- processedDataRef[, processedDataRef$timepoint == time]$sample
        if (length(idPerTime) < 2)
          remove1sampleT <- append(remove1sampleT, time)
      }
      timepointRef <- timepointRef[!timepointRef %in% remove1sampleT]
      # select only the intersection
      allTimepoint <- intersect(allTimepoint, timepointRef)
    }

    allTime(allTimepoint)

    checkboxGroupInput("seleTimeRange", "Time points to include",
                       allTimepoint, allTimepoint, inline = TRUE)
  })

  # UI for adding zero timepoint
  output$zeroTime <- renderUI({
    if (!('0min' %in% allTime())) {
      checkboxInput("addZero","Add zero timepoint", FALSE)
    }
  })

  # UI for conditions with zero timepoint
  output$zeroTreat <- renderUI({
    if (input$addZero) {
      cd <- colData(processedData())
      treatment <- unique(cd[[input$seleMetaColTime]][cd$timepoint == "0min"])

      output$zeroTreatInfo <- renderUI({
        if (input$addZero) {
          HTML(sprintf("Found %s treatments with 0min timepoints",
                       length(treatment)))
        }
      })

      selectInput("seleZeroTreat", "Select treatment for copying the values:",
                  treatment,selected = NULL)
    }
  })

  exprMatObj <- reactive({
    # Processing with selecting expression (i.e. only 1 condition)
    if (input$clusterFor == "expression") {
      inputsValue$clusterFor <- input$clusterFor
      inputsValue$seleTreat_cluster <- input$seleTreat_cluster
      inputsValue$seleTimeRange <- input$seleTimeRange
      if (!is.null(input$seleZeroTreat) && input$addZero) {
        inputsValue$addZero <- input$addZero
        inputsValue$seleZeroTreat <- input$seleZeroTreat
        processedDataSub <- addZeroTime(processedData(),
                                        input$seleMetaColTime,
                                        input$seleTreat_cluster,
                                        input$seleZeroTreat,
                                        input$seleTimeRange)
      }
      else {
        processedDataSub <- processedData()[, processedData()[[input$seleMetaColTime]] == input$seleTreat_cluster &
                                              processedData()$timepoint %in% input$seleTimeRange]
      }
      assayMat <- assay(processedDataSub)
      inputsValue$ifFilterFit <- input$ifFilterFit
      # Filtering genes with p values from the spline fit test
      if (input$ifFilterFit) {
        inputsValue$pSpline <- input$pSpline
        inputsValue$ifSplineFdr <- input$ifSplineFdr
        if (!is.null(processedDataSub$subjectID)) {
          assayMat <- splineFilter(assayMat,
                                   subjectID = processedDataSub$subjectID,
                                   time = processedDataSub$timepoint,
                                   df = length(unique(processedDataSub$timepoint))-1,
                                   pCut = as.numeric(input$pSpline),
                                   ifFDR = input$ifSplineFdr)
        } else {
          assayMat <- splineFilter(assayMat, subjectID = NULL,
                                   time = processedDataSub$timepoint,
                                   df = length(unique(processedDataSub$timepoint))-1,
                                   pCut = as.numeric(input$pSpline),
                                   ifFDR = input$ifSplineFdr)
        }
        processedDataSub <- processedDataSub[rownames(assayMat),]
      }
      exprMat <- lapply(unique(processedDataSub$timepoint), function(tp) {
        rowMedians(assayMat[,processedDataSub$timepoint == tp])
      }) %>% bind_cols() %>% as.matrix()

      rownames(exprMat) <- rownames(processedDataSub)
      colnames(exprMat) <- unique(processedDataSub$timepoint)

    }
    else if (input$clusterFor == "logFC") {
      inputsValue$clusterFor <- input$clusterFor
      inputsValue$seleTreat_cluster <- input$seleTreat_cluster
      inputsValue$seleTimeRange <- input$seleTimeRange
      if (!is.null(input$seleZeroTreat) && input$addZero) {
        inputsValue$addZero <- input$addZero
        inputsValue$seleZeroTreat <- input$seleZeroTreat
        processedDataSub <- processedData()[, processedData()[[input$seleMetaColTime]] == input$seleTreat_cluster]
        allTimepoint <- unique(processedDataSub$timepoint)
        processedDataRef <- processedData()[, processedData()[[input$seleMetaColTime]] == input$seleTreat_clusterRef]
        timepointRef <- unique(processedDataRef$timepoint)
        # add zero timepoint samples if missing
        if (!("0min" %in% allTimepoint) && !("0min" %in% timepointRef)) {
          processedDataSub <- addZeroTime(processedData(), input$seleMetaColTime,
                                          input$seleTreat_cluster,
                                          input$seleZeroTreat, input$seleTimeRange)
          processedDataRef <- addZeroTime(processedData(), input$seleMetaColTime,
                                          input$seleTreat_clusterRef,
                                          input$seleZeroTreat, input$seleTimeRange)
        }
        else if (!("0min" %in% allTimepoint) && ("0min" %in% timepointRef)) {
          processedDataSub <- addZeroTime(processedData(), input$seleMetaColTime,
                                          input$seleTreat_cluster,
                                          input$seleZeroTreat, input$seleTimeRange)
          processedDataRef <- processedData()[,processedData()[[input$seleMetaColTime]] == input$seleTreat_clusterRef &
                                                processedData()$timepoint %in% input$seleTimeRange]
        }
        else if (("0min" %in% allTimepoint) && !("0min" %in% timepointRef)) {
          processedDataSub <- processedData()[,processedData()[[input$seleMetaColTime]] == input$seleTreat_cluster &
                                                processedData()$timepoint %in% input$seleTimeRange]
          processedDataRef <- addZeroTime(processedData(), input$seleMetaColTime,
                                          input$seleTreat_clusterRef,
                                          input$seleZeroTreat, input$seleTimeRange)
        }
      }
      else {
        processedDataSub <- processedData()[,processedData()[[input$seleMetaColTime]] == input$seleTreat_cluster &
                                              processedData()$timepoint %in% input$seleTimeRange]
        processedDataRef <- processedData()[,processedData()[[input$seleMetaColTime]] == input$seleTreat_clusterRef &
                                              processedData()$timepoint %in% input$seleTimeRange]
      }
      assayMat <- assay(processedDataSub)
      RefMat <- assay(processedDataRef)
      # calculate fold change by subtracting assayMat to mean intensities of RefMat
      # here the mean intensities in RefMat are calculated per time point or per time point and subject ID.
      if (!is.null(processedData()$subjectID)) {
        fcMat <- lapply(unique(processedDataSub$timepoint), function(tp) {
          lapply(unique(processedDataSub$subjectID), function(id) {
            RefMean = rowMeans(RefMat[,processedDataRef$timepoint == tp &
                                        processedDataRef$subjectID == id])
            assayMat[,processedDataSub$timepoint == tp & processedDataSub$subjectID == id] - RefMean
          })
        }) %>% bind_cols() %>% as.matrix()
      } else {
        fcMat <- lapply(unique(processedDataSub$timepoint), function(tp) {
          RefMean = rowMeans(RefMat[,processedDataRef$timepoint == tp])
          assayMat[,processedDataSub$timepoint == tp] - RefMean
        }) %>% bind_cols() %>% as.matrix()
      }
      rownames(fcMat) <- rownames(assayMat)
      # rearrange columns in processedDataSub to match with fcMat for splineFilter and calculating mean logFC later
      processedDataSub <- processedDataSub[,colnames(fcMat)]
      #  apply spline filter
      inputsValue$ifFilterFit <- input$ifFilterFit
      if (input$ifFilterFit) {
        inputsValue$pSpline <- input$pSpline
        inputsValue$ifSplineFdr <- input$ifSplineFdr
        if (!is.null(processedDataSub$subjectID)) {
          fcMat <- splineFilter(fcMat, subjectID = processedDataSub$subjectID,
                                time = processedDataSub$timepoint,
                                df = length(unique(processedDataSub$timepoint))-1,
                                pCut = as.numeric(input$pSpline),
                                ifFDR = input$ifSplineFdr)
        } else { # i.e. if subjectID is not present
          fcMat <- splineFilter(fcMat, subjectID = NULL,
                                time = processedDataSub$timepoint,
                                df = length(unique(processedDataSub$timepoint))-1,
                                pCut = as.numeric(input$pSpline),
                                ifFDR = input$ifSplineFdr)
        }
        processedDataSub <- processedDataSub[rownames(fcMat),]
        exprMat <- lapply(unique(processedDataSub$timepoint), function(tp) {
          rowMeans(fcMat[,processedDataSub$timepoint == tp])
        }) %>% bind_cols() %>% as.matrix()
      } else { # i.e. if choose not to spline filter

        # Calculate the mean intensities per time point, THEN calculate the logFC

        assayMatMean <- lapply(unique(processedDataSub$timepoint), function(tp) {
          rowMeans(assayMat[,processedDataSub$timepoint == tp])
        }) %>% bind_cols() %>% as.matrix()
        RefMatMean <- lapply(unique(processedDataSub$timepoint), function(tp) {
          rowMeans(RefMat[,processedDataRef$timepoint == tp])
        }) %>% bind_cols() %>% as.matrix()
        exprMat <- assayMatMean -  RefMatMean
      }
      rownames(exprMat) <- rownames(processedDataSub)
      colnames(exprMat) <- unique(processedDataSub$timepoint)

    }
    else if (input$clusterFor == "two-condition expression") {
      inputsValue$clusterFor <- input$clusterFor
      inputsValue$seleTreat_cluster <- input$seleTreat_cluster
      inputsValue$seleTreat_clusterRef <- input$seleTreat_clusterRef
      inputsValue$seleTimeRange <- input$seleTimeRange
      if (!is.null(input$seleZeroTreat) && input$addZero) {
        inputsValue$addZero <- input$addZero
        inputsValue$seleZeroTreat <- input$seleZeroTreat
        processedDataSub <- processedData()[, processedData()[[input$seleMetaColTime]] == input$seleTreat_cluster]
        allTimepoint <- unique(processedDataSub$timepoint)
        processedDataRef <- processedData()[, processedData()[[input$seleMetaColTime]] == input$seleTreat_clusterRef]
        timepointRef <- unique(processedDataRef$timepoint)
        # add zero timepoint samples if missing
        if (!("0min" %in% allTimepoint) && !("0min" %in% timepointRef)) {
          processedData1 <- addZeroTime(processedData(), input$seleMetaColTime,
                                        input$seleTreat_cluster,
                                        input$seleZeroTreat, input$seleTimeRange)
          processedData2 <- addZeroTime(processedData(), input$seleMetaColTime,
                                        input$seleTreat_clusterRef,
                                        input$seleZeroTreat, input$seleTimeRange)
          assay <- cbind(assay(processedData1), assay(processedData2))
          cd <- rbind(colData(processedData1), colData(processedData2))
          emeta <- elementMetadata(processedData())

          processedDataSub <- SummarizedExperiment(assays=SimpleList(intensity=assay), colData = cd, rowData = emeta)
        }
        else if (!("0min" %in% allTimepoint) && ("0min" %in% timepointRef)) {
          processedData1 <- addZeroTime(processedData(), input$seleMetaColTime,
                                        input$seleTreat_cluster,
                                        input$seleZeroTreat, input$seleTimeRange)
          processedData2 <- processedData()[, processedData()[[input$seleMetaColTime]] == input$seleTreat_clusterRef &
                                              processedData()$timepoint %in% input$seleTimeRange]
          assay <- cbind(assay(processedData1), assay(processedData2))
          cd <- rbind(colData(processedData1), colData(processedData2))
          emeta <- elementMetadata(processedData())

          processedDataSub <- SummarizedExperiment(assays=SimpleList(intensity=assay), colData = cd, rowData = emeta)
        }
        else if (("0min" %in% allTimepoint) && !("0min" %in% timepointRef)) {
          processedData1 <- processedData()[, processedData()[[input$seleMetaColTime]] == input$seleTreat_cluster &
                                              processedData()$timepoint %in% input$seleTimeRange]
          processedData2 <- addZeroTime(processedData(), input$seleMetaColTime,
                                        input$seleTreat_clusterRef,
                                        input$seleZeroTreat, input$seleTimeRange)
          assay <- cbind(assay(processedData1), assay(processedData2))
          cd <- rbind(colData(processedData1), colData(processedData2))
          emeta <- elementMetadata(processedData())

          processedDataSub <- SummarizedExperiment(assays=SimpleList(intensity=assay), colData = cd, rowData = emeta)
        }
      }
      else {
        processedDataSub <- processedData()[,processedData()[[input$seleMetaColTime]] %in% c(input$seleTreat_cluster, input$seleTreat_clusterRef) &
                                              processedData()$timepoint %in% input$seleTimeRange]
      }
      assayMat <- assay(processedDataSub)
      inputsValue$ifFilterFit <- input$ifFilterFit
      if (input$ifFilterFit) {
        inputsValue$pSpline <- input$pSpline
        inputsValue$ifSplineFdr <- input$ifSplineFdr
        if (!is.null(processedDataSub$subjectID)) {
          assayMat <- splineFilter(assayMat, subjectID = processedDataSub$subjectID,
                                   time = processedDataSub$timepoint,
                                   treatment = processedDataSub[[input$seleMetaColTime]],
                                   refTreatment = input$seleTreat_clusterRef,
                                   df = length(unique(processedDataSub$timepoint))-1,
                                   pCut = as.numeric(input$pSpline),
                                   ifFDR = input$ifSplineFdr)
        } else {
          assayMat <- splineFilter(assayMat, subjectID = NULL,
                                   time = processedDataSub$timepoint,
                                   treatment = processedDataSub[[input$seleMetaColTime]],
                                   refTreatment = input$seleTreat_clusterRef,
                                   df = length(unique(processedDataSub$timepoint))-1,
                                   pCut = as.numeric(input$pSpline),
                                   ifFDR = input$ifSplineFdr)
        }
        processedDataSub <- processedDataSub[rownames(assayMat),]
      }
      processedDataSub$timeTreat <- paste0(processedDataSub$timepoint,"_",
                                           processedDataSub[[input$seleMetaColTime]])
      exprMat <- lapply(unique(processedDataSub$timeTreat), function(tt) {
        rowMedians(assayMat[,processedDataSub$timeTreat == tt])
      }) %>% bind_cols() %>% as.matrix()

      rownames(exprMat) <- rownames(processedDataSub)
      colnames(exprMat) <- unique(processedDataSub$timeTreat)
    }
    # filter based on variance
    sds <- apply(exprMat,1,sd)
    varPer <- as.numeric(input$topVarTime)
    inputsValue$topVarTime <- input$topVarTime
    exprMat <- exprMat[order(sds, decreasing = TRUE)[seq(1, varPer/100*nrow(exprMat))], ]

    # only center when it's for expression
    if (input$clusterFor != "logFC") exprMat <- mscale(exprMat)

    # remove NA values
    exprMat <- exprMat[complete.cases(exprMat), ]
    exprMat
  })

  clusterPlotVal <- reactiveVal()
  clusterTabVal <- reactiveVal()

  # (currently not in use) plot to find out optimal number of clusters
  observeEvent(input$plotSilhouette, {
    withProgress(message = "Calculating Silhouette and WSS scores, please wait...", {
      d <- exprMatObj()
      p1 <- factoextra::fviz_nbclust(d,  kmeans, c("silhouette"), k.max = 20)
      p2 <- factoextra::fviz_nbclust(d,  kmeans, c("wss"), k.max = 20)
    })

    clusterPlotVal(cowplot::plot_grid(p1,p2,NULL,NULL, ncol=2))
    clustNum(6)
  })

  # Subset se object for heatmap visualization
  processedDataSubClust <- reactive({
    if (input$clusterFor == "expression") {
      processedData.sub <- processedData()[,processedData()[[input$seleMetaColTime]] == input$seleTreat_cluster
                                           & processedData()$timepoint %in% input$seleTimeRange]
    }
    else {
      processedData.sub <- processedData()[,processedData()[[input$seleMetaColTime]] %in% c(input$seleTreat_cluster, input$seleTreat_clusterRef)
                                           & processedData()$timepoint %in% input$seleTimeRange]
    }
    processedData.sub
  })

  # plot time-series clustering result
  observeEvent(input$runCluster, {
    withProgress(message = "Performing cmeans clustering, please wait...", {
      # plot clustering result and error check
      tryCatch({
        inputsValue$seleNumCluster <- input$seleNumCluster
        inputsValue$seleProbCut <- input$seleProbCut
        if (input$clusterFor != "two-condition expression") {
          clusterRes <- clusterTS(exprMatObj(),
                                  as.numeric(input$seleNumCluster),
                                  pCut = as.numeric(input$seleProbCut))
        }
        else {
          clusterRes <- clusterTS(exprMatObj(),
                                  as.numeric(input$seleNumCluster),
                                  pCut = as.numeric(input$seleProbCut),
                                  twoCondition = TRUE)
        }
        clusterPlotVal(clusterRes$plot)
        clusterTabVal(clusterRes$cluster)
        clustNum(input$seleNumCluster)
        if (nrow(data.frame(clusterRes$cluster)) == 0) {
          showModal(modalDialog(
            title = "No cluster found...",
            "Try again with a different setting (e.g., time points included, number of clusters, top % variant).",
            easyClose = TRUE,
            footer = NULL
          ))
        }
      },
      error = function(e) {
        showModal(modalDialog(
          title = "Running Time series clustering failed...",
          "Try again with a different setting (e.g., time points included, number of clusters, top % variant).",
          easyClose = TRUE,
          footer = NULL
        ))})
    })
  })

  # table showing the proteins or phosphopeptides in a selected cluster
  selectedCluster <- reactive({
    if (is.null(clusterTabVal())) {
      NULL
    }
    else {
      inputsValue$seleCluster <- input$seleCluster
      selectedTab <- filter(clusterTabVal(), cluster == input$seleCluster) %>%
        distinct(feature, .keep_all = TRUE) %>%
        mutate(feature = as.character(feature))
      clusterData <- rowData(processedData()[selectedTab$feature,])
      # If the data is phosphoproteomic: include the columns indicating
      # phosphorylation site and peptide sequence
      if ((!is.null(clusterData$site)) && (!is.null(clusterData$Sequence))) {
        selectedTab <- selectedTab %>%
          mutate(Sequence = clusterData$Sequence,
                 site = clusterData$site)
      }
      selectedTab <- selectedTab %>%
        mutate(Gene = clusterData$Gene,
               UniprotID = clusterData$UniprotID,
               prob = formatC(prob, digits=1)) %>%
        select(any_of(c("feature", "Gene", "site", "prob", "cluster", "UniprotID", "Sequence"))) %>%
        arrange(desc(prob)) %>%
        dplyr::rename(ID = feature, probability = prob)
      selectedTab
    }
  })

  # selection of a cluster
  output$seleClusterBox <- renderUI({
    if(!is.null(clusterTabVal())) {
      selectInput("seleCluster", "Select a cluster",
                  sort(unique(clusterTabVal()$cluster)))
    }
  })

  # report number of genes for clustering
  output$numGeneCluster <- renderText({
    sprintf("Number of genes for clustering: %s", nrow(exprMatObj()))
  })

  clustNum <- reactiveVal()

  output$clusterPlotUI <- renderUI({
    if (!is.null(clustNum())) {
      plotOutput("clusterPlot",
                 height = 250*ceiling(as.numeric(clustNum())/3),
                 width = 800)
    }
  })

  output$clusterPlot <- renderPlot({
    clusterPlotVal()
  })

  output$eachClusterTab <- DT::renderDataTable({
    DT::datatable(selectedCluster(), selection = 'single')
  })

  # Plot the time-series of individual proteins
  output$clusterTimePlot <- renderPlot({
    lastClicked <- input$eachClusterTab_row_last_clicked
    if (!is.null(lastClicked)) {
      geneID <- selectedCluster()[lastClicked,]$ID
      if (input$assay == "Phosphoproteome") {
        geneSymbol <- selectedCluster()[lastClicked,]$site
      } else {
        geneSymbol <- selectedCluster()[lastClicked,]$Gene
      }
      inputsValue$geneIDclust <- geneID
      inputsValue$geneSymbolclust <- geneSymbol

      p <- plotTimeSeries(se = processedData(), type = input$clusterFor,
                          geneID = geneID, symbol = geneSymbol,
                          condition = input$seleMetaColTime,
                          treatment = input$seleTreat_cluster,
                          refTreat = input$seleTreat_clusterRef,
                          addZero = input$addZero,
                          zeroTreat = input$seleZeroTreat,
                          timerange = input$seleTimeRange)
      p
    } else {
      NULL
    }
  })

  ############################ Enrichment analysis #############################

  # reactiveVal for file path
  filePath <- reactiveVal(value = NULL)

  # to store the color values used to show gene sets that contain a certain gene
  colorList <- reactiveValues(cols = c())

  # a reactive variable to store GSE result
  GSEres <- reactiveValues(resTab = NULL, resObj = NULL, enrPlot = NULL)

  # a value to check whether enrichment tab is clicked
  clickRecord <- reactiveValues(enrich = FALSE, gene = FALSE, kinase = FALSE)

  # a reactive variable to save gene table of a clicked dot in a
  # cluster enrichment profile
  geneTable <- reactiveVal()

  # if user upload a gene set database
  observeEvent(input$uploadGeneSet, {
    file <- input$uploadGeneSet
    ext <- tools::file_ext(file$datapath)
    req(file)
    validate(need(ext %in% c("gmt"), "Please upload a .gmt file"))
    filePath(file$datapath)
  })
  # if user upload an PTM set database
  observeEvent(input$uploadPTMSet, {
    file <- input$uploadPTMSet
    ext <- tools::file_ext(file$datapath)
    req(file)
    validate(need(ext %in% c("txt"), "Please upload a .txt file"))
    filePath(file$datapath)
  })

  # function to run GSEA
  # Note: the analysisMethod is default to Pathway enrichment if the Proteome
  # assay is selected (might change in the future!)
  resGSEA <- observeEvent(input$RunEnrich, {
    if (input$seleSourceEnrich == "Differential expression" && input$analysisMethod == "Pathway enrichment") {
      # Pathway enrichment for differential expression
      if (!is.null(filterDE())) {
        output$errMsg <- renderText("")
        withProgress(message = "Running enrichment analysis, please wait...", {
          # set color list to empty
          colorList$cols <- NULL
          # reading geneset database
          if(input$seleGeneSet == "select from available gene set databases") {
            inGMT <- loadGSC(paste0("geneset/",input$sigSet),type="gmt")
            inputsValue$sigSet <- input$sigSet
          }
          else {
            inGMT <- loadGSC(filePath(),type="gmt")
          }
          # method
          inputsValue$enrichMethod <- input$enrichMethod
          gseMethod <- input$enrichMethod
          # parameters for GSEA
          if (gseMethod == "GSEA") {
            nPerm <- input$permNum
            inputsValue$permNum <- input$permNum
          }
          # processing differential expression
          corTab <- filterDE() %>%
            arrange(pvalue) %>%
            filter(!duplicated(Gene)) %>%
            arrange(stat)
          # gene level statistics based on user input
          inputsValue$statType <- input$statType
          if(input$statType == "stat") {
            myCoef <- data.frame(row.names = corTab$Gene,
                                 stat = corTab$stat,
                                 stringsAsFactors = FALSE)
          }
          else {
            myCoef <- data.frame(row.names = corTab$Gene,
                                 stat = corTab$log2FC,
                                 stringsAsFactors = FALSE)
          }
          # perform gene set analysis based on selected statistical GSA method
          if (gseMethod == "PAGE") {
            res <- runGSA(geneLevelStats = myCoef,
                          geneSetStat = "page",
                          adjMethod = "fdr",
                          gsc = inGMT,
                          signifMethod = 'nullDist')
          }
          else if (gseMethod == "GSEA") {
            res <- runGSA(geneLevelStats = myCoef,
                          geneSetStat = "gsea",
                          adjMethod = "fdr",
                          gsc = inGMT,
                          signifMethod = 'geneSampling',
                          nPerm = nPerm)
          }

          resTab <- GSAsummaryTable(res)
          colnames(resTab) <- c("Name", "Gene Number", "Stat", "p.up", "p.up.adj",
                                "p.down", "p.down.adj", "Number up", "Number down")
          if(input$ifEnrichFDR) {
            resTab <- filter(resTab,
                             p.up.adj <= input$sigLevel | p.down.adj <= input$sigLevel) %>%
              arrange(desc(Stat))
          } else {
            resTab <- filter(resTab,
                             p.up <= input$sigLevel | p.down <= input$sigLevel) %>%
              arrange(desc(Stat))
          }
          # if all genes are filtered out: show a notification and set
          # GSEres$resTab to NULL to avoid causing error later on
          if (nrow(resTab) > 0) {
            GSEres$resTab <- resTab
          } else {
            GSEres$resTab <- NULL
            showModal(modalDialog(
              title = "No enrichment found...",
              "Try again with a higher p-value cut-off.",
              easyClose = TRUE,
              footer = NULL
            ))}
          GSEres$resTab <- resTab
          GSEres$resObj <- res
          clickRecord$enrich <- FALSE
          clickRecord$gene <- FALSE
        })
      } else {
        output$errMsg <- renderText("Please perform differential expression analysis first or load a previous result!")
      }
    }
    else if (input$seleSourceEnrich =="Selected time-series cluster" && input$analysisMethod =="Pathway enrichment") {
      # Pathway enrichment for time series cluster
      if (nrow(selectedCluster()) > 0) {
        output$errMsg <- renderText("")
        withProgress(message = "Running enrichment analysis, please wait...", {

          # set color list to empty
          colorList$cols <- NULL
          # reading geneset database
          if(input$seleGeneSet == "select from available gene set databases") {
            inGMT <- loadGSC(paste0("geneset/",input$sigSet),type="gmt")
            inputsValue$sigSet <- input$sigSet
          }
          else {
            inGMT <- loadGSC(filePath(),type="gmt")
          }
          # Applying the Fisher's exact test with the runFisher function
          # Other enrichment methods can be added here
          if (input$enrichMethod1 == "Fisher's exact test")
            inputsValue$enrichMethod1 <- input$enrichMethod1
            resTab <- runFisher(unique(selectedCluster()$Gene),
                                reference = unique(rowData(processedData())$Gene),
                                inputSet = inGMT) %>%
            arrange(pval)
          # Filter by the p-value threshold (input$sigLevel)
          if (input$ifEnrichFDR) {
            resTab <- filter(resTab, padj <= input$sigLevel) %>% arrange(padj)
          } else {
            resTab <- filter(resTab, pval <= input$sigLevel) %>% arrange(pval)
          }
          # if all genes are filtered out: show a notification and set
          # GSEres$resTab to NULL to avoid causing error later on
          if (nrow(resTab) > 0) {
            GSEres$resTab <- resTab
          } else {
            GSEres$resTab <- NULL
            showModal(modalDialog(
              title = "No enrichment found...",
              "Try again with a higher p-value cut-off.",
              easyClose = TRUE,
              footer = NULL
            ))}
          GSEres$resObj <- NULL
          clickRecord$enrich <- FALSE
          clickRecord$gene <- FALSE
        })
      } else {
        output$errMsg <- renderText("Please perform time series clustering first!")
      }
    }
    else if (input$seleSourceEnrich == "All time-series clusters" && input$analysisMethod =="Pathway enrichment") {
      if (nrow(clusterTabVal()) > 0) {
        output$errMsg <- renderText("")
        # withProgress("Running enrichment analysis, please wait..", {
        #
        # })
        if(input$seleGeneSet == "select from available gene set databases") {
          database <- loadGSC(paste0("geneset/",input$sigSet),type="gmt")
          inputsValue$sigSet <- input$sigSet
        }
        else {
          database <- loadGSC(filePath(), type="gmt")
        }
        ptm <- FALSE
        clustEnr <- clusterEnrich(clusterTab = clusterTabVal(),
                                  se = processedData(), inputSet = database,
                                  ptm = ptm, filterP = input$sigLevel,
                                  ifFDR = input$ifEnrichFDR)
        if (nrow(clustEnr$table) > 0) {
          GSEres$resTab <- clustEnr$table
          GSEres$enrPlot <- clustEnr$plot
        }
        else {
          GSEres$resTab <- NULL
          showModal(modalDialog(
            title = "No enrichment found...",
            "Try again with a higher p-value cut-off.",
            easyClose = TRUE,
            footer = NULL
          ))}
      }
      else {
        output$errMsg <- renderText("Please perform time series clustering first!")
      }
    }
    else if ((input$seleSourceEnrich == "Differential expression") && (input$analysisMethod == "Phospho-signature enrichment")){
      # Phospho-signature enrichment for differential expression
      if (!is.null(filterDE())) {
        output$errMsg <- renderText("")
        withProgress(message = "Running enrichment analysis, please wait...", {
          # preprocess the input dataframe
          inputTab <- filterDE() %>%
            arrange(pvalue) %>%
            filter(!duplicated(site)) %>%
            arrange(desc(stat))
          # select to use either t-statistic or log2FC
          if (input$statType == "stat") {
            myCoef <- data.frame(row.names = inputTab$site,
                                 stat = inputTab$stat,
                                 stringsAsFactors = FALSE)
          } else if (input$statType == "log2FC") {
            myCoef <- data.frame(row.names = inputTab$site,
                                 stat = inputTab$log2FC,
                                 stringsAsFactors = FALSE)
          }
          # retrieve the phosphodatabase
          if(input$selePTMSet == "select from available PTM set databases") {
            load(paste0("ptmset/", input$sigSetPTM))
            inputsValue$sigSetPTM <- input$sigSetPTM
          }
          else {
            ptmSetDb <- read.table(filePath(), header = TRUE, sep = "\t",
                                   stringsAsFactors = FALSE)
          }

          # perform GSEA
          inputsValue$permNum <- input$permNum
          resTab <- runGSEAforPhospho(geneStat = myCoef, ptmSetDb = ptmSetDb,
                                      nPerm =  input$permNum,
                                      weight = 1, correl.type = "rank",
                                      statistic = "Kolmogorov-Smirnov",
                                      min.overlap = 5) %>%
            as.data.frame()
          colnames(resTab) <- c("Name", "Site.number", "Stat", "Number.pSite.Db", "Number.PTM.site.Db",
                                "pvalue", "Number.up", "Number.down", "padj")
          resTab <- resTab %>%
            select(Name,Site.number,Stat,Number.up,Number.down,Number.pSite.Db,
                   Number.PTM.site.Db, pvalue, padj) # rearrange column order
          if (input$ifEnrichFDR) {
            resTab <- filter(resTab, padj <= input$sigLevel) %>%
              arrange(desc(Stat))
          } else {
            resTab <- filter(resTab, pvalue <= input$sigLevel) %>%
              arrange(desc(Stat))
          }
          # if all genes are filtered out: show a notification and set
          # GSEres$resTab to NULL to avoid causing error later on
          if (nrow(resTab) > 0) {
            GSEres$resTab <- resTab
          } else {
            GSEres$resTab <- NULL
            showModal(modalDialog(
              title = "No enrichment found...",
              "Try again with a higher p-value cut-off.",
              easyClose = TRUE,
              footer = NULL
            ))}
          GSEres$resTab <- resTab
          #GSEres$resObj <- res
          clickRecord$enrich <- FALSE
          clickRecord$gene <- FALSE
        })
      } else output$errMsg <- renderText("Please perform differential expression analysis first and make sure you have selected the Phosphoproteome assay!")
    }
    else if ((input$seleSourceEnrich == "Selected time-series cluster") && (input$analysisMethod == "Phospho-signature enrichment")) {
      # Phospho-signature enrichment for Time-series clustering
      if (nrow(selectedCluster()) > 0) {
        output$errMsg <- renderText("")
        withProgress(message = "Running enrichment analysis, please wait...", {
          # set color list to empty
          colorList$cols <- NULL
          # get PTM set database
          if(input$selePTMSet == "select from available PTM set databases") {
            load(paste0("ptmset/", input$sigSetPTM))
            inputsValue$sigSetPTM <- input$sigSetPTM
          }
          else {
            ptmSetDb <- read.table(filePath(), header = TRUE, sep = "\t",
                                   stringsAsFactors = FALSE)
          }
          # Run the Fisher's exact test for sites in the cluster
          resTab <- runFisher(genes = selectedCluster()$site,
                              reference = rowData(processedData())$site,
                              inputSet = ptmSetDb,
                              ptm = TRUE) %>%
            rename(Site.number = "Gene.number")
          # Filter by the p-value threshold (input$sigLevel)
          if (input$ifEnrichFDR) {
            resTab <- filter(resTab, padj <= input$sigLevel) %>% arrange(padj)
          } else {
            resTab <- filter(resTab, pval <= input$sigLevel) %>% arrange(pval)
          }
          # if all genes are filtered out: show a notification and set
          # GSEres$resTab to NULL to avoid causing error later on
          if (nrow(resTab) > 0) {
            GSEres$resTab <- resTab
          } else {
            GSEres$resTab <- NULL
            showModal(modalDialog(
              title = "No enrichment found...",
              "Try again with a higher p-value cut-off.",
              easyClose = TRUE,
              footer = NULL
            ))}
          GSEres$resObj <- NULL
          clickRecord$enrich <- FALSE
          clickRecord$gene <- FALSE
        })
      } else output$errMsg <- renderText("Please perform time series clustering first and make sure you have selected the Phopshoproteome assay!")
    }
    else if ((input$seleSourceEnrich =="All time-series clusters") && (input$analysisMethod == "Phospho-signature enrichment")) {
      if (nrow(clusterTabVal()) > 0) {
        output$errMsg <- renderText("")
        # withProgress("Running enrichment analysis, please wait..", {
        #
        # })
        if(input$selePTMSet == "select from available PTM set databases") {
          load(paste0("ptmset/", input$sigSetPTM))
          inputsValue$sigSetPTM <- input$sigSetPTM
        }
        else {
            ptmSetDb <- read.table(filePath(), header = TRUE,
                                 sep = "\t",stringsAsFactors = FALSE)
        }
        ptm <- TRUE
        clustEnr <- clusterEnrich(clusterTab = clusterTabVal(),
                                  se = processedData(), inputSet = ptmSetDb,
                                  ptm = ptm, filterP = input$sigLevel,
                                  ifFDR = input$ifEnrichFDR)
        if (nrow(clustEnr$table) > 0) {
          GSEres$resTab <- clustEnr$table
          GSEres$enrPlot <- clustEnr$plot
        }
        else {
          GSEres$resTab <- NULL
          showModal(modalDialog(
            title = "No enrichment found...",
            "Try again with a higher p-value cut-off.",
            easyClose = TRUE,
            footer = NULL
          ))}
      }
      else {
        output$errMsg <- renderText("Please perform time series clustering first!")
      }
    }
  })

  # show differential expressed genes or enrichment results on table 1 (up table)
  output$enrichTab <- DT::renderDataTable({
    if( !is.null (GSEres$resTab)) {
      inputsValue$seleSourceEnrich <- input$seleSourceEnrich
      inputsValue$analysisMethod <- input$analysisMethod
      inputsValue$ifEnrichFDR <- input$ifEnrichFDR
      inputsValue$sigLevel <- input$sigLevel
      resTab <- GSEres$resTab

      if (input$seleSourceEnrich == "Differential expression") {
        if (!is.null(colorList$cols)) {
          # when the bottom table is clicked, color the pathways according to
          # whether it contains the clicked gene
          datatable(resTab,selection = 'single', caption = "Enriched gene/PTM signature sets") %>%
            formatStyle('Stat',background = styleInterval(c(0), c("lightblue", "pink"))) %>%
            formatStyle('Name', color = styleEqual(resTab$Name, colorList$cols)) %>%
            formatRound(which(colnames(resTab) %in% c("Stat", "pvalue", "padj", "p.up","p.up.adj","p.down","p.down","p.down.adj")), digits=3)
        } else {
          # when bottom table was not clicked, do not show colors
          datatable(resTab,selection = 'single', caption = "Enriched gene/PTM signature sets") %>%
            formatStyle('Stat', background = styleInterval(c(0), c("lightblue", "pink"))) %>%
            formatRound(which(colnames(resTab) %in% c("Stat", "pvalue", "padj", "p.up","p.up.adj","p.down","p.down","p.down.adj")), digits = 3)
        }
      } else {
        if (!is.null(colorList$cols)) {
          # when the bottom table is clicked, color the pathways according to whether it contains the clicked gene
          datatable(resTab,selection = 'single', caption = "Enriched gene/PTM signature sets") %>%
            formatStyle('Name', color = styleEqual(resTab$Name, colorList$cols)) %>%
            formatRound(c('pval', 'padj'), digits=3)
        } else {
          # when bottom table was not clicked, do not show colors
          datatable(resTab,selection = 'single', caption = "Enriched gene/PTM signature sets") %>%
            formatRound(c('pval', 'padj'), digits = 3)
        }
      }
    }
  })

  # find gene sets that contain a certain gene
  setGene <- reactive({
    resTab <- GSEres$resTab

    if (input$analysisMethod == "Pathway enrichment" | input$assay == "Proteome") {
      if(input$seleGeneSet == "select from available gene set databases") {
        setList <- loadGSC(paste0("geneset/",input$sigSet),type="gmt")$gsc[resTab$Name]
      }
      else {
        setList <- loadGSC(filePath(),type="gmt")$gsc[resTab$Name]
      }
    }
    else {
      if(input$selePTMSet == "select from available PTM set databases") {
        load(paste0("ptmset/", input$sigSetPTM))
        setList <- ptmSetDb
      }
      else {
        setList <- read.table(filePath(), sep = "\t", header = TRUE,
                              stringsAsFactors = FALSE)
      }
      if (input$seleSourceEnrich == "Selected time-series cluster")
        setList <- setList %>% mutate(signature = ifelse(site.direction == "u",
                                                         paste0(signature,"_upregulated"), paste0(signature, "_downregulated")))
      setList <- setList %>%
        filter(signature %in% resTab$Name, site.ptm == "p") %>%
        separate(site.annotation, sep=":", into = c("site", "PubMedID"),
                 extra = "merge", fill="right")
    }

    if (input$seleSourceEnrich == "Differential expression") {
      if (input$analysisMethod == "Pathway enrichment" | input$assay == "Proteome")
        genes <- filterDE()$Gene else
          sites <- filterDE()$site
    } else { # aka if input$seleSourceEnrich == "Time series clustering"
      if (input$analysisMethod == "Pathway enrichment" | input$assay == "Proteome")
        genes <- unique(selectedCluster()$Gene) else
          sites <- selectedCluster()$site
    }

    if (input$analysisMethod == "Pathway enrichment" | input$assay == "Proteome")
      allSets <- sapply(genes, function(geneName) names(setList)[sapply(setList, function(x) geneName %in% x)]) else {
        allSets <- list()
        for (site in sites) {
          if (site %in% setList$site)
            allSets[[site]] <- setList[setList$site == site,"signature"]}
      }
    allSets <- allSets[sapply(allSets, function(x) length(x) != 0)]
    allSets
  })

  # list genes that enriched in a certain gene set
  gseaList <- reactive({
    if(!is.null(GSEres$resTab)) {
      setName <- GSEres$resTab[as.integer(input$enrichTab_row_last_clicked), "Name"]

      if (input$assay == "Proteome" | input$analysisMethod == "Pathway enrichment") {
        if(input$seleGeneSet == "select from available gene set databases") {
          geneList <- loadGSC(paste0("geneset/",input$sigSet),type="gmt")$gsc[[setName]]
        }
        else {
          geneList <- loadGSC(filePath(),type="gmt")$gsc[[setName]]
        }
      }

      else {
        if(input$selePTMSet == "select from available PTM set databases") {
          load(paste0("ptmset/", input$sigSetPTM))
          geneList <- ptmSetDb
        }
        else {
          geneList <- read.table(filePath(), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        }
        if (input$seleSourceEnrich == "Selected time-series cluster")
          geneList <- geneList %>% mutate(signature = ifelse(site.direction == "u", paste0(signature,"_upregulated"), paste0(signature, "_downregulated")))
        geneList <- geneList %>%
          filter(signature == setName, site.ptm == "p") %>%
          separate(site.annotation, sep=":", into = c("site", "PubMedID"), extra = "merge", fill="right")
      }

      if (input$seleSourceEnrich == "Differential expression") {
        # Differential expression
        corGene <- filterDE()
        if (input$analysisMethod == "Pathway enrichment" | input$assay == "Proteome") {
          geneTab <- corGene[corGene$Gene %in% geneList,]
          # setNum is the number of gene sets containing the gene in question
          geneTab$setNum <- sapply(geneTab$Gene, function(x) length(setGene()[[x]]))
        }
        else {
          geneTab <- corGene[corGene$site %in% geneList$site,]
          geneTab <- merge(x = geneTab, y = geneList[,c("site","site.direction","PubMedID")], by = "site", all.x = TRUE)
          # only select sites whose direction of regulation (up- or down-regulated) matches with the database
          geneTab <- geneTab[((geneTab$log2FC >=0) & (geneTab$site.direction == "u")) | ((geneTab$log2FC < 0) & (geneTab$site.direction ==  "d")),]
          geneTab$setNum <- sapply(geneTab$site, function(x) length(setGene()[[x]]))
        }
        geneTab <- select(geneTab, any_of(c("site", "Gene","log2FC", "pvalue", "padj", "setNum", "Sequence", "PubMedID", "ID")))
      }
      else {
        # time series cluster
        corGene <- selectedCluster()
        if (input$analysisMethod == "Pathway enrichment" | input$assay == "Proteome") {
          geneTab <- corGene[corGene$Gene %in% geneList,]
          # setNum is the number of gene sets containing the gene in question
          geneTab$setNum <- sapply(geneTab$Gene, function(x) length(setGene()[[x]]))
        } else {
          geneTab <- corGene[corGene$site %in% geneList$site,]
          geneTab <- merge(x = geneTab, y = geneList[,c("site","PubMedID")], by = "site", all.x = TRUE)
          geneTab$setNum <- sapply(geneTab$site, function(x) length(setGene()[[x]]))
        }
        geneTab <- select(geneTab, any_of(c("site", "Gene", "cluster", "probability", "setNum", "Sequence", "PubMedID", "ID")))
      }
      geneTab <- geneTab %>% mutate_if(is.numeric, formatC, digits = 2)
      geneTab
    }
  })

  # if the enrichtab is clicked, cancel the color
  observeEvent(input$enrichTab_row_last_clicked, {
    colorList$cols <- NULL
    clickRecord$enrich <- TRUE
  })

  # get the gene ID as well as set color of the gene sets when click the bottom table
  observeEvent(input$geneTab_row_last_clicked, {
    if (input$analysisMethod == "Pathway enrichment" | input$assay == "Proteome")
      clickSym <- gseaList()[as.integer(input$geneTab_row_last_clicked),"Gene"][[1]] else
        clickSym <- gseaList()[as.integer(input$geneTab_row_last_clicked), "site"][[1]]
    colorList$cols <- sapply(GSEres$resTab$Name, function(x) ifelse(x %in% setGene()[[clickSym]], "red", "black"))
    clickRecord$gene <- TRUE
  })

  # using the bottom right table to show genes enriched in the selected set
  output$geneTab <- DT::renderDataTable({
    if (clickRecord$enrich) {
      datatable(select(gseaList(), -ID),
                selection = 'single', rownames = FALSE,
                caption = "Genes/Phosphosites in the selected set")
    }
  })

  output$plotEnr <- renderPlot({
    if (clickRecord$gene) {
      lastClicked <- input$geneTab_row_last_clicked
      geneID <- gseaList()[lastClicked,]$ID
      inputsValue$geneIDenrich <- geneID
      if (input$assay == "Phosphoproteome")
        geneSymbol <- gseaList()[lastClicked,]$site else
          geneSymbol <- gseaList()[lastClicked,]$Gene
      inputsValue$geneSymbolenrich <- geneSymbol

      if (!is.null(lastClicked)) {

        if (input$seleSourceEnrich == "Differential expression") {
          p <- plotBox(se = processedDataSub(), id = geneID, symbol = geneSymbol)
          p
        }
        else if (input$seleSourceEnrich == "Selected time-series cluster") {
          p <- plotTimeSeries(se = processedData(), type = input$clusterFor,
                              geneID = geneID, symbol = geneSymbol,
                              condition = input$seleMetaColTime,
                              treatment = input$seleTreat_cluster,
                              refTreat = input$seleTreat_clusterRef,
                              addZero = input$addZero,
                              zeroTreat = input$seleZeroTreat,
                              timerange = input$seleTimeRange)
          p
        }
      } else { NULL }
    } else { NULL }
  })

  # render the cluster enrichment profile plot as interactively plotly plot
  output$clustEnrPlot <- renderPlotly({
    if(!is.null(GSEres$resTab)) {
      p <- ggplotly(GSEres$enrPlot, source = "enrich",
                    height = 60*length(unique(GSEres$resTab$Name)))
      p %>% event_register("plotly_click")
      #ggplotly(GSEres$enrPlot)
      p
    }
  })

  # if a dot on the plot is clicked
  observeEvent(event_data("plotly_click", source = "enrich"),{
    d <- event_data("plotly_click", source = "enrich")
    tab <- GSEres$resTab
    name <- d$key
    clusterSele <- d$customdata
    genes <- tab[tab$Name == name & tab$cluster == clusterSele, ]$Genes

    selectedTab <- filter(clusterTabVal(), cluster == clusterSele) %>%
      distinct(feature, .keep_all = TRUE) %>%
      mutate(feature = as.character(feature))
    clusterData <- rowData(processedData()[selectedTab$feature,])

    # If the data is phosphoproteomic: include the columns indicating
    # phosphorylation site and peptide sequence
    if ((!is.null(clusterData$site)) && (!is.null(clusterData$Sequence))) {
      selectedTab <- selectedTab %>%
        mutate(Sequence = clusterData$Sequence,
               site = clusterData$site)
    }
    selectedTab <- selectedTab %>%
      mutate(Gene = clusterData$Gene,
             UniprotID = clusterData$UniprotID,
             prob = formatC(prob, digits=1)) %>%
      select(any_of(c("feature", "Gene", "site", "prob", "cluster", "UniprotID", "Sequence"))) %>%
      arrange(desc(prob)) %>%
      dplyr::rename(ID = feature, probability = prob)
    if (input$analysisMethod == "Pathway enrichment") {
      geneTable <- selectedTab[selectedTab$Gene %in% genes[[1]],]
    }
    else {
      geneTable <- selectedTab[selectedTab$site %in% genes[[1]],]
    }

    geneTable(geneTable)

    output$geneTabClicked <- DT::renderDataTable({
      datatable(geneTable, selection = "single")
    })
  })

  # when a row from gene table is clicked
  observeEvent(input$geneTabClicked_row_last_clicked,{
    lastClicked <- input$geneTabClicked_row_last_clicked
    geneID <- geneTable()[lastClicked,]$ID
    inputsValue$geneIDenrich <- geneID
    if (input$assay == "Phosphoproteome")
      geneSymbol <- geneTable()[lastClicked,]$site else
        geneSymbol <- geneTable()[lastClicked,]$Gene
    inputsValue$geneSymbolenrich <- geneSymbol

    p <- plotTimeSeries(se = processedData(), type = input$clusterFor,
                        geneID = geneID, symbol = geneSymbol,
                        condition = input$seleMetaColTime,
                        treatment = input$seleTreat_cluster,
                        refTreat = input$seleTreat_clusterRef,
                        addZero = input$addZero,
                        zeroTreat = input$seleZeroTreat,
                        timerange = input$seleTimeRange)
    output$plotGeneEnr <- renderPlot(p)
  })

  # link for download the enriched set list
  output$downloadUI1 <- renderUI({
    if(!is.null(GSEres$resTab)) {
      downloadLink("downloadSet", "Download enriched set list")
    }
  })

  output$downloadUI2 <- renderUI({
    if (clickRecord$enrich) {
      downloadLink("downloadGene", "Download gene/phosphosite list")
    }
  })

  # a link to download enrichment results as csv
  output$downloadSet <- downloadHandler(
    filename = function() { paste('geneSetList', '.tsv', sep='') },
    content = function(file) {
      write.csv2(GSEres$resTab, file)
    }
  )

  # a link to download genes in the clicked set as csv
  output$downloadGene <- downloadHandler(
    filename = function() { paste('geneList', '.tsv', sep='') },
    content = function(file) {
      write.csv2(gseaList(), file)
    }
  )

  ###### Not in use: Plot the enrichment map for enrichment analysis ###########

  output$plot4 <- renderPlot({
    if (!is.null(GSEres$resObj)) {
      output$errMsg2 <- renderText("")
      networkPlot(GSEres$resObj, class="distinct", direction = input$setDirection,
                  significance = input$pCut, adjusted = input$ifAdjCut, overlap = input$setOverlap,
                  lay = as.numeric(input$layMethod), label = input$labelMethod, cexLabel = input$sizeLable,
                  ncharLabel = 100)

    } else output$errMsg2 <- renderText("Please run enrichment analysis first! (Enrichment on time series clusters not supported yet)")
  })

  ########################## Kinase activity inference #########################

  ####Widgets

  # reactive variable to store the decoupler network
  # (update itself if change organism)
  decoupler_network <- reactive({
    network <- getDecouplerNetwork(input$speciesRef)
    inputsValue$speciesRef <- input$speciesRef
    network
  })
  # reactive variable to store the kinase activity result
  kinaseRes <- reactiveVal()

  # Perform the analysis when click the runKinase button
  runKinaseAnalysis <- observeEvent(input$runKinase, {
    if (input$assay == "Phosphoproteome") {
      if (input$seleSourceKinase == "Differential expression") {
        inputsValue$seleSourceKinase <- input$seleSourceKinase
        if (!is.null(filterDE())) {
          withProgress(message = "Running kinase activity inference, please wait...", {
            output$errMsgKinase <- renderText("")
            # compute the kinase score or report error
            # (usually because the selected organism was not correct)
            tryCatch({
              scoreTab <- calcKinaseScore(filterDE(), decoupler_network(),
                                          statType = input$statTypeKinase,
                                          nPerm = input$nPermKinase)
              inputsValue$statTypeKinase <- input$statTypeKinase
              inputsValue$nPermKinase <- input$nPermKinase
              scoreTab <- scoreTab %>% mutate(padj = p.adjust(p_value, method = "BH")) %>% arrange(p_value)
              kinaseRes(scoreTab)
            }, error = function(e) {
              showModal(modalDialog(
                title = "No kinase found...",
                "Try again with a different organism or a larger phosphosite list.",
                easyClose = TRUE,
                footer = NULL
              ))})
          })
        } else output$errMsgKinase <- renderText("Please perform Differential expression analysis first")
      } else {

        # Perform the analysis for time-series cluster
        if (!is.null(clusterTabVal())) {
          withProgress(message = "Running kinase activity inference, please wait...",  {
            output$errMsgKinase <- renderText("")
            # Add a column showing phosphosites based on the ID
            clusterData <- clusterTabVal()
            clusterData <- clusterData[clusterData$cluster == input$seleCluster,]
            inputsValue$seleClusterKinase <- input$seleCluster
            allClusterFeature <- clusterData %>%
                distinct(feature, .keep_all = TRUE) %>% .$feature
            allClusterSite <- data.frame(rowData(processedData())[allClusterFeature, "site"])
            allClusterSite$feature <- allClusterFeature
            clusterData <- clusterData %>%
              left_join(allClusterSite, by = "feature") %>%
              rename(site = "rowData.processedData....allClusterFeature...site..")
            if (input$seleKinaseTimeMethod == "activity") {
              inputsValue$seleKinaseTimeMethod <- input$seleKinaseTimeMethod
              # compute the kinase ACTIVITY score
              # initiate an empty dataframe to store score result
              scoreTab <- data.frame(source = c(), score = c(),
                                     p_value = c(), timepoint = c())
              # get the order of time points
              timeVector <- input$seleTimeRange
              inputsValue$seleTimeRangeKinase <- input$seleTimeRange
              timeUnit <- suppressWarnings(str_extract(timeVector, "h|min"))
              timeUnit <- ifelse(is.na(timeUnit), "", timeUnit)
              # If both h and min are present, divide the min time points by 60
              if ((any(timeUnit == "h")) && (any(timeUnit == "min"))) {
                timeValue <- timeVector
                timeValue[timeUnit == "min"] <- 1/60 * as.numeric(gsub("min", "", timeValue[timeUnit == "min"]))
                timeRank <- rank(as.numeric(gsub("h", "", timeValue)))
              } else {
                timeRank <- rank(as.numeric(gsub("h|min", "", timeVector)))
              }
              timeRank <- as.character(timeRank)
              # try computing the kinase score. If error occurs (usually due to wrong organism selected) then inform the user to try again
              tryCatch({
                if (input$clusterFor == "expression") {
                  # Handling for expression case
                  siteTab <- processedData()[rowData(processedData())$site %in% unique(clusterData$site),
                                             processedData()$treatment == input$seleTreat_cluster & processedData()$timepoint %in% input$seleTimeRange]
                  # an empty vector to store order of timepoints in the heatmap
                  timeOrder = c()
                  # Compute the fold change of a time point with respect to the previous time point
                  for (i in 2:length(timeRank)) {
                    time2 <- timeVector[timeRank == as.character(i)]
                    time1 <- timeVector[timeRank == as.character(i-1)]
                    siteTime1 <- siteTab[,siteTab$timepoint == time1]
                    siteTime2 <- siteTab[,siteTab$timepoint == time2]
                    # match samples by subjectID if provided
                    if (!is.null(siteTab$subjectID)) {
                      siteTime2 <- siteTime2[,match(siteTime1$subjectID,siteTime2$subjectID)]
                      fc <- assay(siteTime2) - assay(siteTime1)
                      clusterDataTime <- data.frame(site = as.character(rowData(siteTab)$site),
                                                    log2FC = rowMeans(fc))
                    } else {
                      fc <- rowMeans(assay(siteTime2)) - rowMeans(assay(siteTime1))
                      clusterDataTime <- data.frame(site = as.character(rowData(siteTab)$site),
                                                    log2FC = fc)
                    }
                    clusterDataTime <- na.omit(clusterDataTime)
                    scoreTabTime <- calcKinaseScore(clusterDataTime, decoupler_network(),statType = "log2FC", nPerm = input$nPermKinase)
                    inputsValue$nPermKinase <- input$nPermKinase
                    inputsValue$statTypeKinase <- input$statTypeKinase
                    scoreTabTime$timepoint <- paste0(time2,"_",time1)
                    timeOrder <- append(timeOrder, paste0(time2,"_",time1))
                    scoreTab <- rbind(scoreTab, scoreTabTime)
                  }
                  # make sure the order of timepoints on the heatmap will be correct
                  scoreTab$timepoint <- factor(scoreTab$timepoint, levels = timeOrder)
                } else if (input$clusterFor == "logFC") {
                  # Handling for logFC case
                  for (time in input$seleTimeRange) {
                    clusterDataTime <- clusterData[clusterData$time == time,] %>%
                      select(value, site) %>% rename(log2FC = "value")
                    scoreTabTime <- calcKinaseScore(clusterDataTime, decoupler_network(), statType = "log2FC", nPerm = input$nPermKinase)
                    inputsValue$statTypeKinase <- input$statTypeKinase
                    inputsValue$nPermKinase <- input$nPermKinase
                    scoreTabTime$timepoint <- time
                    scoreTab <- rbind(scoreTab, scoreTabTime)
                  }
                  # make sure the time points are in correct order in the heatmap
                  scoreTab$timepoint <- factor(scoreTab$timepoint, levels = timeVector[order(match(timeRank, sort(timeRank)))])
                } else if (input$clusterFor == "two-condition expression") {
                  # Handling for two-condition expression case
                  siteTab <- processedData()[rowData(processedData())$site %in% unique(clusterData$site),
                                             processedData()$treatment == input$seleTreat_cluster & processedData()$timepoint %in% input$seleTimeRange]
                  refTab <- processedData()[rowData(processedData())$site %in% unique(clusterData$site),
                                            processedData()$treatment == input$seleTreat_clusterRef & processedData()$timepoint %in% input$seleTimeRange]
                  # compute the logFC for each time point
                  for (time in input$seleTimeRange) {
                    # if subjectID is present then match samples by subject ID
                    siteTabTime <- siteTab[,siteTab$timepoint == time]
                    refTabTime <- refTab[,refTab$timepoint == time]
                    if ((!is.null(siteTab$subjectID)) && (!is.null(refTab$subjectID))) {
                      siteTabTime <- siteTabTime[,match(refTabTime$subjectID, siteTabTime$subjectID)]
                      fc <- assay(siteTabTime) - assay(refTabTime)
                      clusterDataTime <- data.frame(site = as.character(rowData(siteTabTime)$site),
                                                    log2FC = rowMeans(fc))
                    } else{
                      fc <- rowMeans(assay(siteTabTime)) - rowMeans(assay(refTabTime))
                      clusterDataTime <- data.frame(site = as.character(rowData(siteTabTime)$site),
                                                    log2FC = fc)
                    }
                    clusterDataTime <- na.omit(clusterDataTime)
                    scoreTabTime <- calcKinaseScore(clusterDataTime, decoupler_network(), statType = "log2FC", nPerm = input$nPermKinase)
                    inputsValue$statTypeKinase <- input$statTypeKinase
                    inputsValue$nPermKinase <- input$nPermKinase
                    scoreTabTime$timepoint <- time
                    scoreTab <- rbind(scoreTab, scoreTabTime)
                  }
                  # make sure the time points are in correct order in the heatmap
                  scoreTab$timepoint <- factor(scoreTab$timepoint, levels = timeVector[order(match(timeRank, sort(timeRank)))])
                }
                scoreTab <- scoreTab %>% mutate(padj = p.adjust(p_value, method = "BH")) %>%
                  arrange(p_value)
                kinaseRes(scoreTab)
              }, error = function(e) {
                kinaseRes(NULL)
                showModal(modalDialog(
                  title = "No kinase found...",
                  "Try again with a different organism or a larger phosphosite list.",
                  easyClose = TRUE,
                  footer = NULL
                ))})
            } else {
              inputsValue$seleKinaseTimeMethod <- input$seleKinaseTimeMethod
              # compute how likely the kinases are associated with the cluster
              # An error is induced if no kinase is found (wrong organism or list too small)
              tryCatch({
                inputsValue$seleAssoMethod <- input$seleAssoMethod
                if (input$seleAssoMethod == "Fisher's exact test") {
                  # Fisher's exact test to test kinase association with cluster
                  pSiteCluster  <- unique(selectedCluster()$site)
                  pSiteRefList <- unique(rowData(processedData())$site)
                  pSiteRefList <- pSiteRefList[!pSiteRefList %in% pSiteCluster]
                  rtab <- lapply(unique(decoupler_network()$source), function(kinase) {
                    kinaseTarget = as.character(decoupler_network()[decoupler_network()$source == kinase,"target"])
                    RinSet = sum(pSiteRefList %in% kinaseTarget)
                    RninSet = length(pSiteRefList) - RinSet
                    GinSet = sum(pSiteCluster %in% kinaseTarget)
                    GninSet = length(pSiteCluster) - GinSet
                    fmat = matrix(c(GinSet, RinSet, GninSet, RninSet), nrow = 2,
                                  ncol = 2, byrow = FALSE)
                    colnames(fmat) = c("inSet", "ninSet")
                    rownames(fmat) = c("genes", "reference")
                    fish = fisher.test(fmat, alternative = "greater")
                    pval = fish$p.value
                    inSet = RinSet + GinSet
                    tibble(source = kinase,
                           `number.pSite.in.cluster`= GinSet,
                           `number.pSite.by.kinase` = inSet,
                           p_value = pval)
                  }) %>% bind_rows() %>%
                    filter(number.pSite.in.cluster>0)
                } else { # i.e. if seleAssoMethod == 'FGSEA'
                  # GSEA to test kinase association with cluster using FGSEA
                  # from decoupleR. Data is taken from clusterTabVal() instead
                  # of selectedCluster() since the former has unrounded
                  # probability values
                  selectedTab <- filter(clusterTabVal(), cluster == input$seleCluster) %>%
                    distinct(feature, .keep_all = TRUE) %>%
                    mutate(feature = as.character(feature))
                  clusterData <- rowData(processedData()[selectedTab$feature,])
                  inputTab <- selectedTab %>%
                    mutate(site = clusterData$site) %>%
                    select(site, prob) %>%
                    column_to_rownames(var = "site") %>%
                    arrange(desc(prob))
                  rtab <- decoupleR::run_fgsea(mat = inputTab,
                                               network = decoupler_network(),
                                               minsize = 1,
                                               times = input$nPermKinase) %>%
                    filter(statistic == "fgsea") %>%
                    select(-condition, -statistic) %>%
                    rename(enrich.score = "score")
                  inputsValue$nPermKinase <- input$nPermKinase
                }
                # adjusting p-values and filtering based on (adjusted) p-values
                rtab <- rtab %>%
                  mutate(padj = p.adjust(p_value, method = "BH")) %>%
                  arrange(p_value)
                if (input$ifKinaseFDR)
                  rtab <- rtab %>% filter(padj <= input$pKinase) else
                    rtab <- rtab %>% filter(p_value <= input$pKinase)
                inputsValue$ifKinaseFDR <- input$ifKinaseFDR
                inputsValue$pKinase <- input$pKinase
                # if no kinase: induce an error to show the pop-up window
                if (nrow(rtab) == 0)
                  stop("No kinase found... perhaps the wrong organism was chosen or the p-value threshold was too low") else
                    kinaseRes(as.data.frame(rtab))
              }, error = function(e) {
                kinaseRes(NULL)
                showModal(modalDialog(
                  title = "No kinase found...",
                  "Try again with a different organism or a larger phosphosite list.",
                  easyClose = TRUE,
                  footer = NULL
                ))})
              }
          })
        } else output$errMsgKinase <- renderText("Please do a time-series clustering (logFC) first")
      }}  else output$errMsgKinase <- renderText("This feature only support phosphoproteome data. Please double check the Preprocessing step!")
  })

  # output to display plot object
  output$plotKinase <- renderPlot({
    if (!is.null(kinaseRes())) {
      scoreTab <- kinaseRes()
      if (input$ifKinaseFDR)
        scoreTab$p_value <- scoreTab$padj
      if (input$seleSourceKinase == "Differential expression") {
        plot <- plotKinaseDE(scoreTab, nTop = input$nTopKinase, pCut = input$pKinase)
        plot
      } else {
        plot <- plotKinaseTimeSeries(scoreTab, pCut = input$pKinase, clusterName = input$seleCluster)
        plot
      }
    }
  })

  # output to display table of kinase score
  output$kinaseTab <- DT::renderDataTable({
    if (!is.null(kinaseRes())) {
      if ((input$seleSourceKinase == "Differential expression") |
          ((input$seleSourceKinase == "Time-series cluster") &
           (input$seleKinaseTimeMethod == "activity"))){
        # table for kinase activity
        resTab <- kinaseRes() %>% rename(Kinase = "source", Activity = "score")
        datatable(resTab, selection = 'single', rownames = FALSE,
                  caption = "Kinase activity") %>%
          formatStyle('Activity',background=styleInterval(c(0),c("lightblue","pink"))) %>%
          formatRound(c("Activity", "p_value", "padj"), digits = 3)
      } else { # i.e. if doing kinase association for time-series cluster
        # table for kinase association
        resTab <- kinaseRes() %>% rename(Kinase = "source")
        if (input$seleAssoMethod == "Fisher's exact test") {
        datatable(resTab, selection = 'single', rownames = FALSE,
                  caption = paste0("Kinases associated to ", input$seleCluster, ", Fisher's exact test")) %>%
          formatRound(c("p_value", "padj"), digits = 3)
        } else { # i.e. if use FGSEA to estimate kinase association
          datatable(resTab, selection = 'single', rownames = FALSE,
                    caption = paste0("Kinases associated to ", input$seleCluster, ", FGSEA")) %>%
            formatStyle('enrich.score',background=styleInterval(c(0),c("lightblue","pink"))) %>%
            formatRound(c("enrich.score","p_value", "padj"), digits = 3)
        }
      }
    }
  })

  # list of phosphosites targeted by the selected kinase
  pSiteList <- reactive({
    kinase <- as.character(kinaseRes()[as.integer(input$kinaseTab_row_last_clicked), "source"])
    pSite <- as.character(decoupler_network()[decoupler_network()$source == kinase, "target"])
    inputsValue$kinaseLastClicked <- kinase
    if (input$seleSourceKinase == "Differential expression")
      pSiteTab <- filterDE()[filterDE()$site %in% pSite,] else
        pSiteTab <- selectedCluster()[selectedCluster()$site %in% pSite,]
    pSiteTab <- pSiteTab %>% select(any_of(c("Gene","site","log2FC","stat","pvalue","padj","cluster","probability","Sequence")))
    pSiteTab
  })
  # get phosphosites that are targets of the selected kinase in kinaseTab
  observeEvent(input$kinaseTab_row_last_clicked, {
    clickRecord$kinase <- TRUE
  })

  # table to show level of phosphosite corresponding to the selected kinase
  output$pSiteTab <- DT::renderDataTable({
    if (clickRecord$kinase) {
      tryCatch({
        datatable(pSiteList(),
                  selection = 'single', rownames = FALSE,
                  caption = "Phosphosites targeted by the selected kinase") %>%
          formatRound(c("log2FC", "stat", "pvalue",  "padj"), digits = 3)
      }, error = function(e) {
        datatable(pSiteList(),
                  selection = 'single', rownames = FALSE,
                  caption = "Phosphosites targeted by the selected kinase")
      })
    }
  })
  #link for download the kinase score result
  output$downloadUI3 <- renderUI({
    if( !is.null (kinaseRes())) {
      downloadLink("downloadKinase","Download kinase activity inference result")
    }
  })
  #a link to download enrichment results as csv
  output$downloadKinase <- downloadHandler(
    filename = function() { paste('kinaseActivity', '.tsv', sep='') },
    content = function(file) {
      write.csv2(kinaseRes(), file)
    }
  )

  #link for download the phosphosite result
  output$downloadUI4 <- renderUI({
    if( !is.null(pSiteList())) {
      downloadLink("downloadPhosphoSite","Download phosphopeptide list")
    }
  })
  #a link to download phosphopeptide result as tsv
  output$downloadPhosphoSite <- downloadHandler(
    filename = function() { paste('phosphopeptide', '.tsv', sep='') },
    content = function(file) {
      write.csv2(pSiteList(), file)
    }
  )

  ################################# log info ###################################

  infoTable <- read.delim(file = 'infoTable.tsv', sep = '\t')

  allInputs <- reactive({
    values <- reactiveValuesToList(inputsValue)
    df <- lapply(names(values), function(n) {
      v <- paste0(unlist(values[n]), collapse = ",")
      data.frame(name = n, value = v)
    }) %>% bind_rows()
    df <- left_join(df, infoTable, by = "name")
  })



  output$show_inputs <- renderDataTable({
    datatable(allInputs(), filter = "top", rownames = FALSE,
              selection = "none", style = "bootstrap")
  })

  output$downloadLogValues <- downloadHandler(
    filename = function() { paste('Logvalues', '.tsv', sep='') },
    content = function(file) {
      write.table(allInputs(), file = file, quote = FALSE, sep = '\t',
                  col.names = NA)
    })

})
