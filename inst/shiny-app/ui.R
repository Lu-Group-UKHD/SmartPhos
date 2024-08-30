# This script contains the user-interface (UI) part for the Phosphoproteomics Shiny app. 
# The UI is simple and very easy to navigate. It contains multiple tabs on the
# top, each focusing on one aspect of the phosphoproteomics and proteomics data
# exploration. User simply has to upload a MultiAssayExperiment object generated 
# by the SmartPhos package and then they are ready to explore. 

library(shiny)
library(shinythemes)
library(shinyjs)
library(shinyBS)
library(shinyWidgets)
library(plotly)

# The top level navigation UI for the app
navbarPage("SmartPhos explorer",
           theme =  shinytheme("flatly"),
           useShinyjs(),
           
           tags$style(type="text/css",
                      ".shiny-output-error { visibility: hidden; }",
                      ".shiny-output-error:before { visibility: hidden; }"
           ),
           
           # First tab. This tab contains the preprocessing options before diving
           # into the data. This tab allows the user to upload a MultiAssayExperiment
           # object or upload a zip file containing mass spec files by Spectronaut
           # or MaxQuant. It contains options for selecting between the proteome
           # or phosphoproteome assay. Options for choosing the transformation method,
           # normalization, and displaying the missing values. 
           tabPanel("Preprocesing",
                    h5("Updated on 01.07.2024, please check the Update Notes panel for detailed information"),
                    # h4(a("Instruction on how to use SmartPhosExplorer", 
                    #      href="https://www.huber.embl.de/users/jlu/smartPhos/vignette.html",
                    #      target="_blank")),
                    titlePanel("Data input and preprocessing"),
                    sidebarLayout(
                      sidebarPanel(
                        # Input for R object upload or zip file upload
                        bsCollapse(id = "collapse1", open = "Upload",
                                   bsCollapsePanel("Upload",
                                                   "This panel has options for uploading Mass Spectrometry files or an R MultiAssayExperiment object.",
                                                   radioButtons("upload", "Select the upload method:", 
                                                                c("Upload a zip file", "Upload a R object"),
                                                                inline = TRUE),
                                                   conditionalPanel(
                                                     condition = "input.upload == 'Upload a zip file'",
                                                     fileInput("uploadZip", "Upload the zip file containing raw MS files"),
                                                     radioButtons("tool", "Select tool",
                                                                  c("Spectronaut", "MaxQuant"),
                                                                  selected = "Spectronaut",
                                                                  inline = TRUE),
                                                     uiOutput("colAnnoBoxPreprocess"),
                                                     actionButton("processUpload", "Process uploaded data")
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.upload == 'Upload a R object'",
                                                     fileInput("uploadObject", "Upload a MultiAssayExperiment object")
                                                   ),
                                                   uiOutput("option"),
                                                   # UI panel based on the selected assay
                                                   uiOutput("seleAssayBox"),
                                                   # Output number of samples and features
                                                   uiOutput("dataInfo"),
                                                   h4("Selection for multiAssayExperiment object", style = "color:steelblue;font-weight: bold"),
                                                   textInput("text", "Prefix for the file", "Phospho1"),
                                                   uiOutput("saveListBox"),
                                                   # options for saving, loading and removing of objects from the Shiny app
                                                   fluidRow(
                                                     column(downloadButton("download", "Download",
                                                                           style="color: #FFFFFF; background-color: #3498DB; border-color: #2E86C1;"),
                                                            width = 3),
                                                     column(actionButton("save", "Save",
                                                                         style="color: #FFFFFF; background-color: #58D68D; border-color: #2ECC71;"),
                                                            width = 2),
                                                     column(actionButton("load", "Load",
                                                                         style="color: #FFFFFF; background-color: #F5B041; border-color: #F39C12"),
                                                            width = 2),
                                                     column(actionButton("remove", "Remove",
                                                                         style="color: #FFFFFF; background-color: #DC7633; border-color: #D35400"),
                                                            width = 2),
                                                   ),
                                                   style = "info")),
                        # Input regarding transformation method, normalization and imputation
                        bsCollapsePanel("Preprocessing options", "This panel has different preprocessing options for the selected assay.",
                                        uiOutput("seleNormCorrect"),
                                        uiOutput("ifAlreadyNormBox"),
                                        uiOutput("ifNormByProteinBox"),
                                        selectInput("transform", "Select a transformation method",
                                                    c("none","log2","vst"), selected = "log2"),
                                        radioButtons("normalize", "Normalization",
                                                     c("Yes" = TRUE, "No" = FALSE),
                                                     selected = TRUE,
                                                     inline = TRUE),
                                        sliderInput("missFilter", "Percent of missing values allowed",
                                                    min = 0, max = 100, value = 50, step = 5),
                                        bsTooltip("missFilter",
                                                  "Proteins with missing values above this threshold will be removed from all analyses.",
                                                  placement = "bottom",
                                                  trigger = "hover",
                                                  options = NULL),
                                        selectInput("impute", "Select imputation method",
                                                    c("none" = "none",
                                                      "MinDet (fast method)" = "MinDet",
                                                      "QRILC (fast method)" = "QRILC",
                                                      "MLE" = "MLE",
                                                      "bpca" = "bpca",
                                                      "missForest (slow method)" = "missForest"),
                                                    selected = "QRILC"),
                                        checkboxInput("batch", "Remove batch effects", value = FALSE),
                                        uiOutput("seleColBatch"),
                                        actionButton("processSelection", "Process"),
                                        uiOutput("downloadSE"),
                                        style = "warning"),
                        hr(),
                        actionButton("missingV", "Plot completness of the assay",
                                     style="color: #FFFFFF; background-color: #3498DB; border-color: #2E86C1"),
                        hr(),
                        actionButton("launch_app", "Launch MatrixQCvis", 
                                     style="color: #FFFFFF; background-color: #CD5C5C; border-color: #D35400"),
                      ),
                    mainPanel(
                      textAreaInput('outliers', 'Enter sample outliers (comma delimited) to be removed',
                                    width = "100%", rows = 5, resize = "both"),
                      tableOutput("label"),
                      uiOutput("colorBoxUI"),
                      plotOutput("boxPlot"),
                      plotOutput("boxPlotLogRatio"),
                      DT::dataTableOutput("metaData"),
                      plotOutput("missingPlot")
                ))),
           
           # This tab performs principal component analysis (PCA) on the imputed
           # assay from the first tab and then plot the principal components. The
           # user has options for selecting the principal component for the axes.
           # Options for colouring and shaping the data points based on different
           # features are available. Users can also download the plot as PDF file
           # PNG file and can also download the PCA values as a .tsv file.
           tabPanel("PCA",
                    titlePanel("Principal Component Analysis and visualizing the
                               Principal Components"),
                    sidebarLayout(
                      sidebarPanel(
                        actionButton("RunPCA", "Run PCA"),
                        hr(),
                        downloadLink('downloadPCATable', 'Download PCA values'),
                        hr(),
                        # rendering UI for different plotting features
                        uiOutput("xaxisPCAui"),
                        uiOutput("yaxisPCAui"),
                        uiOutput("colorPCAui"),
                        uiOutput("shapePCAui"),
                        hr(),
                        h4("Download PCA plot as pdf", style = "color:steelblue;font-weight: bold"),
                        numericInput("figWidthPCA", label = "PDF figure width", value = 16),
                        numericInput("figHeightPCA", label = "PDF figure height", value = 15),
                        downloadButton("downPCA", "Download")
                      ),
                      mainPanel(
                        span(textOutput("messagePCA"), style = "color:red; font-size:20px"),
                        plotlyOutput("pcplot", height = 800, width = 1000)
                    ))),
           
           # This tabs allow the user to plot and download heatmap for the imputed
           # assay from the first tab. There is option to choose the top variants
           # genes, differentially expressed genes or genes from the selected time
           # series cluster. The user can also divide the heatmap plot into the 
           # row and/or column clusters based on Hierarchical clustering and can
           # add additional annotations. Heatmap can be dowloaded as a PDF too. 
           tabPanel("Heatmap",
                    titlePanel("Visualizing heatmap"),
                    sidebarLayout(
                      sidebarPanel(
                        radioButtons("chooseType",label = "Genes to plot", 
                                     choices = c("Top variant", "Differentially expressed", "Selected time series cluster"),
                                     selected = "Top variant"),
                        conditionalPanel(condition = "input.chooseType =='Top variant'",
                                         uiOutput("topGenes"),
                                         uiOutput("colCluster"),
                                         uiOutput("rowCluster")),
                        # additional annotations for heatmap
                        uiOutput("colAnnoBoxHM"),
                        h4("Download heatmap as pdf", style = "color:steelblue;font-weight: bold"),
                        numericInput("figWidthHM", label = "PDF figure width", value = 16),
                        numericInput("figHeightHM", label = "PDF figure height", value = 15),
                        downloadButton("downHM", "Download")
                      ),
                      mainPanel(
                        actionButton("doPlot","Plot!", align = "center" ),
                        span(textOutput("errMsg1"), style = "color:red"),
                        plotOutput("plotHM", width = 1000, height = 1000)
                      ))),
           
           # This tab performs differential expression analysis on the transformed
           # and normalized assay from the first tab. The two methods available are
           # limma and ProDA. Users have options to filter the differentially expressed
           # genes table based on p-value and log fold change value. Box plot
           # for the comparison can be done by simply clicking on the row. There
           # is Volcano plot also available for visualization. The user can download
           # the DE table as a .tsv file.
           tabPanel("Differential expression",
                    titlePanel("Identify DE genes between treatments or time points"),
                    sidebarLayout(
                      sidebarPanel(
                        h4("DE tests",style = "color:steelblue;font-weight: bold"),
                        checkboxInput("seleID", "Select samples by sample ID", value = FALSE),
                        conditionalPanel(condition = "input.seleID",
                                         uiOutput("ID1Box"),
                                         uiOutput("ID2Box")),
                        conditionalPanel(condition = "!input.seleID",
                                         uiOutput("seleMetaColBoxDiff"),
                                         uiOutput("treat1Box"),
                                         uiOutput("time1Box"),
                                         uiOutput("treat2Box"),
                                         uiOutput("time2Box")),
                        uiOutput("infoDE"),
                        uiOutput("seleMethodBox"),
                        actionButton("runDE",label = "Run DE analysis"),
                        textInput("pFilter","p value cutoff", value = 0.05),
                        checkboxInput("ifAdjusted","Use adjusted p values", value = FALSE),
                        textInput("fcFilter","log2FC cutoff", value = 0.5)
                      ),
                      mainPanel(
                        DT::dataTableOutput("DEtab"),
                        uiOutput("downloadTableUI"),
                        fluidRow(
                          splitLayout(cellWidths = c("60%", "40%"), uiOutput("ui.plot"), plotlyOutput("plotVolcano"))
                        )
                      ))),
           
           # This tab performs clustering of time series gene expression patterns
           # using fuzzy c-mean clustering. There is option to perform hypothesis 
           # testing (spline fitting) to filter out proteins whose patterns are
           # not consistent. Samples at different time points are linked by patient
           # ID if provided (PatID in fileTable). Currently, time points should 
           # have either no unit or "h" and/or "min" unit. The time points would
           # have no unit in a future installment. Zero timepoint can also be added 
           # when treatments with no zero timepoints are selected.
           tabPanel("Time series clustering",
                    titlePanel("Clustering of time series gene expression patterns"),
                    sidebarLayout(
                      sidebarPanel(
                        uiOutput("seleMetaColBoxTime"),
                        radioButtons("clusterFor", "Perform clustering on:", c("expression", "logFC", "two-condition expression"), inline = TRUE),
                        conditionalPanel(condition = "input.clusterFor == 'logFC' | input.clusterFor == 'two-condition expression'", 
                                         uiOutput("clusterTreatRefBox")),
                        uiOutput("clusterTreatBox"),
                        uiOutput("timerangeBox"),
                        uiOutput("zeroTime"),
                        uiOutput("zeroTreatInfo"),
                        uiOutput("zeroTreat"),
                        sliderInput("topVarTime", "Use top % variant genes along time" , 0, 100, 80, 5),
                        checkboxInput("ifFilterFit","Filter genes based on spline fit test", value = FALSE),
                        conditionalPanel(condition = "input.ifFilterFit == true", 
                                         textInput("pSpline", "P value cut-off", 0.05), 
                                         checkboxInput("ifSplineFdr","Use FDR", FALSE)),
                        textOutput("numGeneCluster"),
                        br(),
                        textInput("seleNumCluster","Number of clusters:", 5),
                        sliderInput("seleProbCut", "Cut-off for cluster memebership probability", 0, 1, 0.6, 0.05),
                        actionButton("runCluster", "Run clustering", style = "background-color:salmon"),
                        br(),
                        br(),
                        # plotSilhouette is currently disabled --------------------------------------------------------
                        # actionButton("plotSilhouette","Plot Silhouette and WSS scores", style = "background-color:salmon"),
                        # h6("(Silhouette and WSS scores may help to select optimal cluster numbers. 
                        # May take long time to calculate when number of genes is high.)")),
                      ),
                      mainPanel(tags$style(type="text/css",
                                           ".shiny-output-error { visibility: hidden; }",
                                           ".shiny-output-error:before { visibility: hidden; }"),
                                uiOutput("clusterPlotUI"),
                                uiOutput("downloadClusterTabUI"),
                                uiOutput("seleClusterBox"),
                                DT::dataTableOutput("eachClusterTab"),
                                plotOutput("clusterTimePlot",height = 350, width = 566)
                      ))),
           
           # This tab performs enrichment analysis on result from hypothesis testing
           # or time-series clustering. User can choose to perform either Pathway
           # enrichment or Phospho-signature enrichment, the latter is available
           # only for Phosphoproteome data.
           # For Pathway enrichment: the method is gene-centric, i.e., for 
           # phosphoproteome, it does not consider multiple phosphorylation sites.
           # For DE result, user can choose to perform either PAGE (Parametric
           # Analysis of Gene Set Enrichment)  or GSEA (Gene Set Enrichment Analysis)
           # using either t-statistic or logFC from the DE test. Geneset databases
           # are stored in the geneset folder (.gmt files) at the same directory
           # as the scripts of shiny app.
           # ---------------------
           # For Phospho-signature enrichment: the method is available only for 
           # phosphoproteome data and use a site-centric database in which each
           # set (signature set) has PTM sites and their direction of regulation
           # instead of just protein names. The database are from Krug et al.,2019
           # (https://doi.org/10.1074%2Fmcp.TIR118.000943) and stored locally in 
           # the ptmset folder. We use the PTM-SEA algorithm from Krug et al. to
           # analyze the result from Differential Exression analysis and Fisher's
           # exact test for the result from Time series cluster. For the latter, 
           # each signature set is split into one containing the upregulated
           # phosphosites and one containing the downregulated phosphosites.
           tabPanel("Enrichment analysis",
                    titlePanel("Enrichment analysis on result from Differential expression or Time-series clustering"),
                    sidebarLayout(
                      sidebarPanel(
                        conditionalPanel(condition = "input.assay == 'Phosphoproteome'",
                                         radioButtons("analysisMethod",
                                                      "Select analysis method",
                                                      c("Pathway enrichment", "Phospho-signature enrichment"),
                                                      inline = FALSE)),
                        conditionalPanel(condition = "input.assay == 'Proteome'",
                                         radioButtons("analysisMethod1", 
                                                      "Select analysis method",
                                                      c("Pathway enrichment"),
                                                      selected = "Pathway enrichment")),
                        radioButtons("seleSourceEnrich", "Source of gene list:",
                                     c("Differential expression", "Selected time-series cluster", "All time-series clusters")),
                        conditionalPanel(condition = "input.seleSourceEnrich == 'Differential expression'",
                                         conditionalPanel(condition = "input.analysisMethod == 'Pathway enrichment' || input.assay == 'Proteome'",
                                                          radioButtons("enrichMethod", 
                                                                       "Select enrichment method",
                                                                       c("PAGE", "GSEA"),
                                                                       inline = TRUE)),
                                         conditionalPanel(condition = "input.enrichMethod == 'GSEA' || input.analysisMethod == 'Phospho-signature enrichment'",
                                                          numericInput("permNum", 
                                                                       "Permutation number",
                                                                       value = 100,
                                                                       min = 10, max = 10000)),
                                         radioButtons("statType",
                                                      "Statistic used for ranking:",
                                                      c("stat", "log2FC"),
                                                      inline = TRUE)),
                        conditionalPanel(condition = "input.seleSourceEnrich == 'Selected time-series cluster' || input.seleSourceEnrich == 'All time-series clusters'",
                                         radioButtons("enrichMethod1",
                                                      "Select enrichment method",
                                                      c("Fisher's exact test"),
                                                      selected = "Fisher's exact test")),
                        # select geneset or PTM set
                        conditionalPanel(condition = "input.analysisMethod == 'Pathway enrichment' | input.assay == 'Proteome'",
                                         radioButtons("seleGeneSet", "Select the Gene set database",
                                                      c("select from available gene set databases", "upload a gene set database"),
                                                      selected = "select from available gene set databases"),
                                         conditionalPanel(condition = "input.seleGeneSet == 'select from available gene set databases'",
                                                          selectInput("sigSet", "Select geneset database",
                                                                      list.files(path = "geneset/", pattern = "\\.gmt$"),
                                                                      selectize = FALSE, 
                                                                      selected = "Cancer_Hallmark.gmt", size = 8)),
                                         conditionalPanel(condition = "input.seleGeneSet == 'upload a gene set database'",
                                                          fileInput("uploadGeneSet", "Upload the desired gene set file (.gmt format)"))
                        ),
                        conditionalPanel(condition = "input.analysisMethod == 'Phospho-signature enrichment' && input.assay == 'Phosphoproteome'",
                                         radioButtons("selePTMSet", "Select the PTM set database",
                                                      c("select from available PTM set databases", "upload a PTM set database"),
                                                      selected = "select from available PTM set databases"),
                                         conditionalPanel(condition = "input.selePTMSet == 'select from available PTM set databases'",
                                                          selectInput("sigSetPTM",  "Select PTM set database",
                                                                      list.files(path = "ptmset/", pattern = "\\.txt$"),
                                                                      selectize = FALSE,
                                                                      selected = "Human_PTM.txt", size = 8)),
                                         conditionalPanel(condition = "input.selePTMSet == 'upload a PTM set database'",
                                                          fileInput("uploadPTMSet", "Upload the desired PTM set file (.txt format)"))
                        ),
                        numericInput("sigLevel", label = "P value cut-off for enrichment results",
                                     min = 0, max = 1, value = 0.05),
                        checkboxInput("ifEnrichFDR", label = "use FDR", value = FALSE),
                        actionButton("RunEnrich", "Start analysis", style = "background-color:salmon"),width = 3),
                      mainPanel(tags$style(type="text/css",
                                           ".shiny-output-error { visibility: hidden; }",
                                           ".shiny-output-error:before { visibility: hidden; }"),
                                span(textOutput("errMsg"), style="color:red"),
                                conditionalPanel(condition = "input.seleSourceEnrich == 'Differential expression' || input.seleSourceEnrich == 'Selected time-series cluster'",
                                                 uiOutput("downloadUI1"),
                                                 DT::dataTableOutput("enrichTab"),
                                                 column(DT::dataTableOutput("geneTab"),
                                                        uiOutput("downloadUI2"),width = 9),
                                                 column(plotOutput("plotEnr",height = 500, width = 600), width =6)),
                                conditionalPanel(condition = " input.seleSourceEnrich == 'All time-series clusters'",
                                                 plotlyOutput("clustEnrPlot", height = "auto"),
                                                 DT::dataTableOutput("geneTabClicked"),
                                                 plotOutput("plotGeneEnr"))
                      ))),
           
           # This tab performs kinase activity inference using decoupleR. From 
           # changes in abundance of phosphorylation sites between 2 conditions,
           # it infers how active the responsible kinases are between two said 
           # conditions. A network of kinase ("source") and their phosphorylation
           # targets were constructed using databases from OmnipathR. Note that 
           # we remove interactions whose source is only ProtMapper. The network 
           # is constructed with data from either homo sapiens or mus musculus.
           # Result from hypothesis testing or time-series clustering (logFC) are
           # used to infer kinase activity with the function decoupler::run_wmean().
           # For hypothesis testing, user can choose to use either the t-statistics
           # or logFC. For time-series clustering, only logFC is used. The result
           # is displayed as barplot (hypothesis testing) or heatmap (time-series
           # clustering). Kinases whose p-values are smaller than a pre-determined
           # threshold (e.g., 0.05) are highlighted.
           tabPanel("Kinase activity inference",
                    titlePanel("Inferring kinase activity from hypothesis testing or time-series clustering result"),
                    sidebarLayout(
                      sidebarPanel(
                        radioButtons("seleSourceKinase", "Select result to perform the analysis on:", c("Differential expression", "Time-series cluster")),
                        conditionalPanel(condition = "input.seleSourceKinase == 'Differential expression'",
                                         radioButtons("statTypeKinase", "Statistic used for computing kinase score:", c("stat", "log2FC"), inline = TRUE),
                                         numericInput("nTopKinase", "Number of top kinases in plot", min = 1, max = 100, value = 10)),
                        conditionalPanel(condition = "input.seleSourceKinase == 'Time-series cluster'",
                                         radioButtons("seleKinaseTimeMethod", "Infer whether how active the kinases are or how likely the kinases are associated to the cluster",
                                                      c("activity", "association"), inline = TRUE),
                                         conditionalPanel(condition = "input.seleKinaseTimeMethod == 'association'",
                                                          radioButtons("seleAssoMethod", "Select method",
                                                                       c("Fisher's exact test", "FGSEA")))),
                        conditionalPanel(condition = "input.seleAssoMethod == 'FGSEA' || input.seleKinaseTimeMethod == 'activity' || input.seleSourceKinase == 'Differential expression'",
                                         numericInput("nPermKinase", "Number of permutations", min = 10, max = 10000, value = 100)),
                        selectInput("speciesRef", "Select reference species", c("Homo sapiens", "Mus musculus")),
                        numericInput("pKinase", label = "Highlight kinases with p-values under:", min = 0, max = 1, value = 0.05),
                        checkboxInput("ifKinaseFDR", label = "use FDR", value = FALSE),
                        actionButton("runKinase", "Start analysis", style = "background-color:salmon"), width = 3),
                      mainPanel(tags$style(type="text/css",
                                           ".shiny-output-error { visibility: hidden; }",
                                           ".shiny-output-error:before { visibility: hidden; }"),
                                span(textOutput("errMsgKinase"), style="color:red"),
                                column(DT::dataTableOutput("kinaseTab"),
                                       uiOutput("downloadUI3"),width = 6),
                                column(plotOutput("plotKinase",height = 800, width = 600),width =6),
                                DT::dataTableOutput("pSiteTab"),
                                uiOutput("downloadUI4")
                                
                      ))),
           
           # This tab displays the information about all the selections and inputs 
           # by the user. It also allows to download all the log information in a 
           # text/tsv file.
           tabPanel("Log Info",
                    titlePanel("Log information of the inputs and selections"),
                    mainPanel(downloadButton('downloadLogValues', 'Download log values'),
                              hr(),
                              DT::dataTableOutput("show_inputs"))
           ),

           # This tab displays the update notes of the SmartPhos explorer app.
           tabPanel("Update Notes",
                    includeMarkdown("./NEWS.md"))
)
