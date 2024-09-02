#' @name generateInputTable
#'
#' @title Generate Input Table for Proteomic and Phosphoproteomic Analysis
#'
#' @description
#' \code{generateInputTable} generates an input table for proteomic and phosphoproteomic analysis by reading files from a specified folder.
#'
#' @param rawFolder A \code{character} string specifying the path to the folder containing the raw files.
#' @param batchAsFolder A \code{logical} value indicating whether to treat subdirectories as separate batches. Default is \code{FALSE}.
#' @return A \code{data.frame} with columns fileName, sample, type, batch, and id that can be used as input for further analysis.
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Optionally treats subdirectories as separate batches.
#'   \item Reads the summary file containing experimental information.
#'   \item Generates unique experimental IDs based on batch and sample names.
#'   \item Processes file paths for full proteome and phosphoproteome data.
#'   \item Creates a combined input table with file names, sample names, types, batches, and IDs.
#' }
#'
#' @examples
#' # inputTable <- generateInputTable(folder, batchAsFolder = FALSE)
#'
#'
#' @importFrom utils read.delim
#' @export
generateInputTable <- function(rawFolder, batchAsFolder = FALSE) {

    if (batchAsFolder) {
        # Get all subdirectories in the rawFolder
        allFolder <- list.dirs(rawFolder, full.names = FALSE)
        # Remove any empty or NA values from the list of folders
        allFolder <- allFolder[!allFolder %in% c(NA,"")]
    } else {
        allFolder <- ""
    }

    # Process each folder
    inputTab <- lapply(allFolder, function(eachFolder) {
        folderPath <- file.path(rawFolder, eachFolder)
        # Read summary file
        sumFile <- read.delim(file.path(folderPath,"summary.txt"))
        # Get unique sample names from the summary file
        sampleName <- unique(sumFile$Experiment)
        # Replace special characters in sample names
        sampleName <- gsub("[+]", ".", sampleName)
        # Remove any empty or NA sample names
        sampleName <- sampleName[!sampleName %in% c(NA,"")]

        # Create a data frame for full proteome files
        fullTab <- data.frame(fileName = rep(file.path(folderPath,"proteinGroups.txt"),length(sampleName)))
        fullTab$sample <- sampleName
        fullTab$type <- "proteome"

        # Create a data frame for phosphoproteome files
        phosphoTab <- data.frame(fileName = rep(file.path(folderPath,"Phospho(STY)Sites.txt"),length(sampleName)))
        phosphoTab$sample <- sampleName
        phosphoTab$type <- "phosphoproteome"

        # Combine the full proteome and phosphoproteome tables
        inputTab <- rbind(fullTab, phosphoTab)
        inputTab$batch <- eachFolder
        # Create a unique experimental ID
        inputTab$id <- paste0(inputTab$batch,"_",inputTab$sample)
        inputTab
    })
    # Combine all the input tables from each folder into one
    inputTab <- do.call(rbind, inputTab)

    return(inputTab)
}


#' @name generateInputTable_DIA
#'
#' @title Generate Input Table for DIA Analysis
#'
#' @description
#' \code{generateInputTable_DIA} generates an input table for DIA analysis by reading files from a specified folder.
#'
#' @param rawFolder A \code{character} string specifying the path to the folder containing the raw files.
#' @return A \code{data.frame} with columns fileName, type, and id that can be used as input for further analysis.
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Reads the summary file containing experimental information.
#'   \item Generates unique experimental IDs based on sample type, treatment, timepoint, and replicate.
#'   \item Processes file paths for full proteome and phosphoproteome data.
#'   \item Creates a combined input table with file names, types, and IDs.
#' }
#'
#' @examples
#' # inputTable <- generateInputTable_DIA(folder)
#' @importFrom utils read.delim
#' @export
generateInputTable_DIA <- function(rawFolder) {
    # List all files in the rawFolder with pattern "*Labels.txt"
    getFile <- list.files(rawFolder, pattern = "*Labels.txt")

    # Read the summary file containing information about the experiments
    sumFile <- read.delim(file.path(rawFolder, getFile))

    # Create a unique experimental ID
    id <- c()
    for (i in 1: length(sumFile$Sample.Type)) {
        id <- c(id, paste(sumFile$Sample.Type[i], sumFile$treatment[i],
                          sumFile$timepoint[i], sumFile$replicate[i], sep = "_"))
    }

    # Replace special characters in IDs
    id <- gsub("[+]", ".", id)

    # List all files in the rawFolder with pattern "*Protein Report.xls"
    getFileFP <- list.files(rawFolder, pattern = "*Protein Report.xls")
    # Create a data frame for full proteome files
    fullTab <- data.frame(fileName = rep(file.path(rawFolder, getFileFP), length(id)))
    fullTab$type <- "proteome"
    fullTab$id <- id

    # List all files in the rawFolder with pattern "*PTM Report2.xls"
    getFilePhospho <- list.files(rawFolder, pattern = "*PTM Report2.xls")
    # Create a data frame for phosphoproteome files
    phosphoTab <- data.frame(fileName = rep(file.path(rawFolder, getFilePhospho), length(id)))
    phosphoTab$type <- "phosphoproteome"
    phosphoTab$id <- id

    # Combine the full proteome and phosphoproteome tables into a single input table
    inputTab <- rbind(fullTab, phosphoTab)

    return(inputTab)
}


#' @name readExperiment
#'
#' @title Read and Process the DDA experiment.
#'
#' @description
#' \code{readExperiment} reads and processes DDA (Data-Dependent Acquisition) phosphoproteomic and proteomic data from a given file table, and returns a \code{MultiAssayExperiment} object.
#'
#' @param fileTable A \code{data.frame} containing information about the input files, including type, id, sample, and other annotations.
#' @param localProbCut \code{Numeric}, local probability cutoff for filtering phosphoproteomic data. Default is 0.75.
#' @param scoreDiffCut \code{Numeric}, score difference cutoff for filtering phosphoproteomic data. Default is 5.
#' @param fdrCut \code{Numeric}, false discovery rate cutoff for filtering proteomic data. Default is 0.1.
#' @param scoreCut \code{Numeric}, score cutoff for filtering proteomic data. Default is 10.
#' @param pepNumCut \code{Numeric}, peptide number cutoff for filtering proteomic data. Default is 1.
#' @param ifLFQ \code{Logical}, whether to use LFQ quantification for proteomic data. Default is \code{TRUE}.
#' @param annotation_col A \code{character} vector specifying additional columns to be included in the sample annotation. Default is an empty vector.
#'
#' @return A \code{MultiAssayExperiment} object containing the processed phosphoproteomic and proteomic data from a DDA experiment.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Reads and processes the phosphoproteomic data using the \code{readPhosphoExperiment} function.
#'   \item Reads and processes the proteomic data using the \code{readProteomeExperiment} function.
#'   \item Prepares the sample annotation table.
#'   \item Constructs and returns a \code{MultiAssayExperiment} object containing the processed data.
#' }
#'
#' @examples
#' # Example usage:
#' fileTable <- data.frame(id = c("sample1", "sample2"), sample = c("A", "B"), type = c("proteome", "proteome"))
#' mae <- readExperiment(fileTable)
#'
#' @import MultiAssayExperiment
#' @export
readExperiment <- function(fileTable, localProbCut = 0.75, scoreDiffCut = 5, fdrCut =0.1, scoreCut = 10, pepNumCut = 1, ifLFQ = TRUE, annotation_col = c()) {
    # Read phosphoproteomic data
    print("Processing phosphoproteomic data")
    ppe <- readPhosphoExperiment(fileTable, localProbCut, scoreDiffCut)

    # Read full proteome data
    print("Processing proteomic data")
    fpe <- readProteomeExperiment(fileTable, fdrCut, scoreCut, pepNumCut, ifLFQ)

    # Prepare sample annotation
    if ("batch" %in% colnames(fileTable)) annotation_col <- c(annotation_col,"batch")
    sampleTab <- fileTable[,c("id","sample", annotation_col)]
    sampleTab <- sampleTab[!duplicated(sampleTab$id),]
    rownames(sampleTab) <- sampleTab$id
    sampleTab$id <- NULL

    # Construct MultiAssayExperiment object
    if (!is.null(ppe) & !is.null(fpe)) {
      mae <- MultiAssayExperiment(list(Phosphoproteome = ppe, Proteome = fpe),
                                  sampleTab)
    }
    else if (is.null(ppe)) {
      mae <- MultiAssayExperiment(list(Proteome = fpe), sampleTab)
    }
    else {
      mae <- MultiAssayExperiment(list(Phosphoproteome = ppe), sampleTab)
    }
    return(mae)
}


#' @name readExperimentDIA
#'
#' @title Read and Process a DIA Experiment
#'
#' @description
#' \code{readExperimentDIA} reads and processes DIA (Data-Independent Acquisition) data for both phosphoproteome and proteome experiments, and constructs a \code{MultiAssayExperiment} object.
#'
#' @param fileTable A \code{data frame} containing metadata about the files to be read. Must contain columns type, fileName, id, and optionally outputID.
#' @param localProbCut \code{Numeric}, the local probability cutoff for phosphoproteomic data. Default is 0.75.
#' @param annotation_col A \code{character} vector specifying the columns in fileTable to be used for sample annotation.
#' @param onlyReviewed A \code{logical} value indicating whether to include only reviewed proteins. Default is \code{TRUE}.
#' @param normalizeByProtein \code{Logical}, whether to normalize the data by protein. Default is \code{FALSE}.
#' @return A \code{MultiAssayExperiment} object containing the processed phosphoproteome and proteome data.
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Reads and processes phosphoproteomic data using \code{readPhosphoExperimentDIA}.
#'   \item Reads and processes proteomic data using \code{readProteomeExperimentDIA}.
#'   \item Prepares sample annotations based on the provided fileTable and annotation_col.
#'   \item Constructs a \code{MultiAssayExperiment} object with the processed data and sample annotations.
#' }
#' The \code{readPhosphoExperimentDIA} and \code{readProteomeExperimentDIA} functions are used to read and filter the data for phosphoproteome and proteome experiments, respectively, and they must be available in the environment.
#'
#' @examples
#' # Example usage:
#' fileTable <- data.frame(type = c("phosphoproteome", "proteome"),
#'                         fileName = c("phos_file.txt", "prot_file.txt"),
#'                         id = c("sample1", "sample2"),
#'                         outputID = c("s1", "s2"))
#' result <- readExperimentDIA(fileTable, localProbCut = 0.75, annotation_col = c("id"), normalizeByProtein = FALSE)
#'
#' @import MultiAssayExperiment
#' @export
readExperimentDIA <- function(fileTable, localProbCut = 0.75, annotation_col = c(), onlyReviewed = TRUE, normalizeByProtein=FALSE) {

    # Read phospho data
    print("Processing phosphoproteomic data")
    ppe <- readPhosphoExperimentDIA(fileTable, localProbCut, onlyReviewed)
    print("Successful!!")

    # Read full proteome data
    print("Processing proteomic data")
    fpe <- readProteomeExperimentDIA(fileTable)
    print("Successful!!")

    # Use user-specified output sample ID if available
    if("outputID" %in% colnames(fileTable)) {
       fileTable$id <- fileTable$outputID
    }

    # Prepare sample annotation
    sampleTab <- fileTable[,c("id", annotation_col),drop=FALSE]
    sampleTab <- sampleTab[!duplicated(sampleTab$id),,drop=FALSE]

    # Convert sample table to data frame
    sampleTab <- as.data.frame(sampleTab)
    rownames(sampleTab) <- sampleTab$id
    sampleTab$sample <- sampleTab$id
    sampleTab$id <- NULL

    # Get sample name without prefix or suffix
    sampleTab$sampleName <- removePreSuffix(sampleTab$sample)

    # Construct MultiAssayExperiment object
    if (!is.null(ppe) & !is.null(fpe)) {
      mae <- MultiAssayExperiment(list(Phosphoproteome = ppe, Proteome = fpe), sampleTab)
    } else if (is.null(ppe)) {
      mae <- MultiAssayExperiment(list(Proteome = fpe), sampleTab)
    } else {
      mae <- MultiAssayExperiment(list(Phosphoproteome = ppe), sampleTab)
    }

    return(mae)
}


#' @name normByFullProteome
#'
#' @title Normalize Phosphoproteome by Full Proteome
#'
#' @description
#' \code{normByFullProteome} normalizes the phosphoproteome data by the corresponding full proteome data in a \code{MultiAssayExperiment} object. The "Phosphoproteome" assay
#' in the MultiAssayExperiment will be replaced by the ratio.
#'
#' @param mae A \code{MultiAssayExperiment} object containing both phosphoproteome and proteome assays.
#' @param replace \code{Logical}, whether to replace the existing phosphoproteome assay in the \code{MultiAssayExperiment} object. If replace = \code{FALSE}, a new assay, phosphoRatio, will be created. Default is \code{TRUE}.
#' @return A \code{MultiAssayExperiment} object with the normalized phosphoproteome data.
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Checks if both phosphoproteome and proteome assays are present in the \code{MultiAssayExperiment} object.
#'   \item Extracts the phosphoproteome and proteome assays along with the sample annotations.
#'   \item Matches the samples between the phosphoproteome and proteome assays.
#'   \item Normalizes the phosphoproteome data by dividing it by the corresponding proteome data.
#'   \item Replaces the phosphoproteome assay in the \code{MultiAssayExperiment} object or adds the normalized data as a new assay, depending on the \code{replace} parameter.
#' }
#'
#'
#' @import MultiAssayExperiment
#' @export
normByFullProteome <- function(mae, replace = TRUE) {

    # Check if both Phosphoproteome and Proteome assays are present
    if (!all(c("Phosphoproteome","Proteome") %in% names(assays(mae)))) {
        stop("Both Phosphoproteome and Proteome assays should be present in the MultiAssayExperiment object")
    }

    # Extract phosphoproteome and proteome assays
    ppe <- mae[["Phosphoproteome"]]
    fpe <- mae[["Proteome"]]
    sampleTab <- colData(mae)
    sampleTab.pp <- sampleTab[sampleTab$sampleType %in% c("Phospho","PP"),]
    sampleTab.fp <- sampleTab[sampleTab$sampleType %in% c("FullProteome","FP"),]

    # Check if proteome assay for the unenriched samples is present
    if (nrow(sampleTab.fp) ==0 ) {
        stop("Proteome assay for the unenriched samples i.e., sampleType with FullProteome should be present")
    }

    # Extract assay data for phosphoproteome and proteome
    ppMat <- assay(ppe[,colnames(ppe) %in% rownames(sampleTab.pp)])
    fpMat <- assay(fpe[,colnames(fpe) %in% rownames(sampleTab.fp)])

    # Find the full proteome (non-enriched) counterparts of the phosphoproteome (enriched) samples
    ppSampleName <- sampleTab.pp[colnames(ppMat),]$sampleName
    fpSampleID <- rownames(sampleTab.fp)[match(ppSampleName, sampleTab.fp$sampleName)]

    # Handle the situation where not every phosphoproteome sample has a full proteome counterpart
    ppMat <- ppMat[,!is.na(fpSampleID)]
    fpSampleID <- fpSampleID[!is.na(fpSampleID)]

    # Get the full proteomic measurement of the corresponding proteins and samples
    fpMat <- fpMat[match(rowData(ppe)$UniprotID, rowData(fpe)$UniprotID),fpSampleID]

    # Ensure the samples match
    stopifnot(all(sampleTab[(colnames(fpMat)),]$sampleName == sampleTab[(colnames(ppMat)),]$sampleName))

    # Normalize phosphoproteome data by the full proteome data
    ppMat.norm <- ppMat/fpMat
    ppMat.norm <- ppMat.norm[rowSums(!is.na(ppMat.norm))>0,]
    ppe.norm <- ppe[rownames(ppMat.norm),colnames(ppMat.norm)]
    assay(ppe.norm) <- ppMat.norm

    # Replace or add the normalized phosphoproteome data to the MultiAssayExperiment object
    if (replace) {
        mae[["Phosphoproteome"]] <- ppe.norm
    } else {
        mae <- MultiAssayExperiment(list(Phosphoproteome = mae[["Phosphoproteome"]],
                                             Proteome = mae[["Proteome"]],
                                             PhosphoRatio = ppe.norm),
                                             colData = colData(mae))
    }

    return(mae)
}

