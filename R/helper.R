#' @name readOnePhos
#'
#' @title Read and Filter Phosphorylation Data for a Specific Sample
#'
#' @description
#' `readOnePhos` reads phosphorylation data from an input table, filters it based on localization probability, score difference, and intensity, and returns the filtered data for a specific sample.
#'
#' @param inputTab A data.table or data.frame containing phosphorylation data with columns for localization probability, score difference, and intensity for various samples.
#' @param sampleName A character string specifying the sample name to filter data for.
#' @param localProbCut A numeric value specifying the cutoff for localization probability. Default is `0.75`.
#' @param scoreDiffCut A numeric value specifying the cutoff for score difference. Default is `5`.
#' @param multiMap A logical value indicating whether to allow multiple mapping (not used in this function but could be relevant for further extensions).
#'
#' @return A data.frame containing the filtered phosphorylation data for the specified sample, with columns for intensity, Uniprot ID, gene name, position within proteins, amino acid residue, and sequence window.
#'
#' @details
#' The function filters the input phosphorylation data based on three criteria: localization probability, score difference, and intensity. Only rows that meet or exceed the specified cutoffs for these criteria and have non-zero intensity are retained. The filtered data is then returned with a unique identifier for each row.
#'
#' @import data.table
#' @examples
#' # inputTab <- data.table::fread("phosphorylation_data.csv")
#' # filteredData <- readOnePhos(inputTab, sampleName = "Sample1", localProbCut = 0.75, scoreDiffCut = 5, multiMap = FALSE)
#'
#' @export
readOnePhos <- function(inputTab, sampleName, localProbCut = 0.75, scoreDiffCut = 5, multiMap) {

  # Define sample specific column names
  colSele <- paste0(c("Localization.prob.","Score.diff.","Intensity."), sampleName)
  # Check if all required columns are present in the input table
  if (!all(colSele %in% colnames(inputTab))) stop("Sample not found in quantification file")

  # Filter rows based on localization probability, score difference, and non-zero intensity
  keepRow <- (!is.na(inputTab[[colSele[1]]]) & inputTab[[colSele[1]]] >= localProbCut) &
    (!is.na(inputTab[[colSele[2]]]) & inputTab[[colSele[2]]] >= scoreDiffCut) &
    (!is.na(inputTab[[colSele[3]]]) & inputTab[[colSele[3]]]>0)

  # If no rows pass the filter, return a warning and NULL
  if (all(!keepRow)) {
    warning(sprintf("sample %s does not contain any records after filtering",sampleName))
    return(NULL)
  }

  # Subset the input table based on the filtered rows and select relevant columns
  outputTab <- inputTab[keepRow,
                        c(colSele[3],"Proteins","Gene.names",
                          "Positions.within.proteins","Amino.acid","Sequence.window"),
                        with=FALSE]

  # Create a unique identifier for each row by concatenating protein ID and position
  outputTab$rowName <- paste0(outputTab$Proteins, "_",
                              outputTab$Positions.within.proteins)

  # Ensure that all identifiers are unique
  stopifnot(all(!duplicated(outputTab$rowName)))

  # Set the row names of the output table and remove the rowName column
  rownames(outputTab) <- outputTab$rowName
  outputTab$rowName <- NULL
  # Rename columns to more meaningful names
  colnames(outputTab) <- c("Intensity","UniprotID","Gene","Position","Residue","Sequence")

  return(outputTab)
}


#' @name readPhosphoExperiment
#'
#' @title Read Phosphorylation Experiment Data
#'
#' @description
#' `readPhosphoExperiment` reads and processes phosphorylation experiment data from multiple files, filtering based on localization probability and score difference, and constructs a `SummarizedExperiment` object.
#'
#' @param fileTable A data.table or data.frame containing information about the files, including columns for file names, sample names, and other relevant metadata. It must include a column named `type` with value "phosphoproteome" for relevant entries.
#' @param localProbCut A numeric value specifying the cutoff for localization probability. Default is `0.75`.
#' @param scoreDiffCut A numeric value specifying the cutoff for score difference. Default is `5`.
#'
#' @return A `SummarizedExperiment` object containing the processed phosphorylation data.
#'
#' @details
#' This function reads phosphorylation data from multiple files as specified in `fileTable`, filters the data based on localization probability and score difference, and removes reverse and potential contaminant entries. It constructs an intensity matrix and annotation data, which are then used to create a `SummarizedExperiment` object.
#'
#' @import data.table SummarizedExperiment
#' @examples
#' # fileTable <- data.table::fread("file_information.csv")
#' # ppe <- readPhosphoExperiment(fileTable, localProbCut = 0.75, scoreDiffCut = 5)
#'
#' @export
readPhosphoExperiment <- function(fileTable, localProbCut = 0.75, scoreDiffCut = 5) {


  # Select only phosphoproteomic entries
  fileTable <- fileTable[fileTable$type == "phosphoproteome",]
  if (nrow(fileTable) == 0) {
    return(NULL)
  }

  # Read in all batches and store them in a list
  expAll <- lapply(unique(fileTable$fileName), function(eachFileName) {
    fileTableSub <- fileTable[fileTable$fileName == eachFileName,]

    # Read input table with tab as delimiter
    inputTab <- data.table::fread(eachFileName, check.names = TRUE)

    # Remove empty features
    inputTab <- inputTab[!inputTab$Proteins %in% c(NA,"") &
                           #!inputTab$Gene.names %in% c("",NA) &
                           !inputTab$Positions.within.proteins %in% c(NA,""),]

    # Remove reverse and potential contaminants from the whole table
    inputTab <- inputTab[!inputTab$Potential.contaminant %in% "+" &
                           !inputTab$Reverse %in% "+",]

    # Stop if none of the records pass the chosen threshold
    if (nrow(inputTab) == 0) stop("No phosphorylation site could pass the specified threshold in any sample!")

    # Process data for each sample
    expSub <- lapply(seq(nrow(fileTableSub)), function(i) {
      eachTab <- readOnePhos(inputTab,
                             fileTableSub[i,]$sample,
                             localProbCut, scoreDiffCut)
      if (!is.null(eachTab)) {
        eachTab$id <- fileTableSub[i,]$id
        eachTab$rowName <- rownames(eachTab)
      }
      eachTab
    })
    expSub <- data.table::rbindlist(expSub)
    expSub
  })

  expAll <- data.table::rbindlist(expAll) # rbindlist is faster than do.call(rbind)

  # Stop if none of the record passed chosen threshold
  if (nrow(expAll) == 0) stop("No phosphorylation site could pass the specified threshold in any sample!")

  # Prepare annotations
  annoTab <- expAll[!duplicated(expAll$rowName),c("rowName","UniprotID",
                                                  "Gene","Position","Residue","Sequence")]
  annoTab$site <- annoTab$rowName
  rownames(annoTab) <- annoTab$rowName
  annoTab$rowName <- NULL

  # Prepare intensity matrix
  phosMat <- matrix(data = rep(NA, nrow(annoTab)*nrow(fileTable)),
                    nrow(annoTab), nrow(fileTable),
                    dimnames = list(rownames(annoTab),fileTable$id))
  for (each_id in fileTable$id) {
    eachTab <- expAll[(expAll$id %in% each_id),]
    phosMat[,each_id] <- as.numeric(eachTab[match(rownames(phosMat),eachTab$rowName),][["Intensity"]])
  }

  # Rename rownames
  rownames(phosMat) <- rownames(annoTab) <- paste0("s",seq(nrow(annoTab)))

  # Construct SummarizedExperiment object
  ppe <- SummarizedExperiment(assays = list(Intensity = phosMat),
                              rowData = annoTab)
  return(ppe)
}


#' @name readOnePhosDIA
#'
#' @title Read Phosphorylation Data for One Sample from DIA
#'
#' @description
#' `readOnePhosDIA` reads and processes phosphorylation data for a single sample from a DIA experiment, applying filters for localization probability and removing duplicates if specified.
#'
#' @param inputTab A data.table or data.frame containing phosphorylation data.
#' @param sampleName A character string specifying the sample name.
#' @param localProbCut A numeric value specifying the cutoff for localization probability. Default is `0.75`.
#' @param removeDup A logical value indicating whether to remove duplicate entries based on `UniprotID` and intensity. Default is FALSE.
#'
#' @return A data.table containing the processed phosphorylation data for the specified sample.
#'
#' @details
#' This function processes phosphorylation data for a single sample by filtering based on localization probability and non-zero intensity. It handles multiplicity by summarizing intensities and optionally removes duplicates. The resulting data is returned as a data.table with unique identifiers.
#'
#' @import data.table
#' @importFrom stats aggregate
#' @examples
#' # inputTab <- data.table::fread("phosphorylation_data.tsv")
#' # result <- readOnePhosDIA(inputTab, "Sample_1", localProbCut = 0.75, removeDup = TRUE)
#'
#' @export
readOnePhosDIA <- function(inputTab, sampleName, localProbCut = 0.75, removeDup = FALSE) {

  sampleName <- make.names(sampleName) # Ensure sample name is syntactically valid

  # Define sample-specific column names
  colSele <- c(NA,NA) # Placeholder for localization probability and quantity columns
  if (length(grep(pattern = paste0("*", sampleName, ".*PTM.SiteProbability"), colnames(inputTab))) > 0) {
    colSele[1] <- colnames(inputTab)[grep(pattern = paste0("*", sampleName, ".*PTM.SiteProbability"), colnames(inputTab))]
  } else {
    stop("Sample not found in quantification file")
  }
  if (length(grep(pattern = paste0("*", sampleName, ".*PTM.Quantity"), colnames(inputTab))) > 0) {
    colSele[2] <- colnames(inputTab)[grep(pattern = paste0("*", sampleName, ".*PTM.Quantity"), colnames(inputTab))]
  } else {
    stop("Sample not found in quantification file")
  }

  # Replace "Filtered" values with NA
  inputTab[[colSele[1]]][inputTab[[colSele[1]]] == "Filtered"] <- NA
  inputTab[[colSele[2]]][inputTab[[colSele[2]]] == "Filtered"] <- NA

  # Convert to numeric values, handling commas as decimal points
  inputTab[[colSele[1]]] <- as.numeric(gsub(",", ".", inputTab[[colSele[1]]]))
  inputTab[[colSele[2]]] <- as.numeric(gsub(",", ".", inputTab[[colSele[2]]])) # Change , to . if any. Sometimes the "," is used as decimal.

  # Get features passing quality filters and non-zero intensity
  keepRow <- (!is.na(inputTab[[colSele[1]]]) & inputTab[[colSele[1]]] >= localProbCut) &
    (!is.na(inputTab[[colSele[2]]]) & inputTab[[colSele[2]]]>0)

  # Check if any rows remain after filtering
  if (all(!keepRow)) {
    warning(sprintf("sample %s does not contain any records after filtering", sampleName))
    return(NULL)
  }

  # Subset the table based on the filters
  if ("PTM.Multiplicity" %in% colnames(inputTab)) { #Multiplicity is reported by Spectronaut
    outputTab <- inputTab[keepRow,
                          c(colSele[2],"PTM.CollapseKey","PG.UniProtIds","PG.Genes","PTM.Multiplicity",
                            "PTM.SiteLocation","PTM.SiteAA","PTM.FlankingRegion"), with=FALSE]
    # Rename columns
    colnames(outputTab) <- c("Intensity","CollapseKey","UniprotID","Gene","Multiplicity","Position","Residue","Sequence")

  } else { #Multiplicity not reported
    outputTab <- inputTab[keepRow,
                          c(colSele[2],"PTM.CollapseKey","PG.UniProtIds","PG.Genes",
                            "PTM.SiteLocation","PTM.SiteAA","PTM.FlankingRegion"), with=FALSE]
    # Rename columns
    colnames(outputTab) <- c("Intensity","CollapseKey","UniprotID","Gene","Position","Residue","Sequence")
  }



  if (length(grep("PTM.Multiplicity", colnames(inputTab))) != 0) {
    # Handle multiplicity
    outputTab$CollapseKey <-  gsub("_M\\d+","",outputTab$CollapseKey) #it's safer to use regular expression to remove the suffix _M together with the numbers.

    # Summarize multiplicity
    outputTab <- outputTab[order(outputTab$Multiplicity, decreasing = TRUE),]
    outputTabNew <- outputTab[!duplicated(outputTab$CollapseKey),]
    intensityTab <- aggregate(Intensity ~ CollapseKey, outputTab, sum)
    outputTabNew$Intensity <- intensityTab[match(outputTabNew$CollapseKey, intensityTab$CollapseKey),]$Intensity
    outputTab <- outputTabNew
  }

  if (removeDup) {
    # Remove duplicates
    outputTab.rev <- outputTab[rev(seq(nrow(outputTab))),]
    outputTab.rev <- outputTab.rev[!duplicated(paste0(outputTab.rev$UniprotID,"_",
                                                      outputTab.rev[[colSele[[2]]]])),]
    outputTab <- outputTab.rev[rev(seq(nrow(outputTab.rev))),]

    # Create a unique identifier for row names
    outputTab$rowName <- paste0(outputTab$UniprotID, "_",
                                outputTab$Position)
  } else {
    outputTab$rowName <- outputTab$CollapseKey
  }

  # Ensure identifiers are unique
  stopifnot(all(!duplicated(outputTab$rowName)))

  # Set row names and remove temporary identifier column
  rownames(outputTab) <- outputTab$rowName
  outputTab$rowName <- NULL

  return(outputTab)
}


#' @name readPhosphoExperimentDIA
#'
#' @title Read Phosphorylation Experiment Data from DIA
#'
#' @description
#' `readPhosphoExperimentDIA` reads and processes phosphorylation data from DIA experiments, applying filters for localization probability, and optionally including only reviewed proteins. It constructs a `SummarizedExperiment` object.
#'
#' @param fileTable A data.table or data.frame containing metadata about the files to read. Must include columns `fileName`, `type`, and optionally `outputID`.
#' @param localProbCut A numeric value specifying the cutoff for localization probability. Default is `0.75`.
#' @param onlyReviewed A logical value indicating whether to include only reviewed proteins. Default is `TRUE`.
#' @param showProgressBar A logical value indicating whether to show a progress bar. Default is `FALSE`.
#'
#' @return A `SummarizedExperiment` object containing the processed phosphorylation data.
#'
#' @details
#' This function processes phosphorylation data from DIA experiments by filtering based on localization probability and non-zero intensity, handling multiplicity, and optionally including only reviewed proteins. The resulting data is returned as a SummarizedExperiment object with annotations and an intensity matrix.
#'
#' @import data.table
#' @import BiocParallel
#' @import SummarizedExperiment
#' @importFrom utils data
#' @examples
#' # fileTable <- data.table::data.table(fileName = c("file1.tsv", "file2.tsv"), type = "phosphoproteome")
#' # result <- readPhosphoExperimentDIA(fileTable, localProbCut = 0.75, onlyReviewed = TRUE, showProgressBar = TRUE)
#'
#' @export
readPhosphoExperimentDIA <- function(fileTable, localProbCut = 0.75, onlyReviewed = TRUE,
                                     showProgressBar = FALSE) {
  # Select phosphoproteomic entries
  fileTable <- fileTable[fileTable$type == "phosphoproteome",]
  if (nrow(fileTable) == 0) {
    return(NULL)
  }

  # Read in all batches and store them in a list
  expAll <- lapply(unique(fileTable$fileName), function(eachFileName) {
    fileTableSub <- fileTable[fileTable$fileName == eachFileName,]

    # Read input table, "\t" as delimiter, fread is faster than read.delim
    inputTab <- data.table::fread(eachFileName, check.names = TRUE)
    # Keep only Phospho (STY) modifications
    inputTab <- inputTab[inputTab$PTM.ModificationTitle == "Phospho (STY)",]

    # Decide whether "PG.ProteinGroups" or "PG.UniProtIds" should be used as protein IDs
    if ("PG.ProteinGroups" %in% colnames(inputTab) & !"PG.UniProtIds" %in% colnames(inputTab)) {
      inputTab$PG.UniProtIds <- inputTab$PG.ProteinGroups
    } else if (!"PG.ProteinGroups" %in% colnames(inputTab) & !"PG.UniProtIds" %in% colnames(inputTab)) {
      stop("Either PG.ProteinGroups or PG.UniProtIds should be in the quantification table")
    }

    # Handle missing PTM.CollapseKey column
    if (!"PTM.CollapseKey" %in% colnames(inputTab)) {
      inputTab$PTM.CollapseKey <- paste0(inputTab$PG.UniProtIds, "_", inputTab$PTM.SiteAA, inputTab$PTM.SiteLocation)
    }

    # Handle missing PG.Genes column
    if (!"PG.Genes" %in% colnames(inputTab)) {
      inputTab$PG.Genes <- NA
    }

    # Remove empty features
    inputTab <- inputTab[!inputTab$PG.UniProtIds %in% c(NA,"") &
                           #!inputTab$Gene.names %in% c("",NA) &
                           !inputTab$PTM.SiteLocation %in% c(NA,""),]

    # Stop if no records passed the chosen threshold
    if (nrow(inputTab) == 0) stop("No phosphorylation site could pass the specified threshold in any sample!")

    # Include only reviewed proteins if specified
    if (onlyReviewed) {
      data("swissProt")
      inputTab <- inputTab[inputTab$PTM.ProteinId %in% swissProt$Entry,]
    }

    # Process each sample
    expSub <- BiocParallel::bplapply(seq(nrow(fileTableSub)), function(i) {
      eachTab <- readOnePhosDIA(inputTab = inputTab,
                                sampleName = fileTableSub[i,]$id,
                                localProbCut = localProbCut)

      if (!is.null(eachTab)) {
        # Use user-specified output sample IDs if available
        if ("outputID" %in% colnames(fileTableSub)) {
          eachTab$id <- fileTableSub[i,]$outputID
        } else {
          eachTab$id <- fileTableSub[i,]$id
        }
        eachTab$rowName <- rownames(eachTab)
      }
      eachTab
    }, BPPARAM = BiocParallel::MulticoreParam(progressbar = showProgressBar))
    expSub <- data.table::rbindlist(expSub) #faster than do.call
    expSub
  })
  expAll <- data.table::rbindlist(expAll)

  # Stop if no records passed chosen threshold
  if (nrow(expAll) == 0) stop("No phosphorylation site could pass the specified threshold in any sample!")

  expAll <- expAll[order(expAll$rowName),]

  # Prepare annotations
  if ("Multiplicity" %in% colnames(expAll)) {
    annoTab <- expAll[!duplicated(expAll$rowName),c("rowName","UniprotID",
                                                    "Gene","Multiplicity","Position","Residue","Sequence")]
  }
  else {
    annoTab <- expAll[!duplicated(expAll$rowName),c("rowName","UniprotID",
                                                    "Gene","Position","Residue","Sequence")]
  }
  annoTab$site <- annoTab$rowName
  rownames(annoTab) <- annoTab$rowName
  annoTab$rowName <- NULL

  # Prepare intensity matrix
  if ("outputID" %in% colnames(fileTable)) { #use user specific sample ID
    sampleID <- fileTable$outputID
  } else {
    sampleID <- fileTable$id
  }
  phosMat <- matrix(data = rep(NA, nrow(annoTab)*nrow(fileTable)),
                    nrow(annoTab), nrow(fileTable),
                    dimnames = list(rownames(annoTab),sampleID))
  for (each_id in sampleID) {
    eachTab <- expAll[(expAll$id %in% each_id),]
    phosMat[,each_id] <- as.numeric(eachTab[match(rownames(phosMat),eachTab$rowName),][["Intensity"]])
  }

  # Rename rownames
  rownames(phosMat) <- rownames(annoTab) <- paste0("s",seq(nrow(annoTab)))

  # Construct SummarizedExperiment object
  ppe <- SummarizedExperiment(assays = list(Intensity = phosMat),
                              rowData = annoTab)
  return(ppe)
}



#' @name readOneProteom
#'
#' @title Read and Process One Proteomics Sample
#'
#' @description
#' `readOneProteom` reads and processes proteomics data for a single sample, applying filters for peptide count and optionally using LFQ quantification. It returns a data.table with useful columns and unique identifiers.
#'
#' @param inputTab A data.table or data.frame containing the input data for the proteomics sample.
#' @param sampleName A character string specifying the name of the sample to be processed.
#' @param pepNumCut A numeric value specifying the minimum number of peptides required for a protein to be included. Default is `1`.
#' @param ifLFQ A logical value indicating whether to use LFQ quantification. Default is `TRUE`.
#'
#' @return A data.table with the processed proteomics data, including columns for intensity, Uniprot ID, peptide counts, and gene names.
#'
#' @details
#' This function processes proteomics data for a single sample by filtering based on the number of peptides and optionally using LFQ quantification. It ensures that unique identifiers are created for each protein, and removes rows with missing or zero quantification values.
#'
#' @import data.table
#' @examples
#' inputTab <- data.table::data.table(
#'      Intensity.sample1 = c(100, 200, NA, 50),
#'      LFQ.intensity.sample1 = c(90, 180, NA, 45),
#'      Razor...unique.peptides.sample1 = c(2, 1, 0, 3),
#'      Protein.IDs = c("P1", "P2", "P3", "P4"),
#'      Peptide.counts..all. = c(2, 1, 0, 3),
#'      Gene.names = c("Gene1", "Gene2", "Gene3", "Gene4")
#' )
#' result <- readOneProteom(inputTab, sampleName = "sample1", pepNumCut = 1, ifLFQ = TRUE)
#'
#' @export
readOneProteom <- function(inputTab, sampleName, pepNumCut = 1, ifLFQ = TRUE) {

  # Define sample specific column names
  colSele <- paste0(c("Intensity.","LFQ.intensity.","Razor...unique.peptides."),
                    sampleName)
  # Peptide count filtering, based on Razor plus unique peptides
  keepRow <- inputTab[[colSele[3]]] >= pepNumCut

  # If no records pass the peptide count filter, return NULL with a warning
  if (all(!keepRow)) {
    warning(sprintf("sample %s does not contain any records after filtering",sampleName))
    return(NULL)
  }

  # Determine whether to use LFQ quantification
  if (ifLFQ) {
    quantCol <- colSele[2]
  } else {
    quantCol <- colSele[1]
  }

  # Extract and filter relevant columns
  outputTab <- inputTab[keepRow,c(quantCol,
                                  "Protein.IDs", "Peptide.counts..all.",
                                  "Gene.names"), with=FALSE]

  # Create a unique identifier for each row
  outputTab$rowName <- outputTab$Protein.IDs

  # Ensure that identifiers are unique
  stopifnot(all(!duplicated(outputTab$rowName)))
  rownames(outputTab) <- outputTab$rowName
  outputTab$rowName <- NULL

  # Remove rows with NA or zero quantification values
  outputTab <- outputTab[!is.na(outputTab[[quantCol]]) &
                           outputTab[[quantCol]]>0,]

  # Rename columns for clarity
  colnames(outputTab) <- c("Intensity","UniprotID","PeptideCounts","Gene")

  return(outputTab)
}


#' @name readProteomeExperiment
#'
#' @title Read and Process Proteomics Experiment Data
#'
#' @description
#' `readProteomeExperiment` reads and processes proteomics data from multiple samples, applying various quality filters, and returns a `SummarizedExperiment` object.
#'
#' @param fileTable A data.table or data.frame containing the file information with columns for file names, sample names, and IDs.
#' @param fdrCut A numeric value specifying the maximum false discovery rate (FDR) threshold. Default is `0.1`.
#' @param scoreCut A numeric value specifying the minimum score threshold. Default is `10`.
#' @param pepNumCut A numeric value specifying the minimum number of peptides required for a protein to be included. Default is `1`.
#' @param ifLFQ A logical value indicating whether to use LFQ quantification. Default is `TRUE`.
#'
#' @return A `SummarizedExperiment` object containing the processed proteomics data.
#'
#' @details
#' This function processes proteomics data by filtering based on FDR, score, and peptide count, and optionally using LFQ quantification. It aggregates the data from multiple samples and constructs a SummarizedExperiment object.
#'
#' @examples
#' # se <- readProteomeExperiment(fileTable)
#' 
#' @import data.table SummarizedExperiment
#' 
#' @export
readProteomeExperiment <- function(fileTable, fdrCut = 0.1, scoreCut = 10, pepNumCut = 1, ifLFQ = TRUE) {
  # Select proteomics entries
  fileTable <- fileTable[fileTable$type == "proteome",]
  if (nrow(fileTable) == 0) {
    return(NULL)
  }

  # Read in all batches and store them in a list
  expAll <- lapply(unique(fileTable$fileName), function(eachFileName) {

    fileTableSub <- fileTable[fileTable$fileName == eachFileName,]

    # Read input table with "\t" as delimiter
    inputTab <- data.table::fread(eachFileName, check.names = TRUE)
    # Remove unnecessary rows based on FDR and score thresholds
    inputTab <- inputTab[!inputTab$Protein.IDs %in% c(NA,"") &
                           #!inputTab$Gene.names %in% c("",NA) &
                           (!is.na(inputTab$Q.value) & inputTab$Q.value <= fdrCut) &
                           (!is.na(inputTab$Score) & inputTab$Score >= scoreCut),]

    # Stop if none of the proteins passed the chosen threshold
    if (nrow(inputTab) == 0) stop("No proteins could pass the specified threshold in any sample!")

    # Process each sample
    expSub <- lapply(seq(nrow(fileTableSub)), function(i) {
      eachTab <- readOneProteom(inputTab,
                                fileTableSub[i,]$sample,
                                pepNumCut,ifLFQ)

      if (!is.null(eachTab)) {
        eachTab$id <- fileTableSub[i,]$id
        eachTab$rowName <- rownames(eachTab)
      }

      eachTab

    })
    expSub <- data.table::rbindlist(expSub)
    expSub
  })

  expAll <- data.table::rbindlist(expAll)

  # Stop if none of the proteins passed the chosen threshold
  if (nrow(expAll) == 0) stop("No proteins could pass the specified threshold in any sample!")

  # Prepare annotations
  annoTab <- expAll[!duplicated(expAll$rowName),c("rowName", "UniprotID",
                                                  "Gene", "PeptideCounts")]
  rownames(annoTab) <- annoTab$rowName
  annoTab$rowName <- NULL

  # Prepare intensity matrix
  protMat <- matrix(data = rep(NA, nrow(annoTab)*nrow(fileTable)),
                    nrow(annoTab), nrow(fileTable),
                    dimnames = list(rownames(annoTab),fileTable$id))

  for (each_id in fileTable$id) {
    eachTab <- expAll[(expAll$id %in% each_id),]
    protMat[,each_id] <- as.numeric(eachTab[match(rownames(protMat),eachTab$rowName),][["Intensity"]])
  }

  # Rename rownames
  rownames(protMat) <- rownames(annoTab) <- paste0("p",seq(nrow(annoTab)))

  # Construct SummarizedExperiment object
  fpe <- SummarizedExperiment(assays = list(Intensity = protMat),
                              rowData = annoTab)
  return(fpe)
}


#' @name readOneProteomDIA
#'
#' @title Read and Process a Single DIA Proteomics Sample
#'
#' @description
#' `readOneProteomDIA` reads and processes data from a single DIA proteomics sample, applying filtering and data transformation steps.
#'
#' @param inputTab A data.table or data.frame containing the input data.
#' @param sampleName A character string specifying the sample name.
#'
#' @return A data.table containing the processed data with columns for intensity, UniProt ID, and gene name.
#'
#' @details
#' This function processes DIA proteomics data for a single sample by filtering out rows with non-quantitative data, converting character values to numeric, and renaming columns for consistency. It also ensures that each protein group has a unique identifier.
#'
#' @examples
#' # result <- readOneProteomDIA(inputTab, sampleName)
#' 
#' @import data.table
#' @export
readOneProteomDIA <- function(inputTab, sampleName) {

  sampleName <- make.names(sampleName)  # Make sample name syntactically valid

  # Define sample-specific column names
  if (length(grep(pattern = paste0("*", sampleName, ".*PG.Quantity"), colnames(inputTab))) > 0) {
    colSele <- colnames(inputTab)[grep(pattern = paste0("*", sampleName, ".*PG.Quantity"), colnames(inputTab))]
  } else {
    stop("Sample not found in quantification file")
  }

  # Replace "Filtered" values with NA
  inputTab[[colSele[1]]][inputTab[[colSele[1]]] == "Filtered"] <- NA

  # Convert character values to numeric
  inputTab[[colSele[1]]] <- as.numeric(gsub(",", ".", inputTab[[colSele[[1]]]])) #also change , to . if present

  # Remove NA or 0 quantification
  keepRow <- (!is.na(inputTab[[colSele[1]]]) & inputTab[[colSele[1]]]>0)

  # Check if any rows remain after filtering
  if (all(!keepRow)) {
    warning(sprintf("sample %s does not contain any records after filtering",sampleName))
    return(NULL)
  }

  # Output useful information
  outputTab <- inputTab[keepRow, c(colSele[1],"PG.ProteinGroups", "PG.Genes"), with=FALSE]

  # Create a unique identifier
  outputTab$rowName <- outputTab$PG.ProteinGroups
  # Ensure identifiers are unique
  stopifnot(all(!duplicated(outputTab$rowName)))
  rownames(outputTab) <- outputTab$rowName
  outputTab$rowName <- NULL

  # Rename columns for consistency
  colnames(outputTab) <- c("Intensity","UniprotID","Gene")

  return(outputTab)
}

#' @name readProteomeExperimentDIA
#'
#' @title Read and Process a DIA Proteome Experiment
#'
#' @description
#' `readProteomeExperimentDIA` reads and processes DIA (Data-Independent Acquisition) proteome data from multiple files and constructs a `SummarizedExperiment` object.
#'
#' @param fileTable A data frame containing metadata about the files to be read. Must contain columns `type`, `fileName`, `id`, and optionally `outputID`.
#' @param showProgressBar Logical, whether to show a progress bar during processing. Default is `FALSE`.
#'
#' @return A `SummarizedExperiment` object containing the processed proteome data.
#' #' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Filters the `fileTable` to include only rows where `type` is "proteome".
#'   \item For each file specified in `fileTable`, reads the data using `data.table::fread`.
#'   \item Removes rows where the `PG.ProteinGroups` column is NA or empty.
#'   \item Processes each sample in parallel using `BiocParallel::bplapply`, applying the `readOneProteomDIA` function to filter and clean the data for each sample.
#'   \item Combines the processed data from all files.
#'   \item Constructs a matrix of intensities with rows corresponding to proteins and columns corresponding to samples.
#'   \item Constructs a `SummarizedExperiment` object with the intensity matrix and protein annotations.
#' }
#' The `readOneProteomDIA` function is used to read and filter the data for each individual sample, and it must be available in the environment.
#'
#' @examples
#' # se <- readProteomeExperimentDIA(fileTable)
#' 
#' @import data.table
#' @import BiocParallel
#' @import SummarizedExperiment
#' @export
readProteomeExperimentDIA <- function(fileTable, showProgressBar = FALSE) {
  # Select proteomics entries
  fileTable <- fileTable[fileTable$type == "proteome",]
  if (nrow(fileTable) == 0) {
    return(NULL)
  }

  expAll <- lapply(unique(fileTable$fileName), function(eachFileName) {

    fileTableSub <- fileTable[fileTable$fileName == eachFileName,]

    # Read input table, "\t" as delimiter
    inputTab <- data.table::fread(eachFileName, check.names = TRUE)
    # Remove unnecessary rows
    inputTab <- inputTab[!inputTab$PG.ProteinGroups %in% c(NA,""),]

    # Process each sample in the file
    expSub <- BiocParallel::bplapply(seq(nrow(fileTableSub)), function(i) {
      eachTab <- readOneProteomDIA(inputTab,
                                   sampleName = fileTableSub[i,]$id)
      if(!is.null(eachTab)) {
        # Use user-specified output sample IDs if available
        if ("outputID" %in% colnames(fileTableSub)) {
          eachTab$id <- fileTableSub[i,]$outputID
        } else {
          eachTab$id <- fileTableSub[i,]$id
        }
        eachTab$rowName <- rownames(eachTab)
      }

      eachTab

    },BPPARAM = BiocParallel::MulticoreParam(progressbar = showProgressBar))
    expSub <- data.table::rbindlist(expSub)
    expSub
  })
  expAll <- data.table::rbindlist(expAll)

  # Stop if none of the proteins passed chosen threshold
  if (nrow(expAll) == 0) stop("No proteins could pass the specified threshold in any sample!")

  # Prepare annotations
  annoTab <- expAll[!duplicated(expAll$rowName),c("rowName","UniprotID",
                                                  "Gene")]
  rownames(annoTab) <- annoTab$rowName
  annoTab$rowName <- NULL

  # Prepare intensity matrix
  if ("outputID" %in% colnames(fileTable)) {
    sampleID <- fileTable$outputID
  } else {
    sampleID <- fileTable$id
  }
  protMat <- matrix(data = rep(NA, nrow(annoTab)*nrow(fileTable)),
                    nrow(annoTab), nrow(fileTable),
                    dimnames = list(rownames(annoTab),sampleID))
  for (each_id in sampleID) {
    eachTab <- expAll[(expAll$id %in% each_id),]
    protMat[,each_id] <- as.numeric(eachTab[match(rownames(protMat),eachTab$rowName),][["Intensity"]])
  }

  # Rename rownames
  rownames(protMat) <- rownames(annoTab) <- paste0("p",seq(nrow(annoTab)))

  # Construct SummarizedExperiment object
  fpe <- SummarizedExperiment(assays = list(Intensity = protMat),
                              rowData = annoTab)
  return(fpe)
}

#function to run matrixQCvis
runMatrixQCvis <- function(mae, index = 1) {
  if (index > length(mae)) {
    print("Index out of range of MultiAssayExperiment object")
  }
  else {
    se <- mae[[index]]
    colData(se) <- colData(mae)
    qc <- shinyQC(se)
  }
}

#function to deal with one peptide mapped to several proteins (not used currently)
dealMultiMap <- function(annoTab, method = "remove") {
  getElement <- function(vec, pos) {
    return(ifelse(pos == "first", vec[1], vec[length(vec)]))
  }
  #sanity check.
  stopifnot(method %in% c("remove","first","last","copy"))
  #get rows with multi-mapping
  multiTab <- annoTab[grepl(";",rownames(annoTab)),]
  uniqueTab <- annoTab[!grepl(";",rownames(annoTab)),]
  if (method == "remove") {
    return(uniqueTab)
  } else if (method == "first") {
    multiTab$Gene<- sapply(strsplit(multiTab$Gene, ";"), getElement, "first")
    multiTab$Position <- sapply(strsplit(multiTab$Position, ";"), getElement, "first")
    multiTab$Sequence <- sapply(strsplit(multiTab$Sequence, ";"), getElement, "first")
  } else if (method == "last") {
    multiTab$Gene<- sapply(strsplit(multiTab$Gene, ";"), getElement, "last")
    multiTab$Position <- sapply(strsplit(multiTab$Position, ";"), getElement, "last")
    multiTab$Sequence <- sapply(strsplit(multiTab$Sequence, ";"), getElement, "last")
  }
  fullTab <- bind_rows(uniqueTab, multiTab)
  fullTab <- fullTab[order(rownames(fullTab)),]
  return(fullTab)
}

#function to remove prefix or suffix of PP or FP samples
removePreSuffix <- function(sampleName) {
  sampleName <- gsub("([._](FullProteome|FP|Phospho|PP)$)|(^(FullProteome|FP|Phospho|PP)[._])","",sampleName)
  return(sampleName)
}
