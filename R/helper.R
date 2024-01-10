#Function to parse one phosphoproteome table
readOnePhos <- function(inputTab, sampleName, localProbCut, scoreDiffCut, multiMap) {


    #define sample specific column names
    colSele <- paste0(c("Localization.prob.","Score.diff.","Intensity."),
                      sampleName)
    if (!all(colSele %in% colnames(inputTab))) stop("Sample not found in quantification file")

    #get features passed quality filter and non zero intensity
    keepRow <- (!is.na(inputTab[[colSele[1]]]) & inputTab[[colSele[1]]] >= localProbCut) &
        (!is.na(inputTab[[colSele[2]]]) & inputTab[[colSele[2]]] >= scoreDiffCut) &
        (!is.na(inputTab[[colSele[3]]]) & inputTab[[colSele[3]]]>0)

    if (all(!keepRow)) {
        warning(sprintf("sample %s does not contain any records after filtering",sampleName))
        return(NULL)
    }

    #subset
    outputTab <- inputTab[keepRow,
                       c(colSele[3],"Proteins","Gene.names",
                         "Positions.within.proteins","Amino.acid","Sequence.window"),
                       with=FALSE]

    #create a uniqfied identifier
    outputTab$rowName <- paste0(outputTab$Proteins, "_",
                                outputTab$Positions.within.proteins)

    #it's important that identifiers are unique
    stopifnot(all(!duplicated(outputTab$rowName)))

    #output useful information
    rownames(outputTab) <- outputTab$rowName
    outputTab$rowName <- NULL
    #rename columns
    colnames(outputTab) <- c("Intensity","UniprotID","Gene","Position","Residue","Sequence")

    return(outputTab)
}

#Read the whole phosphoproteom and create a SummarizedExperiment object
readPhosphoExperiment <- function(fileTable, localProbCut, scoreDiffCut) {
    #select phosphoproteomic entries
    fileTable <- fileTable[fileTable$type == "phosphoproteome",]
    if (nrow(fileTable) == 0) {
      return(NULL)
    }

    #read in all batch and store them in a list
    expAll <- lapply(unique(fileTable$fileName), function(eachFileName) {
        fileTableSub <- fileTable[fileTable$fileName == eachFileName,]

        #read input table, "\t" as delimiter
        inputTab <- data.table::fread(eachFileName, sep = "\t", check.names = TRUE)
        #remove empty features
        inputTab <- inputTab[!inputTab$Proteins %in% c(NA,"") &
                                 #!inputTab$Gene.names %in% c("",NA) &
                                 !inputTab$Positions.within.proteins %in% c(NA,""),]
        #remove reverse and potential contaminants for the whole table
        inputTab <- inputTab[!inputTab$Potential.contaminant %in% "+" &
                                 !inputTab$Reverse %in% "+",]

        # stop if none of the record passed chosen threshold
        if (nrow(inputTab) == 0) stop("No phosphorylation site could pass the specified threshold in any sample!")

        #each data for each sample
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

    expAll <- data.table::rbindlist(expAll) #rbindlist is faster than do.call(rbind)

    # stop if none of the record passed chosen threshold
    if (nrow(expAll) == 0) stop("No phosphorylation site could pass the specified threshold in any sample!")

    #prepare annotations
    annoTab <- expAll[!duplicated(expAll$rowName),c("rowName","UniprotID",
                                                    "Gene","Position","Residue","Sequence")]
    annoTab$site <- annoTab$rowName
    rownames(annoTab) <- annoTab$rowName
    annoTab$rowName <- NULL

    #prepare intensity matrix
    phosMat <- matrix(data = rep(NA, nrow(annoTab)*nrow(fileTable)),
                      nrow(annoTab), nrow(fileTable),
                      dimnames = list(rownames(annoTab),fileTable$id))
    for (each_id in fileTable$id) {
        eachTab <- expAll[(expAll$id %in% each_id),]
        phosMat[,each_id] <- as.numeric(eachTab[match(rownames(phosMat),eachTab$rowName),][["Intensity"]])
    }

    #rename rownames
    rownames(phosMat) <- rownames(annoTab) <- paste0("s",seq(nrow(annoTab)))

    #construct SummariseExperiment object object
    ppe <- SummarizedExperiment(assays = list(Intensity = phosMat),
                                rowData = annoTab)
    return(ppe)
}



#-----------------------------------------------------------------------------------------------------------------

#Function to parse one phosphoproteome table (DIA)
readOnePhosDIA <- function(inputTab, sampleName, localProbCut, removeDup = FALSE) {

    #define sample specific column names
    colSele <- c(NA,NA) # specify the order of localisation probability and quantity, in case the order changes in the output file
    colSele[1] <- colnames(inputTab)[grep(pattern = paste0("*", sampleName, ".*PTM.SiteProbability"), colnames(inputTab))]
    colSele[2] <- colnames(inputTab)[grep(pattern = paste0("*", sampleName, ".*PTM.Quantity"), colnames(inputTab))]

    #replace "filtered" values with NA
    inputTab[[colSele[1]]][inputTab[[colSele[1]]] == "Filtered"] <- NA
    inputTab[[colSele[2]]][inputTab[[colSele[2]]] == "Filtered"] <- NA

    #change to numeric values
    inputTab[[colSele[1]]] <- as.numeric(gsub(",", ".", inputTab[[colSele[1]]]))
    inputTab[[colSele[2]]] <- as.numeric(gsub(",", ".", inputTab[[colSele[2]]])) #change , to . if any. Sometimes the "," is used as decimal.

    if (!all(colSele %in% colnames(inputTab))) stop("Sample not found in quantification file")

    #get features passed quality filter and non zero intensity
    keepRow <- (!is.na(inputTab[[colSele[1]]]) & inputTab[[colSele[1]]] >= localProbCut) &
        (!is.na(inputTab[[colSele[2]]]) & inputTab[[colSele[2]]]>0)

    if (all(!keepRow)) {
        warning(sprintf("sample %s does not contain any records after filtering", sampleName))
        return(NULL)
    }


    #subset
    outputTab <- inputTab[keepRow,
                          c(colSele[2],"PTM.CollapseKey","PG.UniProtIds","PG.Genes","PTM.Multiplicity",
                            "PTM.SiteLocation","PTM.SiteAA","PTM.FlankingRegion"), with=FALSE]
    #rename column names
    colnames(outputTab) <- c("Intensity","CollapseKey","UniprotID","Gene","Multiplicity","Position","Residue","Sequence")


    #deal with multiplicity
    #outputTab$CollapseKeyNew <- substring(outputTab$CollapseKey, 1, nchar(outputTab$CollapseKey)-1)
    outputTab$CollapseKeyNew <-  gsub("_M\\d+","",outputTab$CollapseKey) #it's safer to use regular expression to remove the suffix _M together with the numbers.

    # Summarise multiplicity
    outputTab <- outputTab[order(outputTab$Multiplicity, decreasing = TRUE),]
    outputTabNew <- outputTab[!duplicated(outputTab$CollapseKeyNew),]
    outputTabNew$CollapseKey <- NULL
    intensityTab <- aggregate(Intensity ~ CollapseKeyNew, outputTab, sum)
    outputTabNew$Intensity <- intensityTab[match(outputTabNew$CollapseKeyNew, intensityTab$CollapseKeyNew),]$Intensity
    outputTab <- outputTabNew

    if (removeDup) {
        #remove duplicates
        outputTab.rev <- outputTab[rev(seq(nrow(outputTab))),]
        outputTab.rev <- outputTab.rev[!duplicated(paste0(outputTab.rev$UniprotID,"_",
                                                          outputTab.rev[[colSele[[2]]]])),]
        outputTab <- outputTab.rev[rev(seq(nrow(outputTab.rev))),]

        #create a uniqfied identifier for rownames
        outputTab$rowName <- paste0(outputTab$UniprotID, "_",
                                    outputTab$Position)
    } else {
        outputTab$rowName <- outputTab$CollapseKeyNew
    }

    #delete the PTM.CollapseKeyNew column
    outputTab$CollapseKeyNew <- NULL

    #it's important that identifiers are unique
    stopifnot(all(!duplicated(outputTab$rowName)))

    #output useful information
    rownames(outputTab) <- outputTab$rowName
    outputTab$rowName <- NULL

    return(outputTab)
}


#Read the whole phosphoproteom and create a SummarizedExperiment object (DIA)
readPhosphoExperimentDIA <- function(fileTable, localProbCut, onlyReviewed = TRUE,
                                     showProgressBar = FALSE) {
    #select phosphoproteomic entries
    fileTable <- fileTable[fileTable$type == "phosphoproteome",]
    if (nrow(fileTable) == 0) {
      return(NULL)
    }

    #read in all batch and store them in a list
    expAll <- lapply(unique(fileTable$fileName), function(eachFileName) {
        fileTableSub <- fileTable[fileTable$fileName == eachFileName,]

        #read input table, "\t" as delimiter, fread is faster than read.delim
        inputTab <- data.table::fread(eachFileName, sep = "\t", check.names = TRUE)
        #keep only Phospho (STY) modifications
        inputTab <- inputTab[inputTab$PTM.ModificationTitle == "Phospho (STY)",]
        #remove empty features
        inputTab <- inputTab[!inputTab$PG.UniProtIds %in% c(NA,"") &
                                 #!inputTab$Gene.names %in% c("",NA) &
                                 !inputTab$PTM.SiteLocation %in% c(NA,""),]

        # stop if none of the record passed chosen threshold
        if (nrow(inputTab) == 0) stop("No phosphorylation site could pass the specified threshold in any sample!")

        #if only reviewed protiens are considered
        if (onlyReviewed) {
            data("swissProt")
            inputTab <- inputTab[inputTab$PTM.ProteinId %in% swissProt$Entry,]
        }

        #each data for each sample
        expSub <- BiocParallel::bplapply(seq(nrow(fileTableSub)), function(i) {
            eachTab <- readOnePhosDIA(inputTab = inputTab,
                                   sampleName = fileTableSub[i,]$id,
                                   localProbCut = localProbCut)

            if (!is.null(eachTab)) {
                if ("outputID" %in% colnames(fileTableSub)) { #use user-specified output sample IDs
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

    # stop if none of the record passed chosen threshold
    if (nrow(expAll) == 0) stop("No phosphorylation site could pass the specified threshold in any sample!")

    expAll <- expAll[order(expAll$rowName),]

    #prepare annotations
    annoTab <- expAll[!duplicated(expAll$rowName),c("rowName","UniprotID",
                                                    "Gene","Multiplicity","Position","Residue","Sequence")]
    annoTab$site <- annoTab$rowName
    rownames(annoTab) <- annoTab$rowName
    annoTab$rowName <- NULL

    #prepare intensity matrix
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

    #rename rownames
    rownames(phosMat) <- rownames(annoTab) <- paste0("s",seq(nrow(annoTab)))

    #construct SummariseExperiment object object
    ppe <- SummarizedExperiment(assays = list(Intensity = phosMat),
                                rowData = annoTab)
    return(ppe)
}

#------------------------------------------------------------------------------------------------------------------

#Read one proteome assay
readOneProteom <- function(inputTab, sampleName, pepNumCut, ifLFQ) {

    #define sample specific column names
    colSele <- paste0(c("Intensity.","LFQ.intensity.","Razor...unique.peptides."),
                      sampleName)
    #peptide count filtering, based on Razor plus unique peptide
    keepRow <- inputTab[[colSele[3]]] >= pepNumCut

    if (all(!keepRow)) {
        warning(sprintf("sample %s does not contain any records after filtering",sampleName))
        return(NULL)
    }

    #whether use LFQ quantification, recommended
    if (ifLFQ) quantCol <- colSele[2] else quantCol <- colSele[1]

    #output useful information
    outputTab <- inputTab[keepRow,c(quantCol,
                             "Protein.IDs", "Peptide.counts..all.",
                             "Gene.names"), with=FALSE]

    #create a uniqfied identifier
    outputTab$rowName <- outputTab$Protein.IDs
    #it's important that identifiers are unique
    stopifnot(all(!duplicated(outputTab$rowName)))
    rownames(outputTab) <- outputTab$rowName
    outputTab$rowName <- NULL

    #remove NA or 0 quantification
    outputTab <- outputTab[!is.na(outputTab[[quantCol]]) &
                               outputTab[[quantCol]]>0,]

    #rename columns
    colnames(outputTab) <- c("Intensity","UniprotID","PeptideCounts","Gene")

    return(outputTab)
}

#Read the whole full proteome and create a SummarizedExperiment object
readProteomeExperiment <- function(fileTable, fdrCut, scoreCut, pepNumCut, ifLFQ) {
    #select proteomics entries
    fileTable <- fileTable[fileTable$type == "proteome",]
    if (nrow(fileTable) == 0) {
      return(NULL)
    }

    expAll <- lapply(unique(fileTable$fileName), function(eachFileName) {

        fileTableSub <- fileTable[fileTable$fileName == eachFileName,]

        #read input table, "\t" as delimiter
        inputTab <- data.table::fread(eachFileName, sep = "\t", check.names = TRUE)
        #remove unnecessary rows
        inputTab <- inputTab[!inputTab$Protein.IDs %in% c(NA,"") &
                                 #!inputTab$Gene.names %in% c("",NA) &
                                 (!is.na(inputTab$Q.value) & inputTab$Q.value <= fdrCut) &
                                 (!is.na(inputTab$Score) & inputTab$Score >= scoreCut),]

        # stop if none of the proteins passed chosen threshold
        if (nrow(inputTab) == 0) stop("No proteins could pass the specified threshold in any sample!")

        #each data for each sample
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

    # stop if none of the proteins passed chosen threshold
    if (nrow(expAll) == 0) stop("No proteins could pass the specified threshold in any sample!")

    #prepare annotations
    annoTab <- expAll[!duplicated(expAll$rowName),c("rowName", "UniprotID",
                                                    "Gene", "PeptideCounts")]
    rownames(annoTab) <- annoTab$rowName
    annoTab$rowName <- NULL

    #prepare intensity matrix
    protMat <- matrix(data = rep(NA, nrow(annoTab)*nrow(fileTable)),
                      nrow(annoTab), nrow(fileTable),
                      dimnames = list(rownames(annoTab),fileTable$id))

    for (each_id in fileTable$id) {
        eachTab <- expAll[(expAll$id %in% each_id),]
        protMat[,each_id] <- as.numeric(eachTab[match(rownames(protMat),eachTab$rowName),][["Intensity"]])
    }

    #rename rownames
    rownames(protMat) <- rownames(annoTab) <- paste0("p",seq(nrow(annoTab)))

    #construct SummariseExperiment  object
    fpe <- SummarizedExperiment(assays = list(Intensity = protMat),
                                rowData = annoTab)
    return(fpe)
}

#Read one proteome assay (DIA)
readOneProteomDIA <- function(inputTab, sampleName) {
    #define sample specific column names
    colSele <- colnames(inputTab)[grep(pattern = paste0("*", sampleName, ".*PG.Quantity"), colnames(inputTab))]

    #replace "filtered" values with NA
    inputTab[[colSele[1]]][inputTab[[colSele[1]]] == "Filtered"] <- NA

    #convert character values to numeric
    inputTab[[colSele[1]]] <- as.numeric(gsub(",", ".", inputTab[[colSele[[1]]]])) #also change , to . if present

    #remove NA or 0 quantification
    keepRow <- (!is.na(inputTab[[colSele[1]]]) & inputTab[[colSele[1]]]>0)

    if (all(!keepRow)) {
        warning(sprintf("sample %s does not contain any records after filtering",sampleName))
        return(NULL)
    }

    #output useful information
    outputTab <- inputTab[keepRow, c(colSele[1],"PG.ProteinGroups", "PG.Genes"), with=FALSE]

    #create a uniqfied identifier
    outputTab$rowName <- outputTab$PG.ProteinGroups
    #it's important that identifiers are unique
    stopifnot(all(!duplicated(outputTab$rowName)))
    rownames(outputTab) <- outputTab$rowName
    outputTab$rowName <- NULL

    #rename columns
    colnames(outputTab) <- c("Intensity","UniprotID","Gene")

    return(outputTab)
}

#Read the whole full proteome and create a SummarizedExperiment object (DIA)
readProteomeExperimentDIA <- function(fileTable, showProgressBar = FALSE) {
    #select proteomics entries
    fileTable <- fileTable[fileTable$type == "proteome",]
    if (nrow(fileTable) == 0) {
      return(NULL)
    }

    expAll <- lapply(unique(fileTable$fileName), function(eachFileName) {

        fileTableSub <- fileTable[fileTable$fileName == eachFileName,]

        #read input table, "\t" as delimiter
        inputTab <- data.table::fread(eachFileName, sep = "\t", check.names = TRUE)
        #remove unnecessary rows
        inputTab <- inputTab[!inputTab$PG.ProteinGroups %in% c(NA,""),]

        #each data for each sample
        expSub <- BiocParallel::bplapply(seq(nrow(fileTableSub)), function(i) {
            eachTab <- readOneProteomDIA(inputTab,
                                      sampleName = fileTableSub[i,]$id)
            if(!is.null(eachTab)) {
                if ("outputID" %in% colnames(fileTableSub)) { #use user-specified output sample IDs
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

    # stop if none of the proteins passed chosen threshold
    if (nrow(expAll) == 0) stop("No proteins could pass the specified threshold in any sample!")

    #prepare annotations
    annoTab <- expAll[!duplicated(expAll$rowName),c("rowName","UniprotID",
                                                    "Gene")]
    rownames(annoTab) <- annoTab$rowName
    annoTab$rowName <- NULL

    #prepare intensity matrix
    if ("outputID" %in% colnames(fileTable)) { #use user specific sample ID
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

    #rename rownames
    rownames(protMat) <- rownames(annoTab) <- paste0("p",seq(nrow(annoTab)))

    #construct SummariseExperiment  object
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
