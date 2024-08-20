#' @name mscale
#' 
#' @title Scale and Center a Matrix
#'
#' @description
#' `mscale` scales and centers each row of a matrix, with options for using mean or median, standard deviation or mean absolute deviation, and censoring extreme values.
#'
#' @param x A numeric matrix where rows are features and columns are samples.
#' @param center Logical. If TRUE, the rows are centered by subtracting the mean or median. Default is `TRUE`.
#' @param scale Logical. If TRUE, the rows are scaled by dividing by the standard deviation or mean absolute deviation. Default is `TRUE`.
#' @param censor A numeric vector of length one or two for censoring the scaled values. 
#' If length one, values are censored symmetrically at positive and negative values. 
#' If length two, the first value is the lower limit and the second value is the upper limit. Default is `NULL`.
#' @param useMad Logical. If TRUE, the mean absolute deviation is used for scaling instead of the standard deviation. Default is `FALSE`.
#'
#' @return A scaled and centered numeric matrix with the same dimensions as the input matrix `x`.
#'
#' @details
#' The function allows for flexible scaling and centering of the rows of a matrix:
#' \itemize{
#'   \item If both `center` and `scale` are TRUE, rows are centered and scaled.
#'   \item If only `center` is TRUE, rows are centered but not scaled.
#'   \item If only `scale` is TRUE, rows are scaled but not centered.
#'   \item If neither `center` nor `scale` is TRUE, the original matrix is returned.
#' }
#' The function can also censor extreme values, either symmetrically or asymmetrically, based on the `censor` parameter.
#'
#' @examples
#' # Example usage:
#' # Assuming `dataMatrix` is a numeric matrix with expression data
#' # scaledMatrix <- mscale(dataMatrix, center = TRUE, scale = TRUE, censor = 2)
#'
#' @importFrom stats median sd
#' @importFrom matrixStats rowMads
mscale <- function(x, center = TRUE, scale = TRUE, censor = NULL, useMad = FALSE){
  
  # Check if both scaling and centering are requested
  if (scale & center) {
    # Scale using Mean Absolute Deviation (MAD)
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) (y-median(y,na.rm = T))/meanAD(y))
    } else {
      # Scale using Standard Deviation (SD)
      x.scaled <- apply(x, 1, function(y) (y-mean(y,na.rm=T))/sd(y,na.rm = T))
    }
  } else if (center & !scale) {
    # Only center the data
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) (y-median(y,na.rm=T)))
    } else {
      x.scaled <- apply(x, 1, function(y) (y-mean(y,na.rm=T)))
    }
  } else if (!center & scale) {
    # Only scale the data
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) y/meanAD(y))
    } else {
      x.scaled <- apply(x, 1, function(y) y/sd(y,na.rm = T))
    }
  } else {
    # Neither center nor scale
    x.scaled <- t(x)
  }
  
  # Apply censoring if requested
  if (!is.null(censor)) {
    if (length(censor) == 1) {
      # Symmetric censoring
      x.scaled[x.scaled > censor] <- censor
      x.scaled[x.scaled < -censor] <- -censor
    } else {
      # Asymmetric censoring
      x.scaled[x.scaled > censor[2]] <- censor[2]  # Upper limit
      x.scaled[x.scaled < censor[1]] <- censor[1]  # Lower limit
    }
  }
  return(t(as.matrix(x.scaled)))
}


#' @name addZeroTime
#' 
#' @title Add Zero Timepoint Data to Treatment Subset
#'
#' @description
#' `addZeroTime` adds a zero timepoint to a specific treatment's data subset.
#'
#' @param data A SummarizedExperiment object containing the experimental data.
#' @param condition The condition corresponds to one of the columns from the colData of SE object.
#' @param treat Character string specifying the treatment to which zero timepoint should be added.
#' @param zeroTreat Character string specifying the treatment representing the zero timepoint.
#' @param timeRange Character vector specifying the timepoints to include for the treatment.
#'
#' @return A SummarizedExperiment object with the zero timepoint added to the specified treatment's data.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Subsets the data for the specified treatment and time range.
#'   \item Subsets the data for the zero timepoint of the specified zero treatment.
#'   \item Combines the assays from the treatment and zero timepoint subsets.
#'   \item Updates the column data to reflect the combined treatment.
#'   \item Returns a SummarizedExperiment object with the combined data.
#' }
#'
#' @examples
#' # Example usage:
#' # Assuming data is a SummarizedExperiment object with appropriate data
#' # result <- addZeroTime(data, treat = "TreatmentA", zeroTreat = "Control", timeRange = c("10min", "20min"))
#'
#' @importFrom SummarizedExperiment colData rowData assay assays elementMetadata SummarizedExperiment
addZeroTime <- function(data, condition, treat, zeroTreat, timeRange) {
  # Subset the data for the specified treatment and time range
  subset1 <- data[, data[[condition]] == treat & data$timepoint %in% timeRange]
  # Subset the data for the zero timepoint of the specified zero treatment
  subset2 <- data[, data[[condition]] == zeroTreat & data$timepoint %in% c("0min", "0", "0h")]
  # Combine the assays from the treatment and zero timepoint subsets
  assay <- cbind(assay(subset1), assay(subset2))
  colnames(assay) <- gsub(zeroTreat, treat, colnames(assay))
  
  # Combine the column data from both subsets
  cd1 <- colData(subset1)
  cd2 <- colData(subset2)
  cd <- rbind(cd1, cd2)
  cd[[condition]][cd[[condition]] == zeroTreat] = treat
  rownames(cd) <- gsub(zeroTreat, treat, rownames(cd))
  
  # Retrieve the element metadata from the original data
  emeta <- elementMetadata(data)
  
  return(SummarizedExperiment(assays=SimpleList(intensity=assay), colData = cd, rowData = emeta))
}


#' @name clusterTS
#' 
#' @title Perform Clustering on Time-Series Data
#'
#' @description
#' `clusterTS` performs clustering on time-series data and generates plots for visualization.
#'
#' @param x A numeric matrix with rows as features and columns as time points.
#' @param k An integer specifying the number of clusters.
#' @param pCut A numeric value specifying the probability cutoff for cluster membership. Default is `NULL`.
#' @param twoCondition A logical value indicating if the data contains two conditions. Default is `FALSE`.
#'
#' @return A list containing:
#' \item{cluster}{A tibble with clustering information for each feature.}
#' \item{plot}{A ggplot object for visualizing the clustering results.}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Sets a seed for reproducibility.
#'   \item Removes rows with missing values.
#'   \item Performs clustering using fuzzy C-means.
#'   \item Filters clusters based on the probability cutoff if provided.
#'   \item Generates plots for visualizing clustering results.
#' }
#'
#' @examples
#' # Example usage:
#' # Assuming x is a numeric matrix with time-series data
#' # result <- clusterTS(x, k = 4, pCut = 0.8, twoCondition = FALSE)
#'
#' @importFrom e1071 cmeans
#' @importFrom dplyr filter left_join mutate arrange group_by ungroup
#' @importFrom tidyr pivot_longer separate
#' @importFrom ggplot2 ggplot aes geom_line scale_color_gradient facet_wrap theme_bw theme element_text unit
#' @importFrom tibble as_tibble tibble
#' @importFrom stringr str_extract str_split
#' @importFrom magrittr %>%
#' @export
clusterTS <- function(x, k, pCut = NULL, twoCondition = FALSE) {
  
  # Set seed for reproducible clustering results
  set.seed(12345)
  
  # Remove rows with NA values
  x.center <- x[complete.cases(x),]
  
  # Perform fuzzy C-means clustering
  res <- e1071::cmeans(x.center, k)
  
  # Create a tibble with clustering results
  resCluster <- tibble(feature = names(res$cluster),
                       cluster = res$cluster,
                       prob = rowMax(res$membership))
  
  # Filter clusters based on probability cutoff if provided
  if (!is.null(pCut)) resCluster <- filter(resCluster, prob >= pCut)
  
  if (!twoCondition) {
    # Handle single condition data
    
    # Extract unique time points and determine time unit
    timeVector <- unique(colnames(x.center))
    timeUnit <- str_extract(timeVector, "h|min")
    timeUnit <- ifelse(is.na(timeUnit), "", timeUnit)
    
    # Adjust time values if both hours and minutes are present
    if ((any(timeUnit == "h")) & (any(timeUnit == "min"))) {
      timeValue <- timeVector
      timeValue[timeUnit == "min"] <- 1/60 * as.numeric(gsub("min", "", timeValue[timeUnit == "min"]))
      timeRank <- rank(as.numeric(gsub("h", "", timeValue)))
    } else {
      timeRank <- rank(as.numeric(gsub("h|min", "", timeVector)))
    }
    
    # Determine the order of time points for plotting
    timeOrder <- timeVector[order(match(timeRank, sort(timeRank)))]  
    
    # Prepare data for plotting
    clusterTab <- x.center %>% as_tibble(rownames = "feature") %>%
      pivot_longer(-feature, names_to = "time", values_to = "value") %>%
      left_join(resCluster, by = "feature") %>% filter(!is.na(cluster)) %>%
      arrange(prob) %>%
      mutate(cluster = paste0("cluster",cluster),
             feature = factor(feature, levels = unique(feature))) %>%
      group_by(cluster) %>%
      mutate(cNum = length(unique(feature))) %>%
      mutate(clusterNum = sprintf("%s (%s)",cluster, cNum)) %>%
      ungroup() %>%
      mutate(time = factor(time, levels = timeOrder)) %>%
      mutate(time = droplevels(time))
    
    # Generate the plot
    p <- ggplot(clusterTab, aes( x = time, y = value, group = feature)) +
      geom_line( aes(col = prob), alpha=0.8) +
      scale_color_gradient(low = "navy", high = "green") +
      facet_wrap(~clusterNum,ncol=3, scales = "free") +
      theme_bw() +
      theme(legend.position = "bottom",
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 12),
            legend.key.width = unit(3,"cm"),
            legend.key.height = unit(0.8, "cm"),
            legend.title = element_text(size = 15, face = "bold"),
            legend.text = element_text(size = 15),
            strip.text = element_text(size=15, face = "bold"))
  } else {
    # Handle data with two conditions
    
    # Extract unique time points and determine time unit
    timeVector <- sapply(colnames(x.center), function(X) unlist(str_split(X, "_"))[1])
    timeVector <- unique(timeVector)
    timeUnit <- str_extract(timeVector, "h|min")
    timeUnit <- ifelse(is.na(timeUnit), "", timeUnit)
    
    # Adjust time values if both hours and minutes are present
    if ((any(timeUnit == "h")) & (any(timeUnit == "min"))) {
      timeValue <- timeVector
      timeValue[timeUnit == "min"] <- 1/60 * as.numeric(gsub("min", "", timeValue[timeUnit == "min"]))
      timeRank <- rank(as.numeric(gsub("h", "", timeValue)))
    } else {
      timeRank <- rank(as.numeric(gsub("h|min", "", timeVector)))
    }
    
    # Determine the order of time points for plotting
    timeOrder <- timeVector[order(match(timeRank, sort(timeRank)))]
    
    # Prepare data for plotting
    clusterTab <- x.center %>% as_tibble(rownames = "feature") %>%
      pivot_longer(-feature, names_to = "timeTreat", values_to = "value") %>%
      left_join(resCluster, by = "feature") %>% filter(!is.na(cluster)) %>%
      arrange(prob) %>%
      mutate(cluster = paste0("cluster",cluster),
             feature = factor(feature, levels = unique(feature))) %>%
      group_by(cluster) %>%
      mutate(cNum = length(unique(feature))) %>%
      mutate(clusterNum = sprintf("%s (%s)",cluster, cNum)) %>%
      ungroup() %>%
      separate(timeTreat,c("time","treatment"),"_", extra = "merge") %>%
      mutate(time = factor(time, levels = timeOrder)) %>%
      mutate(time = droplevels(time)) %>%
      mutate(geneGroup = paste0(feature,treatment))
    
    # Generate the plot
    p <- ggplot(clusterTab, aes( x=time, y= value, group = geneGroup)) +
      geom_line( aes(alpha = prob, color= treatment)) +
      facet_wrap(~clusterNum,ncol=3, scales = "free") +
      theme_bw() +
      theme(legend.position = "bottom",
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 12),
            legend.key.size = unit(1,"cm"),
            legend.title = element_text(size = 15, face = "bold"),
            legend.text = element_text(size = 15),
            strip.text = element_text(size=15, face = "bold"))
  }
  
  return(list(cluster = clusterTab, plot = p))
}


#' @name splineFilter
#' 
#' @title Filter Expression Matrix Using Spline Models
#'
#' @description
#' `splineFilter` filters an expression matrix based on spline models fitted to time-series data, optionally considering treatment and subject ID.
#'
#' @param exprMat A numeric matrix of expression data, where rows are features and columns are samples.
#' @param subjectID An optional vector of subject IDs corresponding to columns in `exprMat`. Default is `NULL`.
#' @param time A numeric vector representing the time points corresponding to columns in `exprMat`.
#' @param df An integer specifying the degrees of freedom for the spline basis.
#' @param pCut A numeric value for the p-value cutoff to filter significant features. Default is `0.05`.
#' @param ifFDR A logical value indicating if the false discovery rate (FDR) should be used for filtering. If FALSE, raw p-values are used. Default is `FALSE`.
#' @param treatment An optional vector of treatment labels corresponding to columns in `exprMat`. Default is `NULL`.
#' @param refTreatment An optional reference treatment label for the `treatment` vector. Default is `NULL`.
#'
#' @return A filtered expression matrix containing only the features that meet the significance criteria.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Converts time points from minutes to hours if both units are present.
#'   \item Removes rows with missing values from the expression matrix.
#'   \item Constructs a design matrix for the spline model, optionally including subject IDs and treatments.
#'   \item Fits a linear model using the design matrix and performs empirical Bayes moderation.
#'   \item Extracts significant features based on the specified p-value or FDR cutoff.
#' }
#'
#' @examples
#' # Example usage:
#' # Assuming exprMat is a numeric matrix with expression data
#' # filteredMat <- splineFilter(exprMat, time = timeVec, df = 4, pCut = 0.05, ifFDR = TRUE)
#'
#' @importFrom limma lmFit eBayes topTable
#' @importFrom splines ns
#' @importFrom dplyr filter as_tibble
#' @importFrom stringr str_ends
#' @importFrom tibble rownames_to_column
splineFilter <- function(exprMat, subjectID = NULL, time, df, pCut = 0.5, ifFDR = FALSE, treatment = NULL, refTreatment = NULL) {
  # Convert time points from minutes to hours if both units are present
  if ((all(str_ends(time,"h|min"))) & (!all(str_ends(time,"h"))) & (!all(str_ends(time,"min")))) {
    time[str_ends(time, "min")] <- 1/60 * as.numeric(gsub("min","", time[str_ends(time, "min")]))
  }
  time <- as.numeric(gsub("h|min", "", time))
  
  # Remove rows with NA values from the expression matrix
  exprMat <- exprMat[complete.cases(exprMat), ]
  if (is.null(treatment)) {
    # Handle case with one condition or log fold-change
    if (is.null(subjectID)) {
      # Create design matrix without subject IDs
      designTab <- data.frame(row.names = colnames(exprMat))
      designTab$X <- splines::ns(time, df)
      design <- model.matrix(~ 0 + X, data = designTab)
    } else { 
      # Create design matrix with subject IDs
      designTab <- data.frame(row.names = colnames(exprMat), subjectID = subjectID)
      designTab$X <- splines::ns(time, df)
      design <- model.matrix(~ 0 + X + subjectID, data = designTab)
    }
    
    # Fit linear model and perform empirical Bayes moderation
    fit <- lmFit(exprMat, design = design)
    fit2 <- eBayes(fit)
    # Extract results table and filter by p-value or FDR
    resTab <- topTable(fit2, coef = seq(df), number = Inf) %>%
      as_tibble(rownames = "ID") 
    
    if (ifFDR) resTab$p <- resTab$adj.P.Val else resTab$p <- resTab$P.Value
    
    resTab <- filter(resTab, p <= pCut)
    # Return filtered expression matrix
    
    return(exprMat[resTab$ID,])
  } else {
    # Handle case with two conditions
    if (is.null(subjectID)) {
      # Create design matrix without subject IDs
      designTab <- data.frame(row.names = colnames(exprMat), treatment = treatment)
      designTab$treatment <- factor(designTab$treatment, levels = unique(treatment))
      designTab$treatment <- relevel(designTab$treatment, refTreatment)
      designTab$treatment <- droplevels(designTab$treatment)
      designTab$X <- splines::ns(time, df)
      design <- model.matrix(~ 0 + X*treatment, data = designTab)
    } else {
      # Create design matrix with subject IDs
      designTab <- data.frame(row.names = colnames(exprMat), subjectID = subjectID, treatment = treatment)
      designTab$treatment <- factor(designTab$treatment, levels = unique(treatment))
      designTab$treatment <- relevel(designTab$treatment, refTreatment)
      designTab$treatment <- droplevels(designTab$treatment)
      designTab$X <- splines::ns(time, df)
      design <- model.matrix(~ 0 + subjectID + X*treatment, data = designTab)
    }
    
    # Fit linear model and perform empirical Bayes moderation
    fit <- lmFit(exprMat, design = design)
    fit2 <- eBayes(fit)
    
    # Extract results table and filter by p-value or FDR
    resTab <- topTable(fit2, coef = (ncol(design)-df+1):ncol(design), number = Inf) %>%
      as_tibble(rownames = "ID")
    
    if (ifFDR) resTab$p <- resTab$adj.P.Val else resTab$p <- resTab$P.Value
    
    resTab <- filter(resTab, p <= pCut)
    
    # Return filtered expression matrix
    return(exprMat[resTab$ID,])
  }
}


#' @name plotTimeSeries
#' 
#' @title Plot Time Series Data for a gene or phospho site from SummarizedExperiment object
#'
#' @description
#' `plotTimeSeries` plots time series data for a given gene or phospho site from a given SummarizedExperiment object, allowing different types of plots such as expression, log fold change, or two-condition expression.
#'
#' @param se A SummarizedExperiment object containing the data.
#' @param type The type of plot to generate. Options are `"expression"`, `"logFC"`, or `"two-condition expression"`.
#' @param geneID The identifier of the gene or feature to plot.
#' @param symbol The symbol or name of the gene or feature to use as the plot title.
#' @param condition The condition corresponds to one of the columns from the colData of SE object.
#' @param treatment The treatment to use for filtering the data.
#' @param refTreat The reference treatment to compare against.
#' @param addZero Logical, whether to add a zero time point to the data. Default is `FALSE`.
#' @param zeroTreat The treatment to use for adding the zero time point. Default is `NULL`.
#' @param timerange The range of time points to include in the plot.
#'
#' @return A ggplot object representing the time series plot.
#'
#' @details
#' This function generates time series plots for a specified gene or feature from a SummarizedExperiment (SE) object. The type of plot can be one of the following:
#' - "expression": Plots normalized expression levels over time.
#' - "logFC": Plots log fold change (logFC) over time, comparing a treatment to a reference treatment.
#' - "two-condition expression": Plots normalized expression levels over time for two conditions.
#'
#' The function can add a zero time point if specified and handles data with and without subject-specific information. The plot includes points for each time point and a summary line representing the mean value.
#'
#' The x-axis represents time, and the y-axis represents the selected metric (normalized expression or logFC). The plot is customized with various aesthetic elements, such as point size, line type, axis labels, and title formatting.
#'
#' @importFrom ggplot2 ggplot aes geom_point stat_summary geom_line ylab xlab ggtitle theme_bw theme element_text
#' @importFrom SummarizedExperiment assays
#' @importFrom dplyr %>% bind_cols filter
#' @importFrom stringr str_extract
#' @examples
#' # Assuming 'se' is a SummarizedExperiment object with appropriate data:
#' plotTimeSeries(se, type = "expression", geneID = "gene1", symbol = "Gene 1", condition = "treatment", treatment = "A", refTreat = "B", timerange = c("0h", "1h", "2h"))
#'
#' @export
plotTimeSeries <- function(se, type, geneID, symbol, condition, treatment, refTreat, addZero = FALSE, zeroTreat = NULL, timerange) {
  
  if (type == "expression") {
    # Handle zero time point addition if specified
    if (!is.null(zeroTreat) & addZero) {
      seqMat <- addZeroTime(se, condition, treatment, zeroTreat, timerange)
    }
    else {
      seqMat <- se[,se[[condition]] == treatment & se$timepoint %in% timerange]
    }
    yLabText <- "Normalized expression"
  } 
  else if (type == "logFC"){
    if (!is.null(zeroTreat) & addZero) {
      seSub <- se[, se[[condition]] == treatment]
      allTimepoint <- unique(seSub$timepoint)
      seRef <- se[, se[[condition]] == refTreat]
      timepointRef <- unique(seRef$timepoint)
      
      # Handle zero time point addition for both treatment and reference
      if (!("0min" %in% allTimepoint) & !("0min" %in% timepointRef)) {
        seqMat <- addZeroTime(se, condition, treatment, zeroTreat, timerange)
        refMat <- addZeroTime(se, condition, refTreat, zeroTreat, timerange)
      }
      else if (!("0min" %in% allTimepoint) & ("0min" %in% timepointRef)) {
        seqMat <- addZeroTime(se, condition, treatment, zeroTreat, timerange)
        refMat <- se[,se[[condition]] == refTreat & se$timepoint %in% timerange]
      }
      else if (("0min" %in% allTimepoint) & !("0min" %in% timepointRef)) {
        seqMat <- se[,se[[condition]] == treatment & se$timepoint %in% timerange]
        refMat <- addZeroTime(se, condition, refTreat, zeroTreat, timerange)
      }
    }
    else {
      seqMat <- se[,se[[condition]] == treatment & se$timepoint %in% timerange]
      refMat <- se[,se[[condition]] == refTreat & se$timepoint %in% timerange]
    }
    
    # Calculate log fold change
    # Here the mean intensities in refMat are calculated per time point or per time point and subject ID.
    if (!is.null(se$subjectID)) {
      fcMat <- lapply(unique(seqMat$timepoint), function(tp) {
        lapply(unique(seqMat$subjectID), function(id) {
          if (length(colnames(assay(refMat)[,refMat$timepoint == tp & refMat$subjectID == id])) > 1) {
            refMean = rowMeans(assay(refMat)[,refMat$timepoint == tp & 
                                               refMat$subjectID == id])
          } else {
            refMean = assay(refMat)[,refMat$timepoint == tp & refMat$subjectID == id]
          }
          assay(seqMat)[,seqMat$timepoint == tp & seqMat$subjectID == id] - refMean
        })
      }) %>% bind_cols() %>% as.matrix()
    } else {
      fcMat <- lapply(unique(seqMat$timepoint), function(tp) {
        if (length(colnames(assay(refMat)[,refMat$timepoint == tp])) > 1) {
          refMean = rowMeans(assay(refMat)[,refMat$timepoint == tp]) 
        } else {
          refMean = assay(refMat)[,refMat$timepoint == tp]
        }
        
        assay(seqMat)[,seqMat$timepoint == tp] - refMean
      }) %>% bind_cols() %>% as.matrix()
    } 
    
    rownames(fcMat) <- rownames(assay(seqMat))
    colnames(fcMat) <- colnames(assay(seqMat))
    # Rearrange columns in seqMat to match with fcMat to assign the fold change
    seqMat <- seqMat[,colnames(fcMat)]
    
    assay(seqMat) <- fcMat
    yLabText <- "logFC"
  } 
  else if (type == "two-condition expression") {
    if (!is.null(zeroTreat) & addZero) {
      seSub <- se[, se[[condition]] == treatment]
      allTimepoint <- unique(seSub$timepoint)
      seRef <- se[, se[[condition]] == refTreat]
      timepointRef <- unique(seRef$timepoint)
      
      # Handle zero time point addition for both conditions
      if (!("0min" %in% allTimepoint) & !("0min" %in% timepointRef)) {
        se1 <- addZeroTime(se, condition, treatment, zeroTreat, timerange)
        se2 <- addZeroTime(se, condition, refTreat, zeroTreat, timerange)
        assay <- cbind(assay(se1), assay(se2))
        cd <- rbind(colData(se1), colData(se2))
        emeta <- elementMetadata(se)
        
        seqMat <- SummarizedExperiment(assays=SimpleList(intensity=assay), colData = cd, rowData = emeta)
      }
      else if (!("0min" %in% allTimepoint) & ("0min" %in% timepointRef)) {
        se1 <- addZeroTime(se, condition, treatment, zeroTreat, timerange)
        se2 <- se[, se[[condition]] == refTreat & se$timepoint %in% timerange]
        assay <- cbind(assay(se1), assay(se2))
        cd <- rbind(colData(se1), colData(se2))
        emeta <- elementMetadata(se)
        
        seqMat <- SummarizedExperiment(assays=SimpleList(intensity=assay), colData = cd, rowData = emeta)
      }
      else if (("0min" %in% allTimepoint) & !("0min" %in% timepointRef)) {
        se1 <- se[, se[[condition]] == treatment & se$timepoint %in% timerange]
        se2 <- addZeroTime(se, condition, refTreat, zeroTreat, timerange)
        assay <- cbind(assay(se1), assay(se2))
        cd <- rbind(colData(se1), colData(se2))
        emeta <- elementMetadata(se)
        
        seqMat <- SummarizedExperiment(assays=SimpleList(intensity=assay), colData = cd, rowData = emeta)
      }
    }
    else {
      seqMat <- se[,se[[condition]] %in% c(treatment, refTreat) & se$timepoint %in% timerange]
    }
    yLabText <- "Normalized expression"
  }
  
  # Create data frame for plotting
  plotTab <- data.frame(time = seqMat$timepoint,
                        value = assay(seqMat)[geneID,], 
                        treatment = as.character(seqMat[[condition]])) 
  
  # Convert time to numerical values
  timeUnit <- str_extract(plotTab$time, "h|min")
  timeUnit <- ifelse(is.na(timeUnit), "", timeUnit)
  if ((any(timeUnit == "h")) & (any(timeUnit == "min"))) {
    plotTab[timeUnit == "min","time"] <- 1/60 * as.numeric(gsub("min", "", plotTab[timeUnit == "min","time"]))
  } 
  plotTab$time <- as.numeric(gsub("h|min", "", plotTab$time))
  
  # Create the ggplot object
  p <- ggplot(plotTab, aes(x= time, y = value)) +
    geom_point(aes(color = treatment), size=3) + 
    stat_summary(aes(color=paste("mean",treatment)),fun = mean, geom = "line", linewidth = 2)
  
  # Add subject-specific lines if subjectID is present
  if (!is.null(seqMat$subjectID)) {
    plotTab$subjectID <- seqMat$subjectID
    p <- ggplot(plotTab, aes(x= time, y = value,color = paste0(subjectID,"_",treatment))) +
      geom_point( size=3) + 
      stat_summary(fun = mean, geom = "line", linewidth = 1,linetype = "dashed")
    
  }
  
  # Finalize plot with labels and theme
  p <- p +
    ylab(yLabText) + xlab("time") + 
    ggtitle(symbol) + theme_bw() + 
    theme(text=element_text(size=15),plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5, size=15))
  return(p)
}