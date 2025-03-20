#' @name getDecouplerNetwork
#'
#' @title Load Kinase-Substrate Interaction Network
#'
#' @description
#' \code{getDecouplerNetwork} loads the kinase-substrate interaction network for
#'  a specified species from pre-defined files.
#'
#' @param speciesRef A \code{character} string specifying the species. Supported
#' values are "Homo sapiens" and "Mus musculus". Default is "Homo sapiens".
#'
#' @return A \code{data frame} containing the kinase-substrate interaction
#' network for the specified species.
#'
#' @examples
#' # Load the human kinase-substrate interaction network
#' getDecouplerNetwork("Homo sapiens")
#'
#' # Load the mouse kinase-substrate interaction network
#' getDecouplerNetwork("Mus musculus")
#'
#' @importFrom utils read.table data
#'
#' @export
getDecouplerNetwork <- function(speciesRef = "Homo sapiens") {

  # load network of kinase-substrate interaction from
  # omnipathR_kinase_network folder
  if (speciesRef == "Homo sapiens") {
    data("Homo_sapien_kinase_substrate_network")
    return(Homo_sapien_kinase_substrate_network)
  } else if (speciesRef == "Mus musculus") {
    data("Mus_musculus_kinase_substrate_network")
    return(Mus_musculus_kinase_substrate_network)
  }
}

#' @name calcKinaseScore
#'
#' @title Calculate Kinase Activity Scores using \code{decoupleR}
#'
#' @description
#' \code{calcKinaseScore} calculates kinase activity scores based on input data
#' and a specified network of regulatory relationships
#' (\code{decoupler network}).
#'
#' @param resTab A \code{data frame} containing the input data with columns
#' site, stat, and log2FC.
#' @param decoupler_network A \code{data frame} representing the decoupleR
#' network with columns source and target.
#' @param corrThreshold A \code{numeric} value specifying the correlation
#' threshold for filtering correlated regulons. Default is 0.9.
#' @param statType A \code{character} string specifying the type of statistic
#' to use. Options are "stat" or "log2FC". Default is "stat".
#' @param nPerm A \code{numeric} value specifying the number of permutations
#' for the null distribution. Default is 100.
#'
#' @return A \code{data frame} with kinase activity scores, including columns
#' for `source`, `score`, and `p_value`.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Removes duplicate rows based on the site column.
#'   \item Filters the data to include only those sites present in the target
#'   column of the \code{decoupler network}.
#'   \item Prepares the input table based on the specified statType.
#'   \item Intersects the input table with the \code{decoupler network} to find
#'   common regulons.
#'   \item Checks for correlated regulons and filters out those exceeding the
#'   correlation threshold.
#'   \item Calculates kinase activity using a weighted mean approach.
#'   \item Processes the results to handle NA values and formats the output.
#' }
#'
#' @importFrom dplyr distinct filter select rename mutate
#' @importFrom decoupleR intersect_regulons check_corr run_wmean
#' @importFrom tibble column_to_rownames
#'
#' @examples
#' resTab <- data.frame(
#' site = c("EGFR_Y1172", "EGFR_Y1197", "EGFR_S1166", "ROCK2_S1374",
#' "WASL_Y256", "GAB1_Y259", "ADD1_S586", "EPHA2_Y772", "PRKDC_T2638",
#' "PRKDC_T2609", "PRKDC_S2612"),
#' stat = c(-10.038770, -5.945562, 5.773384, -7.303834, 5.585326, 5.971104,
#' 5.199119, -5.169500, 5.130228, 5.407387, 4.493933),
#' log2FC = c(-2.6113343, -2.4858615, 1.0056629, -1.1561780, 1.6421145,
#' 2.0296634, 1.3766283, -0.8531656, 1.0742881, 1.0042942, 1.0608129)
#' )
#' decoupler_network <- data.frame(
#' source = c(rep("ABL1", 5), rep("CDK2", 6)),
#' mor = c(rep(1, 11)),
#' target = c("EGFR_Y1172", "EGFR_Y1197", "EGFR_S1166", "ROCK2_S1374",
#' "WASL_Y256", "GAB1_Y259", "ADD1_S586", "EPHA2_Y772", "PRKDC_T2638",
#' "PRKDC_T2609", "PRKDC_S2612"),
#' likelihood = c(rep(1, 11))
#' )
#' # Call the function
#' calcKinaseScore(resTab, decoupler_network)
#'
#'
#' @export
calcKinaseScore <- function(resTab, decoupler_network, corrThreshold = 0.9,
                            statType = "stat", nPerm = 100) {

  # Remove duplicate rows based on the 'site' column and keep all other columns
  resTab <- resTab %>%
    distinct(site, .keep_all = TRUE) %>%
    # Filter the rows where 'site' is present in the 'target' column of
    # decoupler_network
    filter(site %in% decoupler_network$target)

  # Prepare the input table based on the specified statType
  if (statType == "stat") {
    inputTab <- resTab %>% select(site, stat) %>% dplyr::rename(t = stat)
  } else if (statType == "log2FC") {
    inputTab <- resTab %>% select(site, log2FC)
  }
  rownames(inputTab) <- NULL
  inputTab <- inputTab %>% data.frame() %>% column_to_rownames("site")
  # Intersect the input table with the decoupler network to find common regulons
  decoupler_network <- decoupleR::intersect_regulons(
      mat = inputTab, network = decoupler_network, .source = source,
      .target = target, minsize = 5)
  # Check for correlated regulons within the decoupler network and
  # remove interactions with correlation >= threshold
  #not necessary for now
  correlated_regulons <- decoupleR::check_corr(decoupler_network) %>%
    dplyr::filter(correlation >= corrThreshold)
  decoupler_network <- decoupler_network %>%
    dplyr::filter(!source %in% correlated_regulons$source.2)
  # Calculate the kinase score by computing the weighted mean
  kinase_activity <- decoupleR::run_wmean(mat = as.matrix(inputTab),
                                          network = decoupler_network,
                                          sparse = FALSE,
                                          times = nPerm)
  # Get the wmean statistics, replace NA scores with 0, and replace
  # NA p_value with 1
  kinase_activity <- kinase_activity %>% dplyr::filter(statistic == "wmean") %>%
    select(-statistic, -condition) %>%
    mutate(score = ifelse(is.na(score), 0, score),
           p_value = ifelse(is.na(p_value), 1, p_value))
  return(kinase_activity)
}


#' @name plotKinaseDE
#'
#' @title Plot Kinase score for Differential Expression data
#'
#' @description
#' `plotKinaseDE` generates a \code{bar plot} of the top kinases associated
#' with the differentially expressed genes based on their scores.
#'
#' @param scoreTab A \code{data frame} containing kinase scores with columns
#' source, score, and p_value.
#' @param nTop A \code{numeric} value specifying the number of top kinases to
#' plot for each direction. Default is 10.
#' @param pCut A \code{numeric} value specifying the p-value cutoff for
#' significance. Default is 0.05.
#'
#' @return A \code{ggplot2} object representing the bar plot of kinase score.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Adds a column for significance based on the p-value cutoff.
#'   \item Adds a column for the sign of the score.
#'   \item Filters out kinases with a score of 0.
#'   \item Selects the top \code{nTop} kinases by absolute score for each sign
#'   of the score.
#'   \item Creates a bar plot with the selected kinases.
#' }
#'
#' @examples
#' # Example usage:
#' scoreTab <- data.frame(
#'  source = c("Kinase1", "Kinase2", "Kinase3", "Kinase4"),
#'  score = c(2.3, -1.5, 0, 3.1),
#'  p_value = c(0.01, 0.2, 0.05, 0.03)
#' )
#' plotKinaseDE(scoreTab, nTop = 3, pCut = 0.05)
#'
#' @importFrom dplyr mutate filter group_by slice_max
#' @importFrom stats reorder
#' @import ggplot2
#'
#' @export
plotKinaseDE <- function(scoreTab, nTop = 10, pCut = 0.05) {
  plotTab <- scoreTab %>% mutate(significance = ifelse(p_value <= pCut,
                                                       paste0("p <= ",pCut),
                                                       paste0("p > ",pCut)),
                                 score_sign = sign(score)) %>%
    # Remove kinases whose scores are 0 in the plot
    filter(score_sign != 0) %>%
    # Group by score sign and select the top nTop kinases by absolute score
    group_by(score_sign) %>% slice_max(abs(score), n = nTop)
  p <- ggplot(plotTab, aes(x = reorder(source, score), y = score)) +
    geom_bar(aes(fill = significance), stat = "identity")

  # Customize the fill color based on the significance levels
  if (length(unique(plotTab$significance)) == 2) {
    p <- p + scale_fill_manual(values = c("indianred", "lightgrey"),
                               labels = c(paste0("p <= ",pCut),
                                          paste0("p > ",pCut)))
  } else if (unique(plotTab$significance) == paste0("p > ",pCut)) {
    p <- p + scale_fill_manual(values = "lightgrey",
                               labels = paste0("p > ",pCut))
  } else if (unique(plotTab$significance) == c(paste0("p <= ",pCut))) {
    p <- p + scale_fill_manual(values = "indianred",
                               labels = paste0("p <= ",pCut))
  }
  p <- p +
    theme_linedraw() +
    theme(axis.title = element_text(face = "bold", size = 15),
          axis.text.y = element_text(size =15),
          axis.text.x = element_text(size = 15),
          plot.title = element_text(size = 17),
          legend.title = element_text(size =15, face= "bold"),
          legend.text = element_text(size =15),
          legend.key.height = unit(0.8, "cm"),
          legend.key.width = unit(0.8, "cm"),
          #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab("Kinases") + ylab("Kinase score") +
    ggtitle(paste0("Kinase activity inference, top ", nTop, " kinases")) +
    coord_flip()
  return(p)
}


#' @name plotKinaseTimeSeries
#'
#' @title Plot Kinase Activity Time Series
#'
#' @description
#' \code{plotKinaseTimeSeries} creates a heatmap to visualize the result of
#' kinase activity inference for time-series clustering, with significant
#' activity changes marked.
#'
#' @param scoreTab A \code{data frame} containing kinase activity scores,
#' p-values, and time points.
#' @param pCut A \code{numeric} value specifying the p-value threshold for
#' significance. Default is 0.05.
#' @param clusterName A \code{character} string specifying the name of the
#' cluster for the plot title. Default is "cluster1".
#'
#' @return A \code{ggplot2} object representing the heatmap of kinase activity
#' score.
#'
#' @details
#' The heatmap shows kinase activity scores over different time points.
#' Significant activities (based on the specified p-value threshold) are
#' marked with an asterisk (*). The color gradient represents the activity
#' score, with blue indicating low activity, red indicating high activity, and
#' white as the midpoint.
#'
#' @examples
#' # Example usage:
#' scoreTab <- data.frame(
#'  timepoint = rep(c("0h", "1h", "2h"), each = 3),
#'  source = rep(c("KinaseA", "KinaseB", "KinaseC"), times = 3),
#'  score = runif(9, -2, 2),
#'  p_value = runif(9, 0, 0.1)
#' )
#' plotKinaseTimeSeries(scoreTab)
#'
#' @importFrom dplyr mutate rename
#' @import ggplot2
#' @export
plotKinaseTimeSeries <- function(scoreTab, pCut = 0.05,
                                 clusterName = "cluster1") {

  # Add a significance marker based on the p-value threshold
  plotTab <- dplyr::mutate(scoreTab, sig = ifelse(p_value<=pCut, "*", ""))
  # Rename the score column for better readability in the plot
  plotTab <- plotTab %>% rename(Activity_score = "score")

  # Create the heatmap plot
  p <- ggplot(plotTab, aes(x=timepoint, y = source,fill = Activity_score)) +
    geom_tile() +
    geom_text(aes(label = sig), vjust = 0.5, hjust = 0.5, size = 10) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_bw() +
    ylab("Kinase") +
    xlab("Time point") +
    ggtitle(paste("Kinase activity infererence,", clusterName)) +
    theme(axis.title = element_text(face = "bold", size = 15),
          axis.text.y = element_text(size =15),
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1,
                                     hjust = 1),
          plot.title = element_text(size = 17),
          legend.title = element_text(size =15, face= "bold"),
          legend.text = element_text(size =15),
          legend.key.height = unit(0.8, "cm"),
          legend.key.width = unit(0.8, "cm"))
  return(p)
}
