#' @name getDecouplerNetwork
#' 
#' @title Load Kinase-Substrate Interaction Network
#'
#' @description
#' `getDecouplerNetwork` loads the kinase-substrate interaction network for a specified species from pre-defined files.
#'
#' @param speciesRef A character string specifying the species. Supported values are "Homo sapiens" and "Mus musculus".
#'
#' @return A data frame containing the kinase-substrate interaction network for the specified species.
#'
#' @details
#' The function reads from tab-separated value (TSV) files located in the `omnipathR_kinase_network` directory.
#' It supports two species: Homo sapiens and Mus musculus.
#' 
#' @examples
#' # Load the human kinase-substrate interaction network
#' human_network <- getDecouplerNetwork("Homo sapiens")
#'
#' # Load the mouse kinase-substrate interaction network
#' mouse_network <- getDecouplerNetwork("Mus musculus")
#'
#' @importFrom utils read.table
#' @export
getDecouplerNetwork <- function(speciesRef) {
  
  # load network of kinase-substrate interaction from omnipathR_kinase_network folder
  if (speciesRef == "Homo sapiens") {
    decoupler_network <- data(read.table("omnipathR_kinase_network/Homo_sapiens.tsv", sep = "\t", stringsAsFactors = FALSE))
  } else if (speciesRef == "Mus musculus") {
    decoupler_network <- data(read.table("omnipathR_kinase_network/Mus_musculus.tsv", sep = "\t", stringsAsFactors = FALSE))
  }
  
  # Return the loaded network data
  return(decoupler_network)
}


#' @name calcKinaseScore
#' 
#' @title Calculate Kinase Activity Scores using `decoupleR`
#'
#' @description
#' `calcKinaseScore` calculates kinase activity scores based on input data and a specified network of regulatory relationships (decoupler network).
#'
#' @param resTab A data frame containing the input data with columns `site`, `stat`, and `log2FC`.
#' @param decoupler_network A data frame representing the decoupleR network with columns `source` and `target`.
#' @param corrThreshold A numeric value specifying the correlation threshold for filtering correlated regulons. Default is `0.9.
#' @param statType A character string specifying the type of statistic to use. Options are `"stat"` or `"log2FC"`. Default is `"stat"`.
#' @param nPerm Number of permutations for the null distribution. Default is `100`.
#'
#' @return A data frame with kinase activity scores, including columns for `source`, `score`, and `p_value`.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Removes duplicate rows based on the `site` column.
#'   \item Filters the data to include only those sites present in the `target` column of the `decoupler_network`.
#'   \item Prepares the input table based on the specified `statType`.
#'   \item Intersects the input table with the decoupler network to find common regulons.
#'   \item Checks for correlated regulons and filters out those exceeding the correlation threshold.
#'   \item Calculates kinase activity using a weighted mean approach.
#'   \item Processes the results to handle `NA` values and formats the output.
#' }
#'
#' @examples
#' # Example usage:
#' resTab <- data.frame(
#'   site = c("S1", "S2", "S3", "S4"),
#'   stat = c(1.5, -2.3, 0.7, 1.2),
#'   log2FC = c(0.5, -1.1, 0.3, 0.9)
#' )
#' decoupler_network <- data.frame(
#'   source = c("Kinase1", "Kinase2", "Kinase3"),
#'   target = c("S1", "S2", "S4")
#' )
#' result <- calcKinaseScore(resTab, decoupler_network, corrThreshold = 0.8, statType = "stat")
#' print(result)
#'
#' @importFrom dplyr distinct filter select rename mutate
#' @importFrom decoupleR intersect_regulons check_corr run_wmean
#' @importFrom tibble column_to_rownames
#' 
#' @export
calcKinaseScore <- function(resTab, decoupler_network, corrThreshold = 0.9, statType = "stat", nPerm = 100) {
  # Remove duplicate rows based on the 'site' column and keep all other columns
  resTab <- resTab %>%
    distinct(site, .keep_all = TRUE) %>%
    # Filter the rows where 'site' is present in the 'target' column of decoupler_network
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
  decoupler_network <- decoupleR::intersect_regulons(mat = inputTab, 
                                                     network = decoupler_network, 
                                                     .source = source, 
                                                     .target = target, 
                                                     minsize = 5)
  # Check for correlated regulons within the decoupler network and remove interactions with correlation >= threshold
  correlated_regulons <- decoupleR::check_corr(decoupler_network) %>%  #not necessary for now
    dplyr::filter(correlation >= corrThreshold)
  decoupler_network <- decoupler_network %>% 
    dplyr::filter(!source %in% correlated_regulons$source.2)
  # Calculate the kinase score by computing the weighted mean
  kinase_activity <- decoupleR::run_wmean(mat = as.matrix(inputTab), 
                                          network = decoupler_network,
                                          sparse = FALSE,
                                          times = nPerm)
  # Get the wmean statistics, replace NA scores with 0, and replace NA p_value with 1
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
#' `plotKinaseDE` generates a bar plot of the top kinases associated with the differentially expressed genes based on their scores.
#'
#' @param scoreTab A data frame containing kinase scores with columns `source`, `score`, and `p_value`.
#' @param nTop An integer specifying the number of top kinases to plot for each direction. Default is `10`.
#' @param pCut A numeric value specifying the p-value cutoff for significance. Default is `0.05`.
#'
#' @return A ggplot object representing the bar plot of kinase score.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Adds a column for significance based on the p-value cutoff.
#'   \item Adds a column for the sign of the score.
#'   \item Filters out kinases with a score of 0.
#'   \item Selects the top `nTop` kinases by absolute score for each sign of the score.
#'   \item Creates a bar plot with the selected kinases.
#' }
#'
#' @examples
#' # Example usage:
#' scoreTab <- data.frame(
#'   source = c("Kinase1", "Kinase2", "Kinase3", "Kinase4"),
#'   score = c(2.3, -1.5, 0, 3.1),
#'   p_value = c(0.01, 0.2, 0.05, 0.03)
#' )
#' plot <- plotKinaseDE(scoreTab, nTop = 3, pCut = 0.05)
#' print(plot)
#'
#' @importFrom dplyr mutate filter group_by slice_max
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual theme_linedraw theme element_text unit coord_flip ggtitle xlab ylab reorder
#' @export
plotKinaseDE <- function(scoreTab, nTop = 10, pCut = 0.05) {
  plotTab <- scoreTab %>% mutate(significance = ifelse(p_value <= pCut, paste0("p <= ",pCut), paste0("p > ",pCut)),
                                 score_sign = sign(score)) %>%
    # Remove kinases whose scores are 0 in the plot
    filter(score_sign != 0) %>%
    # Group by score sign and select the top nTop kinases by absolute score
    group_by(score_sign) %>% slice_max(abs(score), n = nTop)
  p <- ggplot(plotTab, aes(x = reorder(source, score), y = score)) + 
    geom_bar(aes(fill = significance), stat = "identity") 
  
  # Customize the fill color based on the significance levels
  if (length(unique(plotTab$significance)) == 2) {
    p <- p + scale_fill_manual(values = c("indianred", "lightgrey"), labels = c(paste0("p <= ",pCut), paste0("p > ",pCut)))
  } else if (unique(plotTab$significance) == paste0("p > ",pCut)) {
    p <- p + scale_fill_manual(values = "lightgrey", labels = paste0("p > ",pCut)) 
  } else if (unique(plotTab$significance) == c(paste0("p <= ",pCut))) {
    p <- p + scale_fill_manual(values = "indianred", labels = paste0("p <= ",pCut)) 
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
#' `plotKinaseTimeSeries` creates a heatmap to visualize the result of kinase activity inference for time-series clustering, with significant activity changes marked.
#'
#' @param scoreTab A data frame containing kinase activity scores, p-values, and time points.
#' @param pCut A numeric value specifying the p-value threshold for significance. Default is `0.05`.
#' @param clusterName A character string specifying the name of the cluster for the plot title. Default is `"cluster1"`.
#'
#' @return A ggplot2 object representing the heatmap of kinase activity score.
#'
#' @details
#' The heatmap shows kinase activity scores over different time points. Significant activities (based on the specified p-value threshold) are marked with an asterisk (*). The color gradient represents the activity score, with blue indicating low activity, red indicating high activity, and white as the midpoint.
#' 
#' @examples
#' # Example usage:
#' scoreTab <- data.frame(
#'   timepoint = rep(c("0h", "1h", "2h"), each = 3),
#'   source = rep(c("KinaseA", "KinaseB", "KinaseC"), times = 3),
#'   score = runif(9, -2, 2),
#'   p_value = runif(9, 0, 0.1)
#' )
#' p <- plotKinaseTimeSeries(scoreTab)
#' print(p)
#'
#' @importFrom dplyr mutate rename
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient2 scale_x_discrete scale_y_discrete theme_bw ylab xlab ggtitle theme element_text unit
#' @export
plotKinaseTimeSeries <- function(scoreTab, pCut = 0.05, clusterName = "cluster1") {
  
  # Add a significance marker based on the p-value threshold
  plotTab <- dplyr::mutate(scoreTab, sig = ifelse(p_value<=pCut, "*", ""))
  # Rename the score column for better readability in the plot
  plotTab <- plotTab %>% rename(Activity_score = "score")
  
  # Create the heatmap plot
  p <- ggplot(plotTab, aes(x=timepoint, y = source,fill = Activity_score)) +
    geom_tile() +
    geom_text(aes(label = sig), vjust = 0.5, hjust = 0.5, size = 10) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    scale_x_discrete(expand = c(0,0)) +scale_y_discrete(expand = c(0,0)) +
    theme_bw() +
    ylab("Kinase") + xlab("Time point") + ggtitle(paste("Kinase activity infererence,", clusterName)) + 
    theme(axis.title = element_text(face = "bold", size = 15),
          axis.text.y = element_text(size =15),
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
          plot.title = element_text(size = 17),
          legend.title = element_text(size =15, face= "bold"),
          legend.text = element_text(size =15),
          legend.key.height = unit(0.8, "cm"),
          legend.key.width = unit(0.8, "cm"))
  return(p)
}


#### Function to perform PTM-SEA, a modification of the GSEA algorithm to work on databases of site-centric ptm signatures, 
# in which the direction of regulation (up or down) are specified for each site.
# Most of this function was similar to the project.geneset() function from the GitHub page for 
# ssGSEA2.0: https://github.com/broadinstitute/ssGSEA2.0/tree/master. 
#
# The function was modified to be compatible with SmartPhosExplorer and make recognizing phosphosites
# in the signatures stricter. Specifically, a phosphosite is included in a signature if its sign of the 
# test statistics (+ / -) matches the direction of regulation (u / d, abbreviate for `up` and `down`) 
# in the signature. This was not considered in the original ssGSEA2.0 algorithm.
# The function in this script only consider phosphosites (`site.ptm == "p"`). Signatures starting with "KINASE"
# are not considered since they are targets of kinases and hence would be similar to the kinase activity analysis.
# --------------------------------------
# For the publication associated with the algorithm and database please see Krug et al., 2019, at 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6398202/ 
# For more descriptions of the signature sets please see the PTM signature database website:
# https://proteomics.broadapps.org/ptmsigdb/ 
# --------------------------------------
## Parameters
## geneStat: a dataframe with 1 column listing test statistics (t-stat or logFC) named "stat".
##           row.names are names of the phosphosites.
## ptmSetDb: database of PTM signature set.
## nPerm: Number of permutations
## weight: weighting of the sites based on their test statistic. If weight == 0 then the test statistics do not matter
## correl.type: Can be "rank", "z.score", or "symm.rank"
## statistic: Can be "Kolmogorov-Smirnov" or "area.under.RES"
## min.overlap: Minimum number of sites in the set to be considered for the analysis, default is 5.
# -----------------------------------------
## The resulting enrichment score is normalized to account for differences between signature set sizes


#' @name runGSEAforPhospho
#' 
#' @title Run GSEA for Phosphorylation Data
#'
#' @description 
#' `runGSEAforPhospho` performs Gene Set Enrichment Analysis (GSEA) for phosphorylation data.
#'
#' @param geneStat A data frame containing gene statistics, with gene names as row names and a column named 'stat' for the statistics.
#' @param ptmSetDb A data frame of post-translational modification (PTM) signature sets.
#' @param nPerm An integer specifying the number of permutations for the null distribution.
#' @param weight A numeric value for the weight parameter in the GSEA algorithm. If weight == 0 then the test statistics do not matter. Default is `1`.
#' @param correl.type A character string specifying the correlation type. Options are "rank", "symm.rank", and "z.score". Default is `"rank"`.
#' @param statistic A character string specifying the statistic to be used. Options are "Kolmogorov-Smirnov" and "area.under.RES". Default is `"Kolmogorov-Smirnov"`.
#' @param min.overlap An integer specifying the minimum overlap required between gene sets and the input data. Default is `5`.
#'
#' @return A tibble with enrichment scores and associated statistics for each PTM set.
#'
#' @details
#' This function runs GSEA on phosphorylation data to identify enriched PTM sets. It calculates enrichment scores and p-values for each set, normalizes the scores, and adjusts p-values for multiple testing.
#' 
#' @examples
#' # Example usage:
#' geneStat <- data.frame(stat = runif(100, -2, 2))
#' row.names(geneStat) <- paste0("Gene", 1:100)
#' ptmSetDb <- data.frame(signature = sample(letters, 100, replace = TRUE), category = "example", site.ptm = "p", site.direction = sample(c("u", "d"), 100, replace = TRUE))
#' result <- runGSEAforPhospho(geneStat, ptmSetDb, nPerm = 1000)
#' print(result)
#'
#' @importFrom dplyr mutate rename filter count tibble as_tibble group_by ungroup arrange bind_rows
#' @importFrom tidyr separate
#' @importFrom stats p.adjust
#' @export
runGSEAforPhospho <- function(geneStat, ptmSetDb, nPerm, weight = 1, correl.type = "rank",
                              statistic = "Kolmogorov-Smirnov", min.overlap = 5) {
  
  # Internal function to calculate GSEA enrichment score
  gseaScorePTM <- function (ordered.gene.list, data.expr, gene.set2, 
                            weight = 1, correl.type = "rank", gene.set.direction = NULL,
                            statistic = "Kolmogorov-Smirnov", min.overlap = 5) {
    
    # Function to calculate the enrichment score (ES)
    score <- function(max.ES, min.ES, RES, gaps, valleys, statistic){
      # KM
      if( statistic == "Kolmogorov-Smirnov" ){
        if( max.ES > -min.ES ){
          ES <- signif(max.ES, digits=5)
          arg.ES <- which.max(RES)
        } else{
          ES <- signif(min.ES, digits=5)
          arg.ES <- which.min(RES)
        }
      }
      # AUC
      if( statistic == "area.under.RES"){
        if( max.ES > -min.ES ){
          arg.ES <- which.max(RES)
        } else{
          arg.ES <- which.min(RES)
        }
        gaps = gaps+1
        RES = c(valleys,0) * (gaps) + 0.5*( c(0,RES) - c(valleys,0) ) * (gaps)
        ES = sum(RES)
      }
      return(list(RES=RES, ES=ES, arg.ES=arg.ES))
    } 
    
    n.rows = length(ordered.gene.list)
    
    # Apply weighting to the correlation vector
    if (weight == 0) {
      
      correl.vector <- rep(1, n.rows)
      
    } else if (weight > 0) {
      # If weighting is used (weight > 0), bring 'correl.vector' into the same order as the ordered gene list
      if (correl.type == "rank") {
        correl.vector <- data.expr[ordered.gene.list]
        
      } else if (correl.type == "symm.rank") {
        correl.vector <- data.expr[ordered.gene.list]
        
        correl.vector <- ifelse(correl.vector > correl.vector[ceiling(n.rows/2)],
                                correl.vector,
                                correl.vector + correl.vector - correl.vector[ceiling(n.rows/2)])
      } else if (correl.type == "z.score") {
        x <- data.expr[ordered.gene.list]
        correl.vector <- (x - mean(x))/sd(x)
      }
    }
    
    # Length of gene list is same as the number of rows in input matrix
    N = length(ordered.gene.list)
    
    
    # Sirectionality of the gene set
    if(!is.null(gene.set.direction)){
      
      # Number of 'd' features
      d.idx <- which(gene.set.direction=='d')
      Nh.d <- length(d.idx)
      Nm.d <-  N - Nh.d
      
      # Locations of 'd' features
      tag.d <- sign( match(ordered.gene.list, gene.set2[ d.idx ], nomatch=0) )
      if (weight == 0) {
        ind.d = which(tag.d == 1)} else {
          ind.d = which(tag.d == 1 & correl.vector < 0)}
      number.d = length(ind.d)
      
      # Number of 'u' features
      u.idx <- which(gene.set.direction=='u')
      Nh.u <- length(u.idx)
      Nm.u <-  N - Nh.u
      
      # Locations of 'up' features
      tag.u <- sign( match(ordered.gene.list, gene.set2[ u.idx ], nomatch=0) )
      if (weight == 0) {
        ind.u = which(tag.u == 1)} else {
          ind.u = which(tag.u == 1 & correl.vector >= 0)}
      number.u = length(ind.u)
      
      
      # For up-regulated genes/sites
      if(number.u > 1){
        
        # Extract and apply weighting
        correl.vector.u <- correl.vector[ind.u]
        correl.vector.u <- abs(correl.vector.u)^weight           ## weighting
        
        sum.correl.u <- sum(correl.vector.u)
        
        up.u <- correl.vector.u/sum.correl.u         ## steps up in th random walk
        gaps.u <- (c(ind.u-1, N) - c(0, ind.u))      ## gaps between hits
        down.u <- gaps.u/Nm.u                        ## steps down in the random walk
        
        RES.u <- cumsum(up.u-down.u[1:length(up.u)])  
        
        valleys.u = RES.u-up.u
        
        max.ES.u = suppressWarnings(max(RES.u))
        min.ES.u = suppressWarnings(min(valleys.u))
        
        # Calculate final score
        score.res <- score(max.ES.u, min.ES.u, RES.u, gaps.u, valleys.u, statistic)
        ES.u <- score.res$ES
        arg.ES.u <- score.res$arg.ES
        RES.u <- score.res$RES
        
      } else {
        correl.vector.u <- rep(0, N)
        ES.u=0
        RES.u=0
        arg.ES.u=NA
        up.u=0
        down.u=0
      }
      
      # For down-regulated genes/sites
      if(number.d > 1){  
        # Extract and apply weighting
        correl.vector.d <- correl.vector[ind.d]
        correl.vector.d <- abs(correl.vector.d)^weight           ## weighting
        
        sum.correl.d <- sum(correl.vector.d)
        
        up.d <- correl.vector.d/sum.correl.d
        gaps.d <- (c(ind.d-1, N) - c(0, ind.d))
        down.d <- gaps.d/Nm.d
        
        RES.d <- cumsum(up.d-down.d[1:length(up.d)])               
        valleys.d = RES.d-up.d
        
        max.ES.d = suppressWarnings(max(RES.d))
        min.ES.d = suppressWarnings(min(valleys.d))
        
        # Calculate final score
        score.res <- score(max.ES.d, min.ES.d, RES.d, gaps.d, valleys.d, statistic)
        ES.d <- score.res$ES
        arg.ES.d <- score.res$arg.ES
        RES.d <- score.res$RES
        
      } else {
        correl.vector.d <- rep(0, N)
        ES.d=0
        RES.d=0
        ind.d=NA
        number.d=0
        arg.ES.d=NA
        up.d=0
        down.d=0
      }
      
      # Make sure to meet the min.overlap threshold
      if(Nh.d == 1 & Nh.u < min.overlap | Nh.u == 1 & Nh.d < min.overlap){
        ES.u <- ES.d <- RES.u <- RES.d <- 0
        arg.ES <- arg.ES <- NA
        ind.u <- ind.d <- NULL
      }
      
      # Combine the results
      ES <- ES.u - ES.d
      RES <- list(u=RES.u, d=RES.d)
      arg.ES <- c(arg.ES.u, arg.ES.d)
      correl.vector = list(u=correl.vector.u, d=correl.vector.d)
      
      ind <- list(u=ind.u, d=ind.d)
      step.up <- list(u=up.u, d=up.d )
      step.down <- list(u=1/Nm.u, d=1/Nm.d)
      gsea.results = list(ES = ES, ES.all = list(u=ES.u, d=ES.d), arg.ES = arg.ES, RES = RES, indicator = ind, correl.vector = correl.vector, step.up=step.up, step.down=step.down,
                          number.u = number.u, number.d = number.d)
      
      
    } else { ## end  if(!is.null(gene.set.direction))
      
      # Without directionality
      Nh <- length(gene.set2)
      Nm <-  N - Nh
      
      
      # Match gene set to data
      tag.indicator <- sign(match(ordered.gene.list, gene.set2, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag)
      # Positions of gene set in ordered gene list
      ind = which(tag.indicator==1)
      # 'correl.vector' is now the size of 'gene.set2'
      correl.vector <- abs(correl.vector[ind])^weight
      # Sum of weights
      sum.correl = sum(correl.vector)
      
      # Determine peaks and valleys
      # Divide correl vector by sum of weights
      up = correl.vector/sum.correl     # "up" represents the peaks in the mountain plot
      gaps = (c(ind-1, N) - c(0, ind))  # gaps between ranked pathway genes
      down = gaps/Nm
      
      RES = cumsum(c(up,up[Nh])-down)
      valleys = RES[1:Nh]-up
      
      max.ES = max(RES)
      min.ES = min(valleys)
      
      # Calculate final score
      score.res <- score(max.ES, min.ES, RES[1:Nh], gaps, valleys, statistic)
      
      ES <- score.res$ES
      arg.ES <- score.res$arg.ES
      RES <- score.res$RES
      
      gsea.results = list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = ind, correl.vector = correl.vector, step.up=up, step.down=1/Nm)
    } 
    
    return (gsea.results)
  }
  
  
  # Remove KINASE signature since this is analogous to the kinase activity inference part
  ptmSetDbNoKinase <- ptmSetDb %>%
    filter(!grepl("KINASE", category))
  
  # Get the number of PTM sites for each signature
  ptmSiteCount <- ptmSetDbNoKinase %>%
    count(signature) %>%
    rename(no.PTM.site = "n")
  
  # Preprocessing the geneSetDatabase
  phosphoSetDb <- ptmSetDbNoKinase %>%
    dplyr::as_tibble() %>%
    filter(site.ptm == "p") %>%   
    group_by(signature) %>%
    filter(n() >= 5) %>%
    ungroup() %>%
    separate(site.annotation, sep =  ":", into = c("site", "PubMedID"), extra="merge", fill="right")
  
  # Get the number of phospho sites for each signature
  phosphoSiteCount <- phosphoSetDb %>%
    count(signature) %>%
    rename(no.phospho.site = "n")
  
  # Put input data in a format compatible with gseaScorePTM
  ordered.gene.list <- row.names(geneStat)
  data.expr <- geneStat$stat
  names(data.expr) <- ordered.gene.list
  # Run GSEA for each PTM set
  rtab <- lapply(phosphoSiteCount$signature, function(signature) {
    # Get number of PTM site and phospho site in the database
    nPTMsite = as.numeric(ptmSiteCount[ptmSiteCount$signature == signature, "no.PTM.site"])
    nPpSite = as.numeric(phosphoSiteCount[phosphoSiteCount$signature == signature, "no.phospho.site"])
    signatureSet = phosphoSetDb[phosphoSetDb$signature == signature,]
    gene.set2 = signatureSet$site
    gene.set.direction = signatureSet$site.direction
    gene.set.PMID = signatureSet$PubMedID
    # Calculate the gsea score
    if (sum(row.names(geneStat) %in% gene.set2) < min.overlap) {
      enrichScoreNorm <- enrichScore <- pvalue <- number.u <- number.d <- 0
    } else {
      resGSEA = gseaScorePTM(ordered.gene.list, data.expr =  data.expr, gene.set2 = gene.set2, 
                             weight = weight, correl.type = correl.type,
                             gene.set.direction = gene.set.direction, min.overlap =  min.overlap)
      enrichScore = resGSEA$ES
      if (!is.null(gene.set.direction)) {
        number.u  = resGSEA$number.u
        number.d = resGSEA$number.d
      } else {
        number.u <- number.d <- 0 
      }
      # Calculate the null distribution and pvalue
      if (nPerm == 0) {
        enrichScoreNorm = enrichScore
        pvalue = 1
      } else {
        nullDistES = sapply(1:nPerm,  function(x) gseaScorePTM(sample(ordered.gene.list), data.expr=data.expr, gene.set2=gene.set2, 
                                                               weight, correl.type, gene.set.direction = gene.set.direction, min.overlap = min.overlap)$ES)
        nullDistES = unlist(nullDistES)
        if (enrichScore >= 0) {
          nullDistES.pos = nullDistES[nullDistES >= 0]
          if (length(nullDistES.pos) == 0) nullDistES.pos = 0.5
          posMean = mean(nullDistES.pos)
          enrichScoreNorm = enrichScore/posMean
          s = sum(nullDistES.pos >= enrichScore)/length(nullDistES.pos)
          pvalue = ifelse(s == 0, 1/nPerm, s)
        } else {
          nullDistES.neg = nullDistES[nullDistES < 0]
          if (length(nullDistES.neg) == 0) nullDistES.neg = 0.5
          negMean = mean(nullDistES.neg)
          enrichScoreNorm = enrichScore/negMean
          s = sum(nullDistES.neg <= enrichScore)/length(nullDistES.neg)
          pvalue = ifelse(s == 0, 1/nPerm, s)
        }
      }
    }
    tibble(Name = signature,
           nSite = number.u + number.d,                      # number of phosphosites in the input data
           enrichScore = enrichScoreNorm,                    # normalized enrichment score to correct for differences in signature sizes
           n.P.site.in.Db = nPpSite,                         # number of phosphosites in the database
           n.PTM.site.in.Db = nPTMsite,                      # number of PTM sites in the database
           pvalue = pvalue,
           number.u = number.u,
           number.d = number.d)
  }
  ) %>% bind_rows() %>%
    filter(nSite>= min.overlap, enrichScore!=0) %>%
    mutate(p.adj = p.adjust(pvalue, method = "BH")) %>%
    arrange(pvalue) 
  return(rtab)
}