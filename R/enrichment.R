#' @name runFisher
#' 
#' @title Perform Fisher's Exact Test on Gene Sets
#'
#' @description
#' `runFisher` performs Fisher's Exact Test to determine the enrichment of a set of genes within reference gene sets.
#'
#' @param genes A character vector of genes of interest.
#' @param reference A character vector of reference genes.
#' @param inputSet A list containing gene set collections. If `ptm` is TRUE, this should be a data frame with specific columns.
#' @param ptm Logical. If TRUE, perform the test on post-translational modification (PTM) gene sets. Default is `FALSE`.
#'
#' @return A data frame with the results of the Fisher's Exact Test, including the gene set name, the number of genes in the set, set size, p-value, adjusted p-value, and the genes in the set.
#'
#' @details
#' The function can operate in two modes: standard gene sets and PTM-specific gene sets. For PTM-specific gene sets, additional filtering and processing are performed.
#'
#' @examples
#' # Example usage:
#' # Assuming `genesOfInterest` is a vector of genes and `referenceGenes` is a vector of reference genes
#' # and `geneSetCollection` is a list of gene sets
#' # results <- runFisher(genesOfInterest, referenceGenes, geneSetCollection, ptm = FALSE)
#'
#' @importFrom dplyr filter group_by ungroup mutate bind_rows arrange tibble
#' @importFrom tidyr separate
#' @importFrom stats fisher.test p.adjust
#' @export
runFisher <- function (genes, reference, inputSet, ptm = FALSE) {
  
  # Retrieve the gene sets
  if (!ptm) {
    genesets <- inputSet$gsc
    setList <- 1:length(genesets)
  } else {
    # Filter and process the PTM-specific gene sets
    genesets <- inputSet %>%
      filter(!grepl("KINASE", category)) %>%
      dplyr::as_tibble() %>%
      filter(site.ptm == "p") %>%
      group_by(signature) %>%
      filter(n() >= 5) %>%
      ungroup() %>%
      mutate(signature = ifelse(site.direction == "u", paste0(signature,"_upregulated"), paste0(signature, "_downregulated"))) %>%
      separate(site.annotation, sep =  ":", into = c("site", "PubMedID"), extra="merge", fill="right") %>%
      as.data.frame()
    
    if(nrow(genesets) == 0) stop("Genesets empty. Check inputSet argument.")
    setList <- unique(genesets$signature)
  }
  # Remove genes from the reference set that are in the genes of interest
  reference = reference[!reference %in% genes]
  
  # Apply Fisher's Exact Test to each gene set
  rtab = lapply(setList, function(i) { 
    # Identify the geneset and its name
    if (!ptm) {
      geneset = genesets[[i]]
      nameSet = names(genesets)[i]
    } else {
      geneset = genesets[genesets$signature == i, "site"]
      nameSet = i
    }
    
    # Create the contingency table for Fisher's Exact Test
    RinSet = sum(reference %in% geneset)
    RninSet = length(reference) - RinSet
    GinSet = sum(genes %in% geneset)
    GninSet = length(genes) - GinSet
    fmat = matrix(c(GinSet, RinSet, GninSet, RninSet), nrow = 2,
                  ncol = 2, byrow = F)
    colnames(fmat) = c("inSet", "ninSet")
    rownames(fmat) = c("genes", "reference")
    
    # Perform Fisher's Exact Test
    fish = fisher.test(fmat, alternative = "greater")
    pval = fish$p.value
    inSet = RinSet + GinSet
    
    # Return the result as a tibble
    tibble(Name = nameSet,
           `Gene.number`= GinSet, 
           `Set.size` = inSet, 
           pval = pval,
           Genes = list(intersect(genes, geneset)))
  }) %>% bind_rows() %>%
    filter(Set.size>0) %>%
    mutate(padj = p.adjust(pval, method = "BH")) %>%
    arrange(pval)
  
  return(data.frame(rtab))
}


#' @name clusterEnrich
#' 
#' @title Perform Cluster Enrichment Analysis
#'
#' @description
#' `clusterEnrich` performs enrichment analysis on gene clusters, using Fisher's Exact Test to determine the significance of enrichment for each cluster.
#'
#' @param clusterTab A data frame containing cluster information, where each row corresponds to a gene and its assigned cluster.
#' @param se A SummarizedExperiment object containing gene expression data and metadata.
#' @param inputSet A list or data frame of gene sets to be used for enrichment analysis.
#' @param reference A character vector of reference genes. If NULL, it will be extracted from `se`. Default is `NULL`.
#' @param ptm Logical. If TRUE, the function will perform enrichment analysis on post-translational modification (PTM) gene sets. Default is `FALSE`.
#' @param adj Character. The method for adjusting p-values. Default is `"BH"`.
#' @param filterP Numeric. The p-value threshold for filtering significant results. Default is `0.05`.
#' @param ifFDR Logical. If TRUE, the function will use FDR-adjusted p-values for significance filtering. Default is `FALSE`.
#'
#' @return A list containing two elements:
#' \itemize{
#'   \item `table`: A data frame with enrichment results for each cluster and pathway.
#'   \item `plot`: A ggplot object showing the significance of enrichment for each pathway across clusters.
#' }
#'
#' @details
#' The function first retrieves or computes the reference set of genes or PTM sites. It then performs enrichment analysis for each cluster using the `runFisher` function.
#' The results are filtered based on the p-value threshold and adjusted for multiple testing if `ifFDR` is `TRUE`. The function generates a dot plot where the size and color of the points represent the significance of enrichment.
#'
#' @examples
#' # Example usage:
#' # Assuming `clusterTable` is a data frame with cluster assignments,
#' # `summarizedExp` is a SummarizedExperiment object,
#' # and `geneSetCollection` is a list or data frame of gene sets
#' # results <- clusterEnrich(clusterTable, summarizedExp, geneSetCollection, ptm = FALSE)
#'
#' @importFrom dplyr filter mutate group_by summarise ungroup bind_rows arrange
#' @importFrom ggplot2 ggplot geom_point aes scale_fill_gradient xlab ylab theme element_text
#' @importFrom tidyr separate
#' @importFrom stats p.adjust
#' @export
clusterEnrich <- function(clusterTab, se, inputSet, reference = NULL, ptm = FALSE, adj = "BH", filterP = 0.05, ifFDR = FALSE) {
  
  # If reference is not provided, derive it from the SummarizedExperiment object
  if (is.null(reference)) {
    if (ptm) {
      reference <- rowData(se)$site
    }
    else {
      reference <- unique(rowData(se)$Gene)
    }
  } 
  
  # Perform Fisher's Exact Test for each unique cluster
  resTabFisher <- lapply(unique(clusterTab$cluster), function(cc) {
    # Extract gene IDs for the current cluster
    id <- filter(clusterTab, cluster == cc)$feature
    if (ptm) {
      genes <- unique(elementMetadata(se)[id,]$site)
    }
    else {
      genes <- unique(elementMetadata(se)[id,]$Gene)
    }
    
    # Run Fisher's Exact Test and annotate results with the cluster ID
    eachOut <- runFisher(genes, reference, inputSet, ptm) %>%
      mutate(cluster = cc)
  }) %>% bind_rows()
  
  # Filter results based on significance and prepare for plotting
  if (ifFDR) {
    plotTab <- resTabFisher %>%
      filter(padj <= filterP) %>% arrange(padj) %>%
      mutate(ifSig = padj <= filterP) %>%
      group_by(Name) %>% mutate(atLeast1 = sum(ifSig)>0) %>%
      filter(atLeast1) %>% ungroup()
  }
  else {
    plotTab <- resTabFisher %>%
      filter(pval <= filterP) %>% arrange(pval) %>%
      mutate(ifSig = pval <= filterP) %>%
      group_by(Name) %>% mutate(atLeast1 = sum(ifSig)>0) %>%
      filter(atLeast1) %>% ungroup()
  }
  
  # Create a ggplot object for visualization of enrichment results
  p<- ggplot(plotTab, aes(x=cluster, y=Name, customdata = cluster, key = Name)) +
    geom_point(aes(size =-log10(pval),fill=-log10(pval)), shape = 21, color = "black") +
    scale_fill_gradient(low = "white", high = "red") +
    xlab("Cluster") +
    ylab("Pathway") +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  return(list(table = plotTab, plot = p))
}