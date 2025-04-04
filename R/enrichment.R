#' @name runFisher
#'
#' @title Perform Fisher's Exact Test on Gene Sets
#'
#' @description
#' \code{runFisher} performs Fisher's Exact Test to determine the enrichment of
#' a set of genes within reference gene sets.
#'
#' @param genes A \code{character} vector of genes of interest.
#' @param reference A \code{character} vector of reference genes.
#' @param inputSet A \code{list} containing gene set collections. If ptm is
#' \code{TRUE}, this should be a data frame with specific columns.
#' @param ptm \code{Logical}. If \code{TRUE}, perform the test on
#' post-translational modification (PTM) gene sets. Default is \code{FALSE}.
#'
#' @return A \code{data frame} with the results of the Fisher's Exact Test,
#' including the gene set name, the number of genes in the set, set size,
#' p-value, adjusted p-value, and the genes in the set.
#'
#' @details
#' The function can operate in two modes: standard gene sets and PTM-specific
#' gene sets. For PTM-specific gene sets, additional filtering and processing
#' are performed.
#'
#' @examples
#' library(SummarizedExperiment)
#' library(piano)
#' # Load multiAssayExperiment object
#' data("dda_example")
#' # Get SummarizedExperiment object
#' se <- dda_example[["Proteome"]]
#' colData(se) <- colData(dda_example)
#' # Preprocess the proteome assay
#' result <- preprocessProteome(se, normalize = TRUE)
#' # Call the function to perform differential expression analyis
#' de <- performDifferentialExp(se = result, assay = "Intensity",
#' method = "limma", reference = "1stCrtl", target = "EGF",
#' condition = "treatment")
#' genesList <- unique(de$resDE$Gene)
#' referenceList <- unique(SummarizedExperiment::rowData(result)$Gene)
#' genesetPath <- appDir <- system.file("shiny-app/geneset",
#' package = "SmartPhos")
#' inGMT <- loadGSC(paste0(genesetPath,"/Cancer_Hallmark.gmt"),type="gmt")
#' # Run the function
#' runFisher(genes = genesList, reference = referenceList, inputSet = inGMT)
#'
#' @importFrom dplyr filter group_by ungroup mutate bind_rows arrange tibble n
#' @importFrom tidyr separate
#' @importFrom stats fisher.test p.adjust
#'
#' @export
runFisher <- function (genes, reference, inputSet, ptm = FALSE) {

  category <- Set.size <- pval <- NULL
  site.ptm <- site.direction <- site.annotation <- NULL

  # Retrieve the gene sets
  if (!ptm) {
    genesets <- inputSet$gsc
    setList <- 1 : length(genesets)
  } else {
    # Filter and process the PTM-specific gene sets
    genesets <- inputSet %>%
      filter(!grepl("KINASE", category)) %>%
      dplyr::as_tibble() %>%
      filter(site.ptm == "p") %>%
      group_by(signature) %>%
      filter(n() >= 5) %>%
      ungroup() %>%
      mutate(signature = ifelse(site.direction == "u",
                                paste0(signature,"_upregulated"),
                                paste0(signature, "_downregulated"))) %>%
      separate(site.annotation, sep =  ":", into = c("site", "PubMedID"),
               extra="merge", fill="right") %>%
      as.data.frame()

    if(nrow(genesets) == 0) stop("Genesets empty. Check inputSet argument.")
    setList <- unique(genesets$signature)
  }
  # Remove genes from the reference set that are in the genes of interest
  reference <-  reference[!reference %in% genes]

  # Apply Fisher's Exact Test to each gene set
  rtab <-  lapply(setList, function(i) {
    # Identify the geneset and its name
    if (!ptm) {
      geneset <-  genesets[[i]]
      nameSet <-  names(genesets)[i]
    } else {
      geneset <-  genesets[genesets$signature == i, "site"]
      nameSet <-  i
    }

    # Create the contingency table for Fisher's Exact Test
    RinSet <-  sum(reference %in% geneset)
    RninSet <-  length(reference) - RinSet
    GinSet <-  sum(genes %in% geneset)
    GninSet <-  length(genes) - GinSet
    fmat <-  matrix(c(GinSet, RinSet, GninSet, RninSet), nrow = 2,
                  ncol = 2, byrow = FALSE)
    colnames(fmat) <-  c("inSet", "ninSet")
    rownames(fmat) <-  c("genes", "reference")

    # Perform Fisher's Exact Test
    fish <-  fisher.test(fmat, alternative = "greater")
    pval <-  fish$p.value
    inSet <-  RinSet + GinSet

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


#' @name enrichDifferential
#'
#' @title Perform Enrichment analysis on differentially expressed genes
#' or phospho-sites
#'
#' @description
#' \code{enrichDifferential} performs enrichment analysis on differentially
#' expressed genes and phospho-sites for either pathway or phospho-specific
#' enrichment, depending on the input parameters. It supports multiple
#' statistical methods such as PAGE and GSEA for pathway enrichment and a
#' Kolmogorov-Smirnov approach for phospho-enrichment.
#'
#' @param dea A \code{data frame} containing the differential expression
#' analysis results. It should include columns like
#'   `pvalue`, `Gene` (or `site`), `stat`, and `log2FC`.
#' @param type A \code{character} string indicating the type of enrichment.
#' Options are `"Pathway enrichment"` or `"Phospho enrichment"`.
#' @param gsaMethod A \code{character} string specifying the gene set analysis
#' method for pathway enrichment. Options are `"PAGE"` or `"GSEA"`.
#' @param geneSet A gene set collection to use for pathway enrichment.
#' @param ptmSet A post-translational modification (PTM) set database for
#' phospho-enrichment analysis.
#' @param statType A \code{character} string specifying the statistic type to
#' use. Options are `"stat"` or `"log2FC"`.
#' @param nPerm A \code{numeric} specifying the number of permutations for GSEA.
#' Default is 100.
#' @param sigLevel A \code{numeric} value representing the significance
#' threshold for filtering results. Ddefault is 0.05.
#' @param ifFDR A \code{logical} value indicating whether to filter results
#' using FDR-adjusted p-values. Default is `FALSE`.
#'
#' @return A data frame containing the results of the enrichment analysis,
#' including columns such as the gene set name,
#'   statistical significance, and adjusted p-values.
#'
#' @details
#' The `enrichDifferential` function performs either pathway enrichment or
#' phospho-enrichment analysis based on the `type` parameter.
#' For pathway enrichment, it uses either the PAGE or GSEA method with a
#' provided gene set collection. For phospho-enrichment,
#' it uses a Kolmogorov-Smirnov test with a PTM set database. Results can be
#' filtered by significance level and optionally adjusted for FDR.
#'
#' @examples
#' library(SummarizedExperiment)
#' library(piano)
#' # Load multiAssayExperiment object
#' data("dia_example")
#' # Get SummarizedExperiment object
#' se <- dia_example[["Phosphoproteome"]]
#' colData(se) <- colData(dia_example)
#' # Preprocess the proteome assay
#' result <- preprocessPhos(se, normalize = TRUE)
#' # Call the function to perform differential expression analyis
#' dea <- performDifferentialExp(se = result, assay = "Intensity",
#' method = "limma", reference = "1stCrtl", target = "EGF",
#' condition = "treatment")
#' # Load the gene set
#' genesetPath <- appDir <- system.file("shiny-app/geneset",
#' package = "SmartPhos")
#' inGMT <- loadGSC(paste0(genesetPath,"/Cancer_Hallmark.gmt"),type="gmt")
#' # Call the function
#' resTab <- enrichDifferential(dea = dea$resDE, type = "Pathway enrichment",
#' gsaMethod = "PAGE", geneSet = inGMT, statType = "stat", nPerm = 200,
#' sigLevel = 0.05, ifFDR = FALSE)
#'
#' @importFrom dplyr filter arrange select
#' @importFrom piano runGSA GSAsummaryTable
#'
#' @export
enrichDifferential <- function(dea, type = c("Pathway enrichment",
                                             "Phospho enrichment"),
                               gsaMethod = c("PAGE", "GSEA"),
                               geneSet, ptmSet,
                               statType = c("stat", "log2FC"),
                               nPerm = 100, sigLevel = 0.05, ifFDR = FALSE) {

    type <- match.arg(type)
    gsaMethod <- match.arg(gsaMethod)
    statType <- match.arg(statType)

    pvalue <- Gene <- p.up.adj <- p.down.adj <- desc <- p.up <- p.down <- NULL
    site <- Name <- Site.number <- Number.up <- Number.down <- NULL
    Number.pSite.Db <- Number.PTM.site.Db <- padj <- NULL

    if (type == "Pathway enrichment") {
        # Prepare data for pathway enrichment by filtering and sorting
        corTab <- dea %>%
            arrange(pvalue) %>%
            filter(!duplicated(Gene)) %>% # Remove duplicates by gene
            arrange(stat)

        # Select appropriate statistic column for analysis
        if(statType == "stat") {
            myCoef <- data.frame(row.names = corTab$Gene,
                                 stat = corTab$stat,
                                 stringsAsFactors = FALSE)
        }
        else {
            myCoef <- data.frame(row.names = corTab$Gene,
                                 stat = corTab$log2FC,
                                 stringsAsFactors = FALSE)
        }

        # Perform pathway enrichment analysis using the chosen method
        if (gsaMethod == "PAGE") {
            res <- runGSA(geneLevelStats = myCoef,
                          geneSetStat = "page",
                          adjMethod = "fdr",
                          gsc = geneSet,
                          signifMethod = 'nullDist')
        }
        else if (gsaMethod == "GSEA") {
            res <- runGSA(geneLevelStats = myCoef,
                          geneSetStat = "gsea",
                          adjMethod = "fdr",
                          gsc = geneSet,
                          signifMethod = 'geneSampling',
                          nPerm = nPerm)
        }

        # Summarize the results into a table
        resTab <- GSAsummaryTable(res)
        colnames(resTab) <- c("Name", "Gene Number", "Stat", "p.up",
                              "p.up.adj", "p.down", "p.down.adj", "Number up",
                              "Number down")

        # Filter results based on significance threshold
        if(ifFDR) {
            resTab <- filter(resTab,
                             p.up.adj <= sigLevel | p.down.adj <= sigLevel) %>%
                arrange(desc(Stat))
        } else {
            resTab <- filter(resTab,
                             p.up <= sigLevel | p.down <= sigLevel) %>%
                arrange(desc(Stat))
        }
    }
    else {
        # Prepare data for phospho-enrichment by filtering and sorting
        corTab <- dea %>%
            arrange(pvalue) %>%
            filter(!duplicated(site)) %>% # Remove duplicates by site
            arrange(desc(stat))

        # Select appropriate statistic column for analysis
        if (statType == "stat") {
            myCoef <- data.frame(row.names = corTab$site,
                                 stat = corTab$stat,
                                 stringsAsFactors = FALSE)
        } else if (statType == "log2FC") {
            myCoef <- data.frame(row.names = corTab$site,
                                 stat = corTab$log2FC,
                                 stringsAsFactors = FALSE)
        }

        # Perform phospho-enrichment analysis using the Kolmogorov-Smirnov test
        resTab <- runGSEAforPhospho(geneStat = myCoef, ptmSetDb = ptmSet,
                                    nPerm =  nPerm,
                                    weight = 1, correl.type = "rank",
                                    statistic = "Kolmogorov-Smirnov",
                                    min.overlap = 5) %>%
            as.data.frame()

        # Rename columns for clarity
        colnames(resTab) <- c("Name", "Site.number", "Stat", "Number.pSite.Db",
                              "Number.PTM.site.Db", "pvalue", "Number.up",
                              "Number.down", "padj")

        # Rearrange columns in a specific order
        resTab <- resTab %>%
            select(Name,Site.number,Stat,Number.up,Number.down,Number.pSite.Db,
                   Number.PTM.site.Db, pvalue, padj) # rearrange column order

        # Filter results based on significance threshold
        if (ifFDR) {
            resTab <- filter(resTab, padj <= sigLevel) %>%
                arrange(desc(Stat))
        } else {
            resTab <- filter(resTab, pvalue <= sigLevel) %>%
                arrange(desc(Stat))
        }
    }
    return(resTab)
}


#' @name clusterEnrich
#'
#' @title Perform Cluster Enrichment Analysis
#'
#' @description
#' \code{clusterEnrich} performs enrichment analysis on gene clusters, using
#' Fisher's Exact Test to determine the significance of enrichment for each
#' cluster.
#'
#' @param clusterTab A \code{data frame} containing cluster information, where
#' each row corresponds to a gene and its assigned cluster.
#' @param se A \code{SummarizedExperiment} object containing gene expression
#' data and metadata.
#' @param inputSet A \code{list} or \code{data frame} of gene sets to be used
#' for enrichment analysis.
#' @param reference A \code{character} vector of reference genes. If NULL, it
#' will be extracted from se object. Default is \code{NULL}.
#' @param ptm \code{Logical}. If \code{TRUE}, the function will perform
#' enrichment analysis on post-translational modification (PTM) gene sets.
#' Default is \code{FALSE}.
#' @param adj \code{Character}. The method for adjusting p-values. Default is
#' "BH".
#' @param filterP \code{Numeric}. The p-value threshold for filtering
#' significant results. Default is 0.05.
#' @param ifFDR \code{Logical}. If \code{TRUE}, the function will use
#' FDR-adjusted p-values for significance filtering. Default is \code{FALSE}.
#'
#' @return A \code{list} containing two elements:
#' \itemize{
#'   \item `table`: A \code{data frame} with enrichment results for each
#'   cluster and pathway.
#'   \item `plot`: A \code{ggplot2} object showing the significance of
#'   enrichment for each pathway across clusters.
#' }
#'
#' @details
#' The function first retrieves or computes the reference set of genes or PTM
#' sites. It then performs enrichment analysis for each cluster using the
#' \code{runFisher} function.
#' The results are filtered based on the p-value threshold and adjusted for
#' multiple testing if \code{ifFDR} is \code{TRUE}. The function generates a
#' dot plot where the size and color of the points represent the significance
#' of enrichment.
#'
#' @examples
#' library(SummarizedExperiment)
#' # Load multiAssayExperiment object
#' data("dia_example")
#' # Get SummarizedExperiment object
#' se <- dia_example[["Phosphoproteome"]]
#' colData(se) <- colData(dia_example)
#' seProcess <- preprocessPhos(seData = se, normalize = TRUE, impute = "QRILC")
#' result <- addZeroTime(seProcess, condition = "treatment", treat = "EGF",
#' zeroTreat = "1stCrtl", timeRange = c("20min","40min", "6h"))
#' # Get the numeric matrix
#' exprMat <- SummarizedExperiment::assay(result)
#' # Call the clustering function
#' clust <- clusterTS(x = exprMat, k = 3)
#' genesetPath <- appDir <- system.file("shiny-app/geneset",
#' package = "SmartPhos")
#' inGMT <- piano::loadGSC(paste0(genesetPath,
#' "/Cancer_Hallmark.gmt"),type="gmt")
#' # Call the function
#' clusterEnrich(clust$cluster, seProcess, inGMT)
#'
#' @importFrom dplyr filter mutate group_by summarise ungroup bind_rows
#' arrange n
#' @importFrom tidyr separate
#' @importFrom stats p.adjust
#' @import ggplot2
#'
#' @export
clusterEnrich <- function(clusterTab, se, inputSet, reference = NULL,
                          ptm = FALSE, adj = "BH", filterP = 0.05,
                          ifFDR = FALSE) {

  cluster <- padj <- ifSig <- pval <- Name <- atLeast1 <- NULL
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
  p <- ggplot(plotTab,
              aes(x=cluster, y=Name, customdata = cluster, key = Name)) +
    geom_point(aes(size =-log10(pval),fill=-log10(pval)), shape = 21,
               color = "black") +
    scale_fill_gradient(low = "white", high = "red") +
    xlab("Cluster") +
    ylab("Pathway") +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  return(list(table = plotTab, plot = p))
}

#' @name runGSEAforPhospho
#'
#' @title Run GSEA for Phosphorylation Data
#'
#' @description
#' \code{runGSEAforPhospho} performs Gene Set Enrichment Analysis (GSEA) for
#' phosphorylation data.
#'
#' @param geneStat A \code{data frame} containing gene statistics, with gene
#' names as row names and a column named 'stat' for the statistics.
#' @param ptmSetDb A \code{data frame} of post-translational modification (PTM)
#' signature sets.
#' @param nPerm A \code{numeric} value specifying the number of permutations for
#' the null distribution.
#' @param weight A \code{numeric} value for the weight parameter in the GSEA
#' algorithm. If weight == 0 then the test statistics do not matter. Default
#' is 1.
#' @param correl.type A \code{character} string specifying the correlation type.
#' Options are "rank", "symm.rank", and "z.score". Default is "rank".
#' @param statistic A \code{character} string specifying the statistic to be
#' used. Options are "Kolmogorov-Smirnov" and "area.under.RES". Default is
#' "Kolmogorov-Smirnov".
#' @param min.overlap A \code{numeric} specifying the minimum overlap required
#' between gene sets and the input data. Default is 5.
#'
#' @return A \code{tibble} with enrichment scores and associated statistics for
#' each PTM set.
#'
#' @details
#' This function runs GSEA on phosphorylation data to identify enriched PTM
#' sets. It calculates enrichment scores and p-values for each set, normalizes
#' the scores, and adjusts p-values for multiple testing.
#'
#' @importFrom dplyr mutate rename filter count tibble as_tibble group_by
#' ungroup arrange bind_rows
#' @importFrom tidyr separate
#' @importFrom stats p.adjust
#'
#' @export
runGSEAforPhospho <- function(geneStat, ptmSetDb, nPerm, weight = 1,
                              correl.type = c("rank", "symm.rank", "z.score"),
                              statistic = c("Kolmogorov-Smirnov",
                                            "area.under.RES"),
                              min.overlap = 5) {

  correl.type <- match.arg(correl.type)
  statistic <- match.arg(statistic)
  category <- site.ptm <- site.annotation <- nSite <- NULL
  enrichScore <- pvalue <- NULL

  # Internal function to calculate GSEA enrichment score
  gseaScorePTM <- function (ordered.gene.list, data.expr, gene.set2,
                            weight = 1, correl.type = "rank",
                            gene.set.direction = NULL,
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
        gaps <- gaps+1
        RES <-  c(valleys,0) * (gaps) + 0.5*( c(0,RES) - c(valleys,0) ) * (gaps)
        ES <-  sum(RES)
      }
      return(list(RES=RES, ES=ES, arg.ES=arg.ES))
    }

    n.rows <-  length(ordered.gene.list)

    # Apply weighting to the correlation vector
    if (weight == 0) {

      correl.vector <- rep(1, n.rows)

    } else if (weight > 0) {
      # If weighting is used (weight > 0), bring 'correl.vector' into
      # the same order as the ordered gene list
      if (correl.type == "rank") {
        correl.vector <- data.expr[ordered.gene.list]

      } else if (correl.type == "symm.rank") {
        correl.vector <- data.expr[ordered.gene.list]

        correl.vector <- ifelse(correl.vector >
                                    correl.vector[ceiling(n.rows/2)],
                                correl.vector, correl.vector +
                                    correl.vector -
                                    correl.vector[ceiling(n.rows/2)])
      } else if (correl.type == "z.score") {
        x <- data.expr[ordered.gene.list]
        correl.vector <- (x - mean(x))/sd(x)
      }
    }

    # Length of gene list is same as the number of rows in input matrix
    N <-  length(ordered.gene.list)


    # Sirectionality of the gene set
    if(!is.null(gene.set.direction)){

      # Number of 'd' features
      d.idx <- which(gene.set.direction=='d')
      Nh.d <- length(d.idx)
      Nm.d <-  N - Nh.d

      # Locations of 'd' features
      tag.d <- sign( match(ordered.gene.list, gene.set2[ d.idx ], nomatch=0) )
      if (weight == 0) {
        ind.d <-  which(tag.d == 1)} else {
          ind.d <-  which(tag.d == 1 & correl.vector < 0)}
      number.d <-  length(ind.d)

      # Number of 'u' features
      u.idx <- which(gene.set.direction=='u')
      Nh.u <- length(u.idx)
      Nm.u <-  N - Nh.u

      # Locations of 'up' features
      tag.u <- sign( match(ordered.gene.list, gene.set2[ u.idx ], nomatch=0) )
      if (weight == 0) {
        ind.u <-  which(tag.u == 1)} else {
          ind.u <-  which(tag.u == 1 & correl.vector >= 0)}
      number.u <-  length(ind.u)


      # For up-regulated genes/sites
      if(number.u > 1){

        # Extract and apply weighting
        correl.vector.u <- correl.vector[ind.u]
        correl.vector.u <- abs(correl.vector.u)^weight       ## weighting

        sum.correl.u <- sum(correl.vector.u)

        up.u <- correl.vector.u/sum.correl.u     ## steps up in th random walk
        gaps.u <- (c(ind.u-1, N) - c(0, ind.u))  ## gaps between hits
        down.u <- gaps.u/Nm.u                    ## steps down in random walk

        RES.u <- cumsum(up.u-down.u[1:length(up.u)])

        valleys.u <-  RES.u-up.u

        max.ES.u <-  suppressWarnings(max(RES.u))
        min.ES.u <-  suppressWarnings(min(valleys.u))

        # Calculate final score
        score.res <- score(max.ES.u, min.ES.u, RES.u, gaps.u,
                           valleys.u, statistic)
        ES.u <- score.res$ES
        arg.ES.u <- score.res$arg.ES
        RES.u <- score.res$RES

      } else {
        correl.vector.u <- rep(0, N)
        ES.u <- 0
        RES.u <- 0
        arg.ES.u <- NA
        up.u <- 0
        down.u <- 0
      }

      # For down-regulated genes/sites
      if(number.d > 1){
        # Extract and apply weighting
        correl.vector.d <- correl.vector[ind.d]
        correl.vector.d <- abs(correl.vector.d)^weight     ## weighting

        sum.correl.d <- sum(correl.vector.d)

        up.d <- correl.vector.d/sum.correl.d
        gaps.d <- (c(ind.d-1, N) - c(0, ind.d))
        down.d <- gaps.d/Nm.d

        RES.d <- cumsum(up.d-down.d[1:length(up.d)])
        valleys.d <-  RES.d-up.d

        max.ES.d <-  suppressWarnings(max(RES.d))
        min.ES.d <-  suppressWarnings(min(valleys.d))

        # Calculate final score
        score.res <- score(max.ES.d, min.ES.d, RES.d, gaps.d, valleys.d,
                           statistic)
        ES.d <- score.res$ES
        arg.ES.d <- score.res$arg.ES
        RES.d <- score.res$RES

      } else {
        correl.vector.d <- rep(0, N)
        ES.d <- 0
        RES.d <- 0
        ind.d <- NA
        number.d <- 0
        arg.ES.d <- NA
        up.d <- 0
        down.d <- 0
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
      gsea.results <-  list(ES = ES, ES.all = list(u=ES.u, d=ES.d),
                            arg.ES = arg.ES, RES = RES, indicator = ind,
                            correl.vector = correl.vector, step.up=step.up,
                            step.down=step.down, number.u = number.u,
                            number.d = number.d)


    } else { ## end  if(!is.null(gene.set.direction))

      # Without directionality
      Nh <- length(gene.set2)
      Nm <-  N - Nh


      # Match gene set to data
      # notice that the sign is 0 (no tag) or 1 (tag)
      tag.indicator <- sign(match(ordered.gene.list, gene.set2, nomatch=0))
      # Positions of gene set in ordered gene list
      ind <-  which(tag.indicator==1)
      # 'correl.vector' is now the size of 'gene.set2'
      correl.vector <- abs(correl.vector[ind])^weight
      # Sum of weights
      sum.correl <-  sum(correl.vector)

      # Determine peaks and valleys
      # Divide correl vector by sum of weights
      # "up" represents the peaks in the mountain plot
      up <-  correl.vector/sum.correl
      # gaps between ranked pathway genes
      gaps <-  (c(ind-1, N) - c(0, ind))
      down <-  gaps/Nm

      RES <-  cumsum(c(up,up[Nh])-down)
      valleys <-  RES[1:Nh]-up

      max.ES <-  max(RES)
      min.ES <-  min(valleys)

      # Calculate final score
      score.res <- score(max.ES, min.ES, RES[1:Nh], gaps, valleys, statistic)

      ES <- score.res$ES
      arg.ES <- score.res$arg.ES
      RES <- score.res$RES

      gsea.results <- list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = ind,
                           correl.vector = correl.vector, step.up = up,
                           step.down = 1/Nm)
    }

    return (gsea.results)
  }


  # Remove KINASE signature since this is analogous to the kinase
  # activity inference part
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
    separate(site.annotation, sep =  ":", into = c("site", "PubMedID"),
             extra="merge", fill="right")

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
    nPTMsite <-  as.numeric(ptmSiteCount[ptmSiteCount$signature == signature,
                                         "no.PTM.site"])
    nPpSite <-  as.numeric(phosphoSiteCount[phosphoSiteCount$signature ==
                                                signature, "no.phospho.site"])
    signatureSet <-  phosphoSetDb[phosphoSetDb$signature == signature,]
    gene.set2 <-  signatureSet$site
    gene.set.direction <-  signatureSet$site.direction
    gene.set.PMID <-  signatureSet$PubMedID
    # Calculate the gsea score
    if (sum(row.names(geneStat) %in% gene.set2) < min.overlap) {
      enrichScoreNorm <- enrichScore <- pvalue <- number.u <- number.d <- 0
    } else {
      resGSEA <-  gseaScorePTM(ordered.gene.list, data.expr =  data.expr,
                               gene.set2 = gene.set2, weight = weight,
                               correl.type = correl.type,
                               gene.set.direction = gene.set.direction,
                               min.overlap =  min.overlap)
      enrichScore <-  resGSEA$ES
      if (!is.null(gene.set.direction)) {
        number.u  <-  resGSEA$number.u
        number.d <-  resGSEA$number.d
      } else {
        number.u <- number.d <- 0
      }
      # Calculate the null distribution and pvalue
      if (nPerm == 0) {
        enrichScoreNorm <-  enrichScore
        pvalue = 1
      } else {
        nullDistES <-  sapply(seq_len(nPerm),
                              function(x) gseaScorePTM(
                                  sample(ordered.gene.list),
                                  data.expr=data.expr,
                                  gene.set2=gene.set2,
                                  weight,
                                  correl.type,
                                  gene.set.direction = gene.set.direction,
                                  min.overlap = min.overlap)$ES)
        nullDistES <-  unlist(nullDistES)
        if (enrichScore >= 0) {
          nullDistES.pos <-  nullDistES[nullDistES >= 0]
          if (length(nullDistES.pos) == 0) nullDistES.pos <-  0.5
          posMean <-  mean(nullDistES.pos)
          enrichScoreNorm <-  enrichScore/posMean
          s <-  sum(nullDistES.pos >= enrichScore)/length(nullDistES.pos)
          pvalue <-  ifelse(s == 0, 1/nPerm, s)
        } else {
          nullDistES.neg <-  nullDistES[nullDistES < 0]
          if (length(nullDistES.neg) == 0) nullDistES.neg <-  0.5
          negMean <-  mean(nullDistES.neg)
          enrichScoreNorm <-  enrichScore/negMean
          s <-  sum(nullDistES.neg <= enrichScore)/length(nullDistES.neg)
          pvalue <-  ifelse(s == 0, 1/nPerm, s)
        }
      }
    }
    tibble(Name = signature,
           nSite = number.u + number.d,
           enrichScore = enrichScoreNorm, # normalized enrichment score
           n.P.site.in.Db = nPpSite,      # number of phosphosites in database
           n.PTM.site.in.Db = nPTMsite,   # number of PTM sites in database
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
