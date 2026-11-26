# Perform Fisher's Exact Test on Gene Sets

`runFisher` performs Fisher's Exact Test to determine the enrichment of
a set of genes within reference gene sets.

## Usage

``` r
runFisher(genes, reference, inputSet, ptm = FALSE)
```

## Arguments

- genes:

  A `character` vector of genes of interest.

- reference:

  A `character` vector of reference genes.

- inputSet:

  A `list` containing gene set collections. If ptm is `TRUE`, this
  should be a data frame with specific columns.

- ptm:

  `Logical`. If `TRUE`, perform the test on post-translational
  modification (PTM) gene sets. Default is `FALSE`.

## Value

A `data frame` with the results of the Fisher's Exact Test, including
the gene set name, the number of genes in the set, set size, p-value,
adjusted p-value, and the genes in the set.

## Details

The function can operate in two modes: standard gene sets and
PTM-specific gene sets. For PTM-specific gene sets, additional filtering
and processing are performed.

## Examples

``` r
library(SummarizedExperiment)
library(piano)
# Load multiAssayExperiment object
data("dda_example")
# Get SummarizedExperiment object
se <- dda_example[["Proteome"]]
colData(se) <- colData(dda_example)
# Preprocess the proteome assay
result <- preprocessProteome(se, normalize = TRUE)
# Call the function to perform differential expression analyis
de <- performDifferentialExp(se = result, assay = "Intensity",
method = "limma", reference = "1stCrtl", target = "EGF",
condition = "treatment")
genesList <- unique(de$resDE$Gene)
referenceList <- unique(SummarizedExperiment::rowData(result)$Gene)
genesetPath <- appDir <- system.file("shiny-app/geneset",
package = "SmartPhos")
inGMT <- loadGSC(paste0(genesetPath,"/Cancer_Hallmark.gmt"),type="gmt")
# Run the function
runFisher(genes = genesList, reference = referenceList, inputSet = inGMT)
#>                                          Name Gene.number Set.size pval
#> 1            HALLMARK_TNFA_SIGNALING_VIA_NFKB           1        1    1
#> 2                            HALLMARK_HYPOXIA           3        3    1
#> 3                    HALLMARK_MITOTIC_SPINDLE           3        3    1
#> 4            HALLMARK_IL6_JAK_STAT3_SIGNALING           1        1    1
#> 5                         HALLMARK_DNA_REPAIR           2        2    1
#> 6                     HALLMARK_G2M_CHECKPOINT           1        1    1
#> 7                          HALLMARK_APOPTOSIS           2        2    1
#> 8            HALLMARK_ESTROGEN_RESPONSE_EARLY           5        5    1
#> 9             HALLMARK_ESTROGEN_RESPONSE_LATE           2        2    1
#> 10                 HALLMARK_ANDROGEN_RESPONSE           2        2    1
#> 11                        HALLMARK_MYOGENESIS           3        3    1
#> 12                 HALLMARK_PROTEIN_SECRETION           2        2    1
#> 13                   HALLMARK_APICAL_JUNCTION           1        1    1
#> 14                        HALLMARK_COMPLEMENT           2        2    1
#> 15         HALLMARK_UNFOLDED_PROTEIN_RESPONSE           1        1    1
#> 16                  HALLMARK_MTORC1_SIGNALING           4        4    1
#> 17                       HALLMARK_E2F_TARGETS           2        2    1
#> 18                    HALLMARK_MYC_TARGETS_V1           2        2    1
#> 19                    HALLMARK_MYC_TARGETS_V2           2        2    1
#> 20 HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION           2        2    1
#> 21             HALLMARK_XENOBIOTIC_METABOLISM           1        1    1
#> 22         HALLMARK_OXIDATIVE_PHOSPHORYLATION           2        2    1
#> 23                        HALLMARK_GLYCOLYSIS           4        4    1
#> 24   HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY           2        2    1
#> 25                    HALLMARK_UV_RESPONSE_DN           1        1    1
#> 26                   HALLMARK_HEME_METABOLISM           2        2    1
#> 27               HALLMARK_IL2_STAT5_SIGNALING           1        1    1
#> 28              HALLMARK_BILE_ACID_METABOLISM           1        1    1
#> 29                        HALLMARK_PEROXISOME           2        2    1
#> 30                 HALLMARK_KRAS_SIGNALING_DN           2        2    1
#>           Genes padj
#> 1         CYR61    1
#> 2  CYR61, I....    1
#> 3  ARHGAP10....    1
#> 4          HAX1    1
#> 5  SDCBP, S....    1
#> 6        DDX39A    1
#> 7   DFFA, DNM1L    1
#> 8  KRT19, K....    1
#> 9  KRT19, S....    1
#> 10  KRT19, KRT8    1
#> 11 BIN1, MY....    1
#> 12 AP3B1, DNM1L    1
#> 13         WASL    1
#> 14 ATOX1, KIF2A    1
#> 15        NOP56    1
#> 16 PLOD2, P....    1
#> 17 NOP56, D....    1
#> 18 NOP56, P....    1
#> 19  PES1, NOP56    1
#> 20 CYR61, PLOD2    1
#> 21       PGRMC1    1
#> 22 PDHX, NDUFS8    1
#> 23 PLOD2, K....    1
#> 24 PDLIM1, ....    1
#> 25        CYR61    1
#> 26 SDCBP, PSMD9    1
#> 27        MYO1C    1
#> 28      SULT2B1    1
#> 29 NUDT19, ....    1
#> 30  KRT15, KRT5    1
```
