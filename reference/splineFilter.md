# Filter Expression Matrix Using Spline Models

`splineFilter` filters an expression `matrix` based on spline models
fitted to time-series data, optionally considering treatment and subject
ID.

## Usage

``` r
splineFilter(
  exprMat,
  subjectID = NULL,
  time,
  df,
  pCut = 0.5,
  ifFDR = FALSE,
  treatment = NULL,
  refTreatment = NULL
)
```

## Arguments

- exprMat:

  A `numeric` matrix of expression data, where rows are features and
  columns are samples.

- subjectID:

  `Character`. An optional vector of subject IDs corresponding to
  columns in exprMat. Default is `NULL`.

- time:

  A `numeric` vector representing the time points corresponding to
  columns in exprMat.

- df:

  A `numeric` value specifying the degrees of freedom for the spline
  basis.

- pCut:

  A `numeric` value for the p-value cutoff to filter significant
  features. Default is 0.05.

- ifFDR:

  A `logical` value indicating if the false discovery rate (FDR) should
  be used for filtering. If `FALSE`, raw p-values are used. Default is
  `FALSE`.

- treatment:

  `Character`. An optional vector of treatment labels corresponding to
  columns in exprMat. Default is `NULL`.

- refTreatment:

  `Character`. An optional reference treatment label for the treatment
  vector. Default is `NULL`.

## Value

A filtered expression `matrix` containing only the features that meet
the significance criteria.

## Details

The function performs the following steps:

1.  Converts time points from minutes to hours if both units are
    present.

2.  Removes rows with missing values from the expression `matrix`.

3.  Constructs a design matrix for the spline model, optionally
    including subject IDs and treatments.

4.  Fits a linear model using the design matrix and performs empirical
    Bayes moderation.

5.  Extracts significant features based on the specified p-value or FDR
    cutoff.
