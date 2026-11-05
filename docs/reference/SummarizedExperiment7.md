# The `SummarizedExperiment7` class

An S7 implementation of the biocspec `SummarizedExperiment` interface: a
container coordinating one or more matrix-like assays with row
("feature") and column ("sample") metadata tables.

## Usage

``` r
SummarizedExperiment7(
  assays = NULL,
  rowData = NULL,
  colData = NULL,
  metadata = list(),
  rowNames = NULL
)
```

## Arguments

- assays:

  A matrix, list/named list of matrices, or an `Assays7` object. `NULL`
  for a zero-assay container.

- rowData:

  A `data.frame` of feature-level annotations, one row per feature.

- colData:

  A `data.frame` of sample-level annotations, one row per sample;
  `rownames(colData)` become `colnames(x)`.

- metadata:

  A named list of experiment-wide metadata.

- rowNames:

  Optional character vector of feature identifiers (`[[NAMES]]` in the
  specification); independent of `rowData`.

## Examples

``` r
m <- matrix(1:6, nrow = 2, ncol = 3, dimnames = list(c("g1", "g2"), c("s1", "s2", "s3")))
se <- SummarizedExperiment7(
  assays = list(counts = m),
  colData = data.frame(group = c("a", "a", "b"), row.names = c("s1", "s2", "s3"))
)
se
#> <SummarizedExperiment7> 2 x 3
#> assays: counts 
#> rownames: g1, g2 
#> colnames: s1, s2, s3 
dim(se)
#> [1] 2 3
rownames(se)
#> [1] "g1" "g2"
colData7(se)
#>    group
#> s1     a
#> s2     a
#> s3     b
```
