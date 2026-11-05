# The `RangedSummarizedExperiment7` class

An S7 implementation of the biocspec `RangedSummarizedExperiment`
interface: a `SummarizedExperiment7` whose per-feature identity and
metadata are anchored to `rowRanges` rather than an independent
`rowData`/`rowNames` pair (clause `RSE-INV-3`).

## Usage

``` r
RangedSummarizedExperiment7(
  assays = NULL,
  rowRanges,
  colData = NULL,
  metadata = list()
)
```

## Arguments

- assays:

  A matrix, list/named list of matrices, or an `Assays7` object. `NULL`
  for a zero-assay container.

- rowRanges:

  Any object with a [`length()`](https://rdrr.io/r/base/length.html)
  method, one element per feature (row); typically a
  [`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  or
  [`GenomicRanges::GRangesList`](https://rdrr.io/pkg/GenomicRanges/man/GRangesList-class.html),
  but duck-typed so a plain list or vector also works for the purposes
  of this demonstration package.

- colData:

  A `data.frame` of sample-level annotations, one row per sample;
  `rownames(colData)` become `colnames(x)`.

- metadata:

  A named list of experiment-wide metadata.

## Examples

``` r
m <- matrix(1:6, nrow = 2, ncol = 3, dimnames = list(c("g1", "g2"), c("s1", "s2", "s3")))
rse <- RangedSummarizedExperiment7(
  assays = m,
  rowRanges = list(g1 = 1:10, g2 = 11:20)
)
rse
#> <RangedSummarizedExperiment7> 2 x 3
#> assays:  
#> rowRanges class: list 
#> colnames: s1, s2, s3 
dim(rse)
#> [1] 2 3
rowRanges7(rse)
#> $g1
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
#> $g2
#>  [1] 11 12 13 14 15 16 17 18 19 20
#> 
```
