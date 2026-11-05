# Get or set the list of assays

Get or set the list of assays

## Usage

``` r
assays7(x)

assays7(x) <- value
```

## Arguments

- x:

  A `SummarizedExperiment7`.

- value:

  Replacement assays: a matrix, list/named list of matrices, or an
  `Assays7` object.

## Examples

``` r
se <- SummarizedExperiment7(assays = list(counts = matrix(1:6, 2, 3)))
assays7(se)
#> $counts
#>      [,1] [,2] [,3]
#> [1,]    1    3    5
#> [2,]    2    4    6
#> 
assays7(se) <- list(counts = matrix(1:6, 2, 3), logcounts = matrix(7:12, 2, 3))
assayNames7(se)
#> [1] "counts"    "logcounts"
```
